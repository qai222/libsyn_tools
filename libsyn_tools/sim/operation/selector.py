"""
Dynamic participant-selection helpers for `libsyn_tools.sim`.
The goal is to allow the user to specify action participants in a dynamic, non-deterministic way.
E.g., "use a clean vial that is available" or "use a transfer device that has not been used for a particular liquid"
We model the availability in this selection using simpy.FilterStore (FIFO).
User-defined selection criteria are passed to the FilterStore as a `predicate` (single argument lambda).

The selection, based on the given criteria, can be done at different levels, for each level we have a selector class:
Level 1 `LiteralSelector`: exact IRI
Level 2 `AttributeSelector`: chooses the first LabObject in a pool whose *current* in-memory attributes satisfy the
given predicate.
Level 3 `HistorySelector`: predicate may also inspect the LabObject’s recent history.
Level 4 `KgQuerySelector`: predicate is a generic query to the full KG, not implemented yet

All selectors guarantee **atomic selection + locking**: as soon as a
LabObject is picked it is locked (via a per-object `simpy.Resource`)
before the IRI is returned to the caller, preventing race conditions.

Known limitation: multi-resource selection from small pools can deadlock
under contention even with stable ordering (e.g., each operation holds one
resource and waits for another). This is not prevented by the current
implementation. Mitigation: avoid operations that require N resources from
pool size N when concurrent; serialize those operations; or increase pool
size.
"""

from __future__ import annotations

import inspect
from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import Generator, Tuple

import simpy
from loguru import logger
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph import LabObject
from libsyn_tools.sim.knowledge_graph import identifier_from_iri
from libsyn_tools.sim.env_utils import get_effect_engine
from libsyn_tools.sim.operation.runtime import get_runtime_context, get_runtime_state
from libsyn_tools.sim.validation import require_singleton_or_error


class FilterStoreRegistry:
    @classmethod
    def get_filter_store(cls, pool_type: str, env: simpy.Environment | None = None) -> simpy.FilterStore:
        if env is None:
            raise RuntimeError("env required to access a filter store")
        ctx = get_runtime_context(env)
        if pool_type not in ctx.filter_stores:
            ctx.filter_stores[pool_type] = simpy.FilterStore(env)
        return ctx.filter_stores[pool_type]

    @classmethod
    def put_obj_into_filter_store(cls, obj: LabObject, env: simpy.Environment):
        if not obj.has_pool_type:
            return
        pool_type = require_singleton_or_error(
            obj.has_pool_type, "has_pool_type", obj.identifier, context="filter store insert"
        )
        if obj.is_present != {True}:
            return
        rs = get_runtime_state(obj, env)
        if rs.lock.count > 0 or len(rs.lock.queue) > 0:
            return
        store = cls.get_filter_store(pool_type, env)
        identifiers = [item.identifier for item in store.items]
        if obj.identifier in identifiers:
            deduped: list[LabObject] = []
            seen: set[str] = set()
            for item in store.items:
                if item.identifier in seen:
                    continue
                seen.add(item.identifier)
                deduped.append(item)
            store.items[:] = deduped
            return
        store.put(obj)

    @classmethod
    def safe_put_obj_into_filter_store(
        cls,
        obj: LabObject,
        env: simpy.Environment,
        *,
        context: str,
    ) -> None:
        """
        Best-effort reinsertion that never raises and records diagnostics.

        If the object's pool_type is missing/invalid, ensure the object is removed
        from all filter stores to avoid ghost entries.
        """
        try:
            if getattr(obj, "is_present", {False}) != {True}:
                cls.remove_obj_from_filter_store(obj, env)
                return
            if not getattr(obj, "has_pool_type", None):
                cls.remove_obj_from_filter_store(obj, env)
                logger.warning(
                    f"FilterStore reinsertion skipped ({context}) for {obj.identifier}: "
                    "missing pool_type"
                )
                return
            try:
                require_singleton_or_error(
                    obj.has_pool_type,
                    "has_pool_type",
                    obj.identifier,
                    context=f"filter store insert ({context})",
                )
            except Exception as exc:
                cls.remove_obj_from_filter_store(obj, env)
                logger.warning(
                    f"FilterStore reinsertion failed ({context}) for {obj.identifier}: {exc}"
                )
                return
            try:
                cls.put_obj_into_filter_store(obj, env)
            except Exception as exc:
                cls.remove_obj_from_filter_store(obj, env)
                logger.warning(
                    f"FilterStore reinsertion failed ({context}) for {obj.identifier}: {exc}"
                )
        except Exception as exc:
            logger.warning(
                f"FilterStore reinsertion failed ({context}) for {getattr(obj, 'identifier', 'unknown')}: {exc}"
            )

    @classmethod
    def remove_obj_from_filter_store(cls, obj: LabObject, env: simpy.Environment):
        ctx = get_runtime_context(env)
        for store in ctx.filter_stores.values():
            if store and store.items:
                store.items[:] = [item for item in store.items if item.identifier != obj.identifier]


class SelectorCancelled(RuntimeError):
    """Raised when a selector is cancelled so callers can treat it as an interrupt."""


class Selector(ABC):
    """
    Abstract base for all participant selectors.

    A selector behaves like a *SimPy process*.  Call it via
        `yield env.process(selector.resolve(env))`
    to obtain ``(iri, lock_request_event)`` where ``lock_request_event``
    is the successful ``simpy.Resource.request()`` that keeps the object
    reserved until you release it.
    """

    @abstractmethod
    def resolve(
            self, env: simpy.Environment
    ) -> Generator[simpy.events.Event, None, Tuple[str, simpy.events.Event]]:
        ...

    # ................................................................. #
    # Shared helper: atomic get ⟶ lock ⟶ return
    # ................................................................. #
    def _atomic_get_and_lock(
            self,
            env: simpy.Environment,
            store: simpy.FilterStore,
            pred: Callable[[LabObject], bool],
    ) -> Generator[simpy.events.Event, None, Tuple[str, simpy.events.Event]]:
        """
        Wait in the *FilterStore* until an object fulfils `pred`, then
        request its lock.  If the object changes state between `get`
        and `request`, put it back and retry.

        FIX – starvation: when a candidate fails under the lock we
        re-insert it into the store via `store.put(...)` so it can be
        reselected later (FIFO append).

        Cancellation/interrupt safety
        -----------------------------
        This generator can become "orphaned" if a parent pre_act() process is
        interrupted while awaiting selector resolution. In that case we must
        ensure we:
          * release (or cancel) any in-flight lock request, and
          * return any temporarily removed candidate back to the pool store,
        and then exit *without* raising, so the SimPy environment does not crash.
        """
        obj: LabObject | None = None
        req: simpy.events.Event | None = None
        get_ev: simpy.events.Event | None = None

        def _cleanup_selection() -> None:
            try:
                if get_ev is not None and not getattr(get_ev, "triggered", True) and hasattr(get_ev, "cancel"):
                    get_ev.cancel()
            except Exception:
                pass
            try:
                if req is not None:
                    resource = getattr(req, "resource", None)
                    if resource is not None:
                        users = getattr(resource, "users", None)
                        queue = getattr(resource, "queue", None)

                        if users is not None and req in users:
                            resource.release(req)
                        elif hasattr(req, "cancel"):
                            req.cancel()
                        elif queue is not None and req in queue:
                            try:
                                queue.remove(req)
                            except ValueError:
                                pass
            finally:
                if obj is not None:
                    FilterStoreRegistry.put_obj_into_filter_store(obj, env)

        try:
            while True:
                # 1) wait until a suitable candidate appears in the pool
                get_ev = store.get(filter=pred)  # blocking
                obj = yield get_ev
                get_ev = None
                rs = get_runtime_state(obj, env)

                # 2) request a *capacity-1* lock – blocks if already taken
                req = rs.lock.request()
                yield req

                # 3) re-validate predicate under the lock; if it still
                #    passes we are done, otherwise roll back and retry
                try:
                    predicate_ok = pred(obj)
                except Exception:
                    predicate_ok = False
                    raise
                finally:
                    if not predicate_ok:
                        rs.lock.release(req)
                        req = None
                        FilterStoreRegistry.put_obj_into_filter_store(obj, env)
                        obj = None

                if predicate_ok:
                    return obj.identifier, req

        except simpy.Interrupt:
            _cleanup_selection()
            return
        except Exception:
            _cleanup_selection()
            raise

class LiteralSelector(Selector):
    """Always returns the exact IRI given at construction."""

    def __init__(self, iri: str):
        self._iri = iri

    def resolve(
            self, env: simpy.Environment
    ) -> Generator[simpy.events.Event, None, Tuple[str, simpy.events.Event]]:
        obj = KnowledgeGraph.get_object_from_lookup(iri=self._iri)
        if obj is None:
            obj = KnowledgeGraph.get_object_from_lookup(iri=identifier_from_iri(self._iri))
        if obj is None:
            raise ValueError(f"LiteralSelector could not resolve IRI {self._iri!r}")
        obj: LabObject

        store: simpy.FilterStore | None = None
        get_ev: simpy.events.Event | None = None
        removed_from_store = False
        req: simpy.events.Event | None = None
        store_obj: LabObject | None = None
        rs = None

        try:
            pool_type = None
            if obj.has_pool_type:
                pool_type = require_singleton_or_error(
                    obj.has_pool_type, "has_pool_type", obj.identifier, context="literal selector"
                )
            if pool_type is None:
                if getattr(obj, "is_present", {False}) != {True}:
                    raise ValueError(f"LiteralSelector cannot select non-present IRI {self._iri!r}")
                obj_for_lock = obj
            else:
                store = FilterStoreRegistry.get_filter_store(pool_type, env)
                identifier = obj.identifier
                in_store = any(
                    getattr(candidate, "identifier", None) == identifier
                    for candidate in store.items
                )
                if in_store:
                    get_ev = store.get(
                        filter=lambda candidate: getattr(candidate, "identifier", None) == identifier
                    )
                    if getattr(get_ev, "triggered", False):
                        store_obj = get_ev.value
                    else:
                        store_obj = yield get_ev
                    get_ev = None
                    removed_from_store = True
                    obj_for_lock = store_obj
                else:
                    if obj.is_present == {True}:
                        rs = get_runtime_state(obj, env)
                        lock_busy = rs.lock.count > 0 or bool(rs.lock.queue)
                        if not lock_busy:
                            obj_for_lock = obj
                            if store.items:
                                store.items[:] = [
                                    item
                                    for item in store.items
                                    if getattr(item, "identifier", None) != identifier
                                ]
                        else:
                            get_ev = store.get(
                                filter=lambda candidate: getattr(candidate, "identifier", None) == identifier
                            )
                            if getattr(get_ev, "triggered", False):
                                store_obj = get_ev.value
                            else:
                                store_obj = yield get_ev
                            get_ev = None
                            removed_from_store = True
                            obj_for_lock = store_obj
                    else:
                        get_ev = store.get(
                            filter=lambda candidate: getattr(candidate, "identifier", None) == identifier
                        )
                        if getattr(get_ev, "triggered", False):
                            store_obj = get_ev.value
                        else:
                            store_obj = yield get_ev
                        get_ev = None
                        removed_from_store = True
                        obj_for_lock = store_obj

            if rs is None:
                rs = get_runtime_state(obj_for_lock, env)
            req = rs.lock.request()
            yield req
            return obj.identifier, req

        except simpy.Interrupt:
            # If interrupted while blocked on store.get(...), cancel the pending get
            # so a future store.put(obj) doesn't get consumed by a stale request.
            try:
                if get_ev is not None and not getattr(get_ev, "triggered", True) and hasattr(get_ev, "cancel"):
                    get_ev.cancel()
            except Exception:
                pass
            # Roll back any partial reservation and exit cleanly so orphaned
            # selector processes can't crash the environment.
            if req is not None:
                resource = getattr(req, "resource", None)
                if resource is not None:
                    users = getattr(resource, "users", None)
                    queue = getattr(resource, "queue", None)
                    if users is not None and req in users:
                        resource.release(req)
                    elif hasattr(req, "cancel"):
                        req.cancel()
                    elif queue is not None and req in queue:
                        try:
                            queue.remove(req)
                        except ValueError:
                            pass

            if removed_from_store:
                FilterStoreRegistry.put_obj_into_filter_store(store_obj or obj, env)
            return
        except Exception:
            try:
                if get_ev is not None and not getattr(get_ev, "triggered", True) and hasattr(get_ev, "cancel"):
                    get_ev.cancel()
            except Exception:
                pass
            if req is not None:
                resource = getattr(req, "resource", None)
                if resource is not None:
                    users = getattr(resource, "users", None)
                    queue = getattr(resource, "queue", None)
                    if users is not None and req in users:
                        resource.release(req)
                    elif hasattr(req, "cancel"):
                        req.cancel()
                    elif queue is not None and req in queue:
                        try:
                            queue.remove(req)
                        except ValueError:
                            pass
            if removed_from_store:
                FilterStoreRegistry.put_obj_into_filter_store(store_obj or obj, env)
            raise

    def __str__(self) -> str:
            return f'LiteralSelector("{self._iri}")'


class RuntimeSelector(Selector):
    def __init__(self, pool_type: str, predicate: Callable[[LabObject], bool]):
        self.pool_type = pool_type
        self._predicate = predicate
        # store predicate source for replay / provenance
        try:
            self._source = inspect.getsource(predicate).strip()
        except OSError:
            self._source = repr(predicate)

    def resolve(
            self, env: simpy.Environment
    ) -> Generator[simpy.events.Event, None, Tuple[str, simpy.events.Event]]:
        store = FilterStoreRegistry.get_filter_store(self.pool_type, env)
        return (yield from self._atomic_get_and_lock(env, store, self._predicate))


class AttributeSelector(RuntimeSelector):
    """
    Pick any LabObject of *pool_type* for which ``predicate(obj)`` is True.

    Parameters
    ----------
    pool_type:
        The `LabObject.has_pool_type` value that defines the SimPy
        `FilterStore` from which to draw candidates.
    predicate:
        A function that inspects the LabObject’s *current* attributes.
        It **must not perform blocking calls**; keep it cheap.
    """

    def __str__(self) -> str:
        return f"AttributeSelector(pool={self.pool_type}, pred={self._source})"


class HistorySelector(RuntimeSelector):
    """
    Like *AttributeSelector* but allows the predicate to take the
    LabObject’s recent `UnitaryEdit`s into account.

    The selector itself does **not** peek into history – you encode that
    logic in the predicate:

    >>> def cleaned_twice(obj, env):
    ...     rs = get_runtime_state(obj, env)
    ...     return sum(isinstance(a, Clean) for a in rs.recent_actions) >= 2
    """

    def __str__(self) -> str:
        return f"HistorySelector(pool={self.pool_type}, pred={self._source})"

    def __init__(self, pool_type: str, predicate: Callable[..., bool]):
        super().__init__(pool_type, predicate)
        self._predicate_accepts_env = False
        try:
            sig = inspect.signature(predicate)
            has_varargs = any(
                param.kind is inspect.Parameter.VAR_POSITIONAL for param in sig.parameters.values()
            )
            positional = [
                param
                for param in sig.parameters.values()
                if param.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)
            ]
            if has_varargs or len(positional) >= 2:
                self._predicate_accepts_env = True
        except (TypeError, ValueError):
            self._predicate_accepts_env = False

    def resolve(
            self, env: simpy.Environment
    ) -> Generator[simpy.events.Event, None, Tuple[str, simpy.events.Event]]:
        store = FilterStoreRegistry.get_filter_store(self.pool_type, env)
        if self._predicate_accepts_env:
            predicate = lambda obj: self._predicate(obj, env)
        else:
            predicate = self._predicate
        return (yield from self._atomic_get_and_lock(env, store, predicate))


class KgQuerySelector(Selector):
    """
    Run a SPARQL query over the union graph (KG + overlays) and pick a candidate
    from a FilterStore pool.
    """

    def __init__(self, pool_type: str, sparql: str, var: str = "s"):
        self.pool_type = pool_type
        self.sparql = sparql
        self.var = var

    def resolve(
            self, env: simpy.Environment
    ) -> Generator[simpy.events.Event, None, Tuple[str, simpy.events.Event]]:
        engine = get_effect_engine(env)
        query_graph = engine.build_query_graph()
        var_name = self.var.lstrip("?")
        candidates: set[str] = set()
        for row in query_graph.query(self.sparql):
            value = row.asdict().get(var_name)
            if value is None:
                continue
            candidates.add(identifier_from_iri(str(value)))

        if not candidates:
            raise RuntimeError(f"KgQuerySelector query returned no candidates: {self.sparql}")

        store = FilterStoreRegistry.get_filter_store(self.pool_type, env)
        return (yield from self._atomic_get_and_lock(
            env,
            store,
            lambda obj: identifier_from_iri(obj.identifier) in candidates,
        ))  # iri, req

    def __str__(self) -> str:
        return f"KgQuerySelector(pool={self.pool_type}, var={self.var})"


__all__ = [
    "Selector",
    "LiteralSelector",
    "AttributeSelector",
    "HistorySelector",
    "KgQuerySelector",
    "SelectorCancelled",
    "FilterStoreRegistry"
]
