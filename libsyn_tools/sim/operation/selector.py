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
"""

from __future__ import annotations

import inspect
from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Callable
from typing import Generator, Tuple, Dict

import simpy

from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject
from libsyn_tools.sim.operation.runtime import get_runtime_state


class FilterStoreRegistry:
    _stores: Dict[str, simpy.FilterStore] = defaultdict()

    @classmethod
    def get_filter_store(cls, pool_type: str, env: simpy.Environment | None = None) -> simpy.FilterStore:
        if pool_type not in cls._stores:
            if env is None:
                raise RuntimeError("env required on first call for a pool_type")
            cls._stores[pool_type] = simpy.FilterStore(env)
        return cls._stores[pool_type]

    @classmethod
    def put_obj_into_filter_store(cls, obj: LabObject, env: simpy.Environment):
        pool_type = next(iter(obj.has_pool_type), None)
        if pool_type is not None:
            store = cls.get_filter_store(pool_type, env)
            if obj not in store.items:
                store.put(obj)


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
        """
        while True:
            # 1) wait until a suitable candidate appears in the pool
            obj: LabObject = yield store.get(filter=pred)  # blocking
            rs = get_runtime_state(obj, env)

            # 2) request a *capacity-1* lock – blocks if already taken
            req = rs.lock.request()
            yield req

            # 3) re-validate predicate under the lock; if it still
            #    passes we are done, otherwise roll back and retry
            if pred(obj):
                return obj.identifier, req

            rs.lock.release(req)
            store.put(obj)  # put back so others may pick it, note the item loses FIFO position (goes to tail).


class LiteralSelector(Selector):
    """Always returns the exact IRI given at construction."""

    def __init__(self, iri: str):
        self._iri = iri

    def resolve(
            self, env: simpy.Environment
    ) -> Generator[simpy.events.Event, None, Tuple[str, simpy.events.Event]]:
        obj = LabObject.object_lookup[self._iri]
        obj: LabObject
        rs = get_runtime_state(obj, env)
        req = rs.lock.request()
        yield req
        return obj.identifier, req

    # Human-readable representation, logged for provenance
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
        iri, req = yield from self._atomic_get_and_lock(env, store, self._predicate)

        return iri, req


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


__all__ = [
    "Selector",
    "LiteralSelector",
    "AttributeSelector",
    "HistorySelector",
]
