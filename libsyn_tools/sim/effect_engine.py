from __future__ import annotations

from collections.abc import Iterable
from typing import List

import simpy
from loguru import logger
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph.physical_entities import BaseClass
from libsyn_tools.sim.operation import Operation, UnitaryEdit, UnitaryEditType, FilterStoreRegistry, get_runtime_state
from libsyn_tools.sim.operation.runtime import _RESOURCE_MAP, _RUNTIME_CACHE, _needs_runtime_tracking


class EditApplicationError(RuntimeError):
    """Wraps the original exception + the offending edit for richer traceback."""

    def __init__(self, edit: UnitaryEdit, original: Exception):
        self.edit = edit
        self.original = original
        super().__init__(f"Failed applying {edit}: {original!r}")


class EffectEngine:
    """
    This module provides **transaction-like semantics** around
    `UnitaryEdit` application:

    1.  *prepare()* – return a **copy** of the staged edits that an
        `Operation` built in its *pre-act* phase.
    2.  *apply()*  – iterate through the edits **atomically**.
        If any single edit raises, all previously-applied edits are rolled
        back in *reverse* order and the original exception is re-raised.
    3.  *rollback()* – best-effort inverse for each edit (utility).

    The class is intentionally *stateless* so one instance can be shared by
    a whole `Simulation`.
    """

    def prepare(self, action: Operation) -> List[UnitaryEdit]:
        """Return a **defensive copy** of the staged edits for *action*.

        The `Operation` must already have executed its *pre-act* phase
        so the list is fully populated.
        """
        return action.operation_effects.copy()

    def apply(self, edits: Iterable[UnitaryEdit], env: simpy.Environment) -> None:
        """
        Apply each `UnitaryEdit` **in order**.

        If an edit fails, everything applied so far is *rolled back* and
        the original exception bubbles up.
        """
        applied: list[UnitaryEdit] = []
        inverses: list[UnitaryEdit] = []

        try:
            for edit in edits:

                subj = KnowledgeGraph.get_object_from_lookup(edit.instance_1_iri)
                obj2 = (
                    KnowledgeGraph.get_object_from_lookup(edit.instance_2_iri)
                    if edit.type in (UnitaryEditType.ADD_OBJECT_PROPERTY,
                                     UnitaryEditType.REMOVE_OBJECT_PROPERTY)
                       and edit.instance_2_iri
                    else None
                )

                # ensure every referenced LabObject owns a Resource
                if edit.type is UnitaryEditType.CREATE:
                    pass  # postpone until after .apply()
                elif obj2:  # only for ADD, not REMOVE
                    self._register_if_new(obj2, env)

                logger.debug(f"Applying edit: {edit.type} – {subj.__class__.__name__}={edit.instance_1_iri}")
                inverse = edit.compute_inverse()
                edit.apply()
                inverses.append(inverse)
                applied.append(edit)

                if edit.type is UnitaryEditType.CREATE:
                    self._register_if_new(subj, env)

                self._log_and_sync(subj, edit, env)
                self._log_and_sync(obj2, edit, env)

                if edit.type is UnitaryEditType.ANNIHILATE:
                    # remove resource + filter entry only after logging
                    self._unregister_object(subj)

        except Exception as exc:  # pragma: no cover – transaction abort
            logger.error(f"Edit failed, rolling back {len(applied)} edits")
            self.rollback(inverses, env)
            raise EditApplicationError(edit, exc) from exc

    def _log_and_sync(self, obj: BaseClass, edit: UnitaryEdit,
                      env: simpy.Environment) -> None:
        """Append *edit* to obj.runtime history and normalise FilterStore."""
        if obj is None or not _needs_runtime_tracking(obj):
            return
        get_runtime_state(obj, env).recent_edits.append(edit)
        self._sync_filter_stores(obj, env)

    def rollback(self, inverses: List[UnitaryEdit], env: simpy.Environment) -> None:
        """Undo *applied_edits* in **reverse order** (best-effort)."""
        for inv in reversed(inverses):
            try:
                logger.debug(f"Rollback: {inv}")
                subj = KnowledgeGraph.get_object_from_lookup(inv.instance_1_iri)
                obj2 = (
                    KnowledgeGraph.get_object_from_lookup(inv.instance_2_iri)
                    if inv.type in {UnitaryEditType.ADD_OBJECT_PROPERTY,
                                    UnitaryEditType.REMOVE_OBJECT_PROPERTY}
                       and inv.instance_2_iri else None
                )

                inv.apply()

                if subj.is_present == {True}:                 # object exists ⇒ ensure resource & pool
                    self._register_if_new(subj, env)
                else:                                         # vanished ⇒ clean runtime artefacts
                    self._unregister_object(subj)

                self._sync_filter_stores(subj, env)
                if obj2:
                    self._sync_filter_stores(obj2, env)

            except Exception as exc:  # pragma: no cover
                logger.error(f"Rollback failed for {inv}: {exc!r}")

    @staticmethod
    def _register_if_new(obj: BaseClass, env: simpy.Environment):
        if not _needs_runtime_tracking(obj):
            return
        if obj.instance_iri not in _RESOURCE_MAP:
            _RESOURCE_MAP[obj.instance_iri] = simpy.Resource(env, capacity=1)
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)
            logger.debug(f"auto register new object: {obj.__class__.__name__}={obj.instance_iri}")

    def _unregister_object(self, obj: BaseClass):
        """Remove *all* runtime artefacts for a vanished object."""
        if not _needs_runtime_tracking(obj):
            return
        _RESOURCE_MAP.pop(obj.instance_iri, None)
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        _RUNTIME_CACHE.pop(obj.instance_iri, None)
        setattr(obj, "_runtime", None)  # clear pointer

    def _sync_filter_stores(self, obj: BaseClass, env: simpy.Environment):
        """
        After *any* edit we normalise FilterStore membership so rollbacks
        never leave an object in the wrong pool.
        """
        if not _needs_runtime_tracking(obj):
            return
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        if obj.is_present == {True}:  # ← guard
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)
