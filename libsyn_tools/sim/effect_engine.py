from __future__ import annotations

from collections.abc import Iterable
from typing import List

import simpy
from loguru import logger

from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject, BaseClass
from libsyn_tools.sim.operation import Operation, UnitaryEdit, UnitaryEditType, FilterStoreRegistry, get_runtime_state
from libsyn_tools.sim.operation.runtime import _RESOURCE_MAP, _RUNTIME_CACHE


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
                inverse = edit.compute_inverse()

                # auto register simpy resource
                obj = BaseClass.object_lookup[edit.instance_1_iri]
                if edit.type == UnitaryEditType.CREATE:
                    self._register_if_new(obj, env)
                elif edit.type == UnitaryEditType.ANNIHILATE:
                    self._unregister_object(obj)

                logger.debug(f"Applying edit: {edit.type} – {edit.instance_1_iri}")
                edit.apply()
                applied.append(edit)
                inverses.append(inverse)

                self._sync_filter_stores(obj, env)

                rs = get_runtime_state(BaseClass.object_lookup[edit.instance_1_iri], env)
                rs.recent_edits.append(edit)

        except Exception as exc:  # pragma: no cover – transaction abort
            logger.error(f"Edit failed, rolling back {len(applied)} edits")
            self.rollback(inverses, env)
            raise EditApplicationError(edit, exc) from exc

    def rollback(self, inverses: List[UnitaryEdit], env: simpy.Environment) -> None:
        """Undo *applied_edits* in **reverse order** (best-effort)."""
        for inv in reversed(inverses):
            try:
                logger.debug(f"Rollback: {inv}")
                # NEW: `CREATE` during rollback may need registering too
                obj = BaseClass.object_lookup[inv.instance_1_iri]
                if inv.type is UnitaryEditType.CREATE:
                    self._register_if_new(obj, env)
                elif inv.type is UnitaryEditType.ANNIHILATE:
                    self._unregister_object(obj)
                self._sync_filter_stores(obj, env)
                inv.apply()
            except Exception as exc:  # pragma: no cover
                logger.error(f"Rollback failed for {inv}: {exc!r}")

    @staticmethod
    def _register_if_new(obj: LabObject, env: simpy.Environment):
        if obj.instance_iri not in _RESOURCE_MAP:
            _RESOURCE_MAP[obj.instance_iri] = simpy.Resource(env, capacity=1)
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)

    def _unregister_object(self, obj: LabObject):
        """Remove *all* runtime artefacts for a vanished object."""
        _RESOURCE_MAP.pop(obj.instance_iri, None)
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        _RUNTIME_CACHE.pop(obj.instance_iri, None)

    def _sync_filter_stores(self, obj: LabObject, env: simpy.Environment):
        """
        After *any* edit we normalise FilterStore membership so rollbacks
        never leave an object in the wrong pool.
        """
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        if obj.is_present == {True}:  # ← guard
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)
