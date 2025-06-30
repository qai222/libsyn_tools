from __future__ import annotations

from collections.abc import Iterable
from typing import List

from loguru import logger

from libsyn_tools.sim.operation.operation import Operation, UnitaryEdit, UnitaryEditType


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

    def apply(self, edits: Iterable[UnitaryEdit]) -> None:
        """
        Apply each `UnitaryEdit` **in order**.

        If an edit fails, everything applied so far is *rolled back* and
        the original exception bubbles up.
        """
        applied: list[UnitaryEdit] = []
        try:
            for edit in edits:
                logger.debug(f"Applying edit: {edit.type} – {edit.instance_1_iri}")
                edit.apply()
                applied.append(edit)
        except Exception as exc:  # pragma: no cover – transaction abort
            logger.error(f"Edit failed, rolling back {len(applied)} edits")
            self.rollback(applied)
            raise

    def rollback(self, applied_edits: List[UnitaryEdit]) -> None:
        """Undo *applied_edits* in **reverse order** (best-effort)."""
        for edit in reversed(applied_edits):
            try:
                inverse = self._inverse_edit(edit)
                logger.debug(f"Rollback: {inverse.type} – {inverse.instance_1_iri}")
                inverse.apply()
            except Exception as exc:  # pragma: no cover
                # We *never* raise from here – a rollback must not cascade.
                logger.error(f"Rollback failed for {edit}: {exc!r}")

    @staticmethod
    def _inverse_edit(edit: UnitaryEdit) -> UnitaryEdit:
        """Return a **best-effort** inverse of *edit*."""
        t = UnitaryEditType
        match edit.type:
            case t.CREATE:
                return edit.model_copy(update={"type": t.ANNIHILATE})
            case t.ANNIHILATE:
                return edit.model_copy(update={"type": t.CREATE})
            case t.CHANGE_DATA_PROPERTY:
                # Without state diff we cannot deterministically invert –
                # in production you would record the *before* value.
                raise NotImplementedError(
                    "Inverse for data-property change requires original value"
                )
            case t.ADD_OBJECT_PROPERTY:
                return edit.model_copy(update={"type": t.REMOVE_OBJECT_PROPERTY})
            case t.REMOVE_OBJECT_PROPERTY:
                return edit.model_copy(update={"type": t.ADD_OBJECT_PROPERTY})
            case _:
                raise ValueError(f"Unknown edit type: {edit.type}")
