from __future__ import annotations

from collections.abc import Iterable
from typing import List

from loguru import logger

from libsyn_tools.sim.operation.operation import Operation, UnitaryEdit

"""
All effects are centralised in one place so the DES layer remains domain‑agnostic and easier to test.
"""


class EffectEngine:

    def prepare(self, action: Operation) -> List[UnitaryEdit]:
        """Populate and return the *staged* edits for *action*.

        `Action.pre_act()` is intentionally called here – not inside the SimPy
        process – so that **all graph reads happen outside timing logic**.
        """
        action.pre_act()
        return action.operation_effects.copy()

    def apply(self, edits: Iterable[UnitaryEdit]) -> None:
        """Apply each `UnitaryEdit` in order and log progress."""
        for edit in edits:
            logger.debug(f"Applying edit: {edit.type} – {edit.instance_1_iri}")
            edit.apply()

    def finalize(self, action: Operation) -> None:
        """Hook called after *all* edits commit successfully."""
        action.post_act()

    def rollback(self, applied_edits: List[UnitaryEdit]) -> None:
        """
        Undo *applied_edits* in **reverse order** (best‑effort).

        This naïve rollback assumes every edit has a deterministic inverse; for
        production use consider capturing before/after snapshots instead.
        """
        logger.warning("Rolling back %d partially‑applied edits", len(applied_edits))
        for edit in reversed(applied_edits):
            try:
                inverse = self._inverse_edit(edit)
                inverse.apply()
            except Exception as exc:  # pragma: no cover – best‑effort
                logger.error(f"Rollback failed for {edit}: {exc!r}")

    @staticmethod
    def _inverse_edit(edit: UnitaryEdit) -> UnitaryEdit:
        """Return a *best‑effort* inverse of *edit* (utility)."""
        from libsyn_tools.sim.operation.operation import UnitaryEditType as T  # local import avoids cycle

        match edit.type:
            case T.CREATE:
                return edit.model_copy(update={"type": T.ANNIHILATE})
            case T.ANNIHILATE:
                return edit.model_copy(update={"type": T.CREATE})
            case T.CHANGE_DATA_PROPERTY:
                raise NotImplementedError("Inverse for data‑property change unclear")
            case T.ADD_OBJECT_PROPERTY:
                return edit.model_copy(update={"type": T.REMOVE_OBJECT_PROPERTY})
            case T.REMOVE_OBJECT_PROPERTY:
                return edit.model_copy(update={"type": T.ADD_OBJECT_PROPERTY})
            case _:
                raise ValueError(f"Unknown edit type: {edit.type}")
