from __future__ import annotations

from collections.abc import Iterable
from typing import List

from loguru import logger

from .action.core import Action, UnitaryEdit

"""
Effect Engine for libsyn_tools.sim
====================================

All side‑effects are centralised in one place so the DES layer remains domain‑agnostic and easier to test.

Key design goals
----------------
1.  **Single‑responsibility** – `ActionProcess` now orchestrates timing only;
    `EffectEngine` applies or rolls back edits.
2.  **Atomicity** – an action’s edits are staged, then committed; a failure
    triggers `rollback()` on the subset that succeeded.
3.  **Extensibility** – later we can plug in alternate back‑ends (mock KG,
    audit‑only mode, etc.) by subclassing this engine.
"""


class EffectEngine:

    def prepare(self, action: Action) -> List[UnitaryEdit]:
        """Populate and return the *staged* edits for *action*.

        `Action.pre_act()` is intentionally called here – not inside the SimPy
        process – so that **all graph reads happen outside timing logic**.
        """
        action.pre_act()
        return action.action_effects.copy()

    def apply(self, edits: Iterable[UnitaryEdit]) -> None:
        """Apply each `UnitaryEdit` in order and log progress."""
        for edit in edits:
            logger.debug(f"Applying edit: {edit.type} – {edit.instance_1_iri}")
            edit.apply()

    def finalize(self, action: Action) -> None:
        """Hook called after *all* edits commit successfully."""
        action.post_act()

    # ------------------------------------------------------------------
    # Rollback support – simplistic 180° reversal.  A richer implementation
    # could store inverse edits explicitly.
    # ------------------------------------------------------------------
    def rollback(self, applied_edits: List[UnitaryEdit]) -> None:  # noqa: D401
        """Undo *applied_edits* in **reverse order** (best‑effort).

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

    # ---------------------- helper utilities -------------------------
    @staticmethod
    def _inverse_edit(edit: UnitaryEdit) -> UnitaryEdit:
        """Return a *best‑effort* inverse of *edit* (utility)."""
        from .action.core import UnitaryEditType as T  # local import avoids cycle

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
