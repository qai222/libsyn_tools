from __future__ import annotations

import datetime as _dt
from typing import List, Dict

from loguru import logger
from pydantic import BaseModel, Field

from .action.core import UnitaryEdit, UnitaryEditType
from .knowledge_graph import Individual


class GraphEditRecord(BaseModel):
    """A single committed *edit → inverse* pair with timestamp & actor."""

    t_sim: float  # simulation time
    t_wall: _dt.datetime = Field(default_factory=_dt.datetime.utcnow)
    action_id: str
    edit: UnitaryEdit
    inverse: UnitaryEdit | None = None

    class Config:
        arbitrary_types_allowed = True


class Snapshot(BaseModel):
    """Lightweight KG snapshot – maps *IRI → serialised instance*."""

    t_sim: float
    triples: Dict[str, dict]  # Individual JSON representations


class SnapshotManager:
    """Keeps a rolling history of edits + sparse snapshots for time‑travel queries."""

    def __init__(self, snapshot_interval: float = 60.0):
        self.edits: List[GraphEditRecord] = []
        self.snapshots: List[Snapshot] = []
        self._interval = snapshot_interval

    # ------------------------------------------------------------------
    def record_edit(self, t_sim: float, action_id: str, edit: UnitaryEdit, inverse: UnitaryEdit | None):
        self.edits.append(GraphEditRecord(t_sim=t_sim, action_id=action_id, edit=edit, inverse=inverse))

    def maybe_snapshot(self, t_sim: float, root_kg: Dict[str, Individual]):
        if not self.snapshots or (t_sim - self.snapshots[-1].t_sim) >= self._interval:
            triples = {iri: ind.model_dump(mode="json") for iri, ind in root_kg.items()}
            self.snapshots.append(Snapshot(t_sim=t_sim, triples=triples))
            logger.debug(f"[audit] Snapshot taken at t={t_sim:.1f}; size={len(triples):,}")

    # ------------------------------------------------------------------
    def rewind(self, t_sim: float) -> Dict[str, dict]:
        """Return *logical* KG state at **t_sim** (immutable JSON)."""
        # 1) find latest snapshot ≤ t_sim
        snap = max((s for s in self.snapshots if s.t_sim <= t_sim), key=lambda s: s.t_sim, default=None)
        if not snap:
            raise RuntimeError("No snapshot before requested time")
        state: Dict[str, dict] = {**snap.triples}  # shallow copy
        # 2) replay edits between snapshot and target time
        for rec in (e for e in self.edits if snap.t_sim < e.t_sim <= t_sim):
            _apply_json_edit(state, rec.edit)
        return state


def _apply_json_edit(state: Dict[str, dict], edit: UnitaryEdit):
    """Pure‑function JSON patch so audit layer is decoupled from live objects."""
    match edit.type:
        case UnitaryEditType.CREATE:
            state[edit.instance_1_iri] = {"is_present": True}
        case UnitaryEditType.ANNIHILATE:
            if edit.instance_1_iri in state:
                state[edit.instance_1_iri]["is_present"] = False
        # NOTE: Other edit kinds omitted for brevity – implement analogously.
        case _:
            pass
