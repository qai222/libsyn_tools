from typing import Callable, Any

from .knowledge_graph.physical_entities import LabObject, Unit_Registry

Predicate = Callable[[Any], bool]


class Invariant:
    PHASE_PRE = 1
    PHASE_IN = 2
    PHASE_POST = 4

    def __init__(self, name: str, predicate: Predicate, phase_mask: int):
        self.name = name
        self.predicate = predicate
        self.phase_mask = phase_mask

    def check(self, ctx: dict):
        if not self.predicate(ctx):
            raise RuntimeError(f"Invariant '{self.name}' violated – ctx: {ctx}")


class InvariantEngine:
    def __init__(self, invariants: list[Invariant] | None = None):
        self.invariants = invariants or []

    def validate(self, phase: int, ctx: dict):
        for inv in self.invariants:
            if inv.phase_mask & phase:
                inv.check(ctx)


def default_invariants():
    def no_overfill(ctx):
        lb: LabObject | None = ctx.get("lab_object")
        if not lb or lb.capacity is None:
            return True
        return lb.current_volume() <= lb.capacity + 1e-9 * Unit_Registry.millilitre

    return [Invariant("no_overfill", no_overfill, Invariant.PHASE_POST)]
