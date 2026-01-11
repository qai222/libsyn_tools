from __future__ import annotations

from libsyn_tools.sim.knowledge_graph import identifier_from_iri
from libsyn_tools.sim.operation_preset.drain import DrainExcess
from twa.data_model.base_ontology import KnowledgeGraph


_EPS = 1e-6


def make_drain_to_capacity(focus_iri: str, waste_iri: str) -> DrainExcess:
    focus_id = identifier_from_iri(focus_iri)
    waste_id = identifier_from_iri(waste_iri)
    focus = KnowledgeGraph.get_object_from_lookup(focus_iri)
    if focus is None:
        focus = KnowledgeGraph.get_object_from_lookup(focus_id)
    waste = KnowledgeGraph.get_object_from_lookup(waste_iri)
    if waste is None:
        waste = KnowledgeGraph.get_object_from_lookup(waste_id)
    if focus is None or waste is None:
        raise RuntimeError("make_drain_to_capacity requires focus and waste to exist in KG")
    if not getattr(focus, "has_capacity", set()):
        raise RuntimeError(f"{focus.identifier} has no capacity set")
    if len(focus.has_capacity) != 1:
        raise RuntimeError(f"{focus.identifier} has invalid capacity set")
    try:
        capacity = getattr(focus, "capacity", None)
    except ValueError as exc:
        raise RuntimeError(f"{focus.identifier} has invalid capacity") from exc
    if capacity is None:
        raise RuntimeError(f"{focus.identifier} has no capacity set")
    if capacity <= 0:
        raise RuntimeError(f"{focus.identifier} has non-positive capacity")
    target_volume = capacity
    if capacity > _EPS:
        target_volume = max(capacity - _EPS, _EPS)
    target_volume = min(target_volume, capacity)
    return DrainExcess(
        participant_source=focus.identifier,
        participant_destination=waste.identifier,
        participant_device=waste.identifier,
        target_volume=target_volume,
    )
