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
    capacity = getattr(focus, "capacity", None)
    if capacity is None:
        raise RuntimeError(f"{focus.identifier} has no capacity set")
    if capacity <= _EPS:
        raise ValueError(f"{focus.identifier} has non-positive capacity {capacity}")
    target_volume = min(max(capacity - _EPS, _EPS), capacity)
    return DrainExcess(
        participant_source=focus.identifier,
        participant_destination=waste.identifier,
        participant_device=waste.identifier,
        target_volume=target_volume,
    )
