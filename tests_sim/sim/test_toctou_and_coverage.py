# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_toctou_and_coverage.py ###
from __future__ import annotations

import pytest
import simpy
from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph  import (
    LabObject,
    Has_interrupt_events,
    Is_directly_contained_by
)
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.selector import AttributeSelector, FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, AddDataProperty, AddObjectProperty, Create


# ---------- TOCTOU race: re-check fails under lock & selector picks second vial ----------

class MarkDirty(Operation):
    """Locks a specific vial, marks it 'DIRTY' via has_interrupt_events."""

    participant_vial: str = Field(...)

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            AddDataProperty(
                instance_1_iri=self.participant_vial,
                property_iri=Has_interrupt_events.predicate_iri,
                data_value="DIRTY",
            )
        ]


class PickClean(Operation):
    """Picks any vial from pool 'VIAL' whose has_interrupt_events does not contain 'DIRTY'."""

    participant_vial: AttributeSelector = Field(
        default_factory=lambda: AttributeSelector(
            pool_type="VIAL",
            predicate=lambda obj: ("DIRTY" not in getattr(obj, "has_interrupt_events", set()))
        )
    )

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []


# replace the helper with:
def _make_vial(pool="VIAL") -> LabObject:
    v = LabObject()
    v.has_pool_type.add(pool)
    KnowledgeGraph.get_object_from_lookup(v.identifier)
    # Mark the vial as present so post_act will reinsert it into the pool
    Create(instance_1_iri=v.identifier).apply()
    return v


def test_toctou_recheck_fails_then_picks_other(env: simpy.Environment):
    # Two vials in the pool
    v1 = _make_vial()
    v2 = _make_vial()

    # A will lock v1 and mark it DIRTY; B tries to pick a "clean" vial from the pool
    a = MarkDirty(identifier="A", participant_vial=v1.identifier, temporal_cost=0.0)
    b = PickClean(identifier="B", temporal_cost=0.0)

    # Simulation will build resources + pools for all LabObjects in memory
    sim = Simulation([a, b])
    sim.run()

    # B should have resolved to v2 (not v1) because re-check under lock fails on v1
    assert b.resolved_resources.get("vial") == v2.identifier

    # Optional sanity: pool still contains both vials after post_act normalization
    store = FilterStoreRegistry.get_filter_store("VIAL", sim.env)
    assert v1 in store.items and v2 in store.items


# ---------- Coverage enforcement for object-property writes ----------

class BadWrite(Operation):
    """
    Tries to write an object property between two lab objects while
    only one endpoint is locked (and neither is created in this batch).
    Should fail mechanical pre-check (coverage).
    """

    participant_src: str = Field(...)
    # NO participant_dst on purpose — not locked
    dst_iri: str = Field(..., description="Destination IRI (intentionally not locked)")

    def get_operation_effects(self) -> list[UnitaryEdit]:
        # write: src --is_directly_contained_by--> dst (dst not locked/created)
        return [
            AddObjectProperty(
                instance_1_iri=self.participant_src,
                instance_2_iri=self.dst_iri,
                property_iri=Is_directly_contained_by.predicate_iri,
            )
        ]


def test_coverage_check_blocks_unlocked_endpoint(env: simpy.Environment):
    src = LabObject()
    dst = LabObject()
    KnowledgeGraph.get_object_from_lookup(src.identifier)
    KnowledgeGraph.get_object_from_lookup(dst.identifier)

    # Put only src into resources via Simulation build (dst remains unlocked)
    op = BadWrite(participant_src=src.identifier, dst_iri=dst.identifier, temporal_cost=0.0)
    sim = Simulation([op])

    # Engine mechanical pre-check must abort; ensure no hang
    sim.run()

    assert sim.effect_engine._shacl_violations
    types = [r.event_type for r in sim.history_log]
    assert "OPERATION_START" in types and "OPERATION_END" not in types
    assert "OPERATION_ABORT" in types
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_toctou_and_coverage.py ###
