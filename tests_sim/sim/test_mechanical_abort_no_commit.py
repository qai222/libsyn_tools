# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_mechanical_abort_no_commit.py ###
from __future__ import annotations

import simpy
from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph  import Has_interrupt_events
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, AddDataProperty


class BadOp(Operation):
    """Operation that tries to write to a non-existent IRI (mechanical abort)."""

    name: str = Field("bad")

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return [
            AddDataProperty(
                instance_1_iri="LabObject_FAKE",  # not created in this batch
                property_iri=Has_interrupt_events.predicate_iri,
                data_value="boom",
            )
        ]


def test_mechanical_abort_raises_and_no_commit(env: simpy.Environment):
    sim = Simulation([BadOp(identifier="bad1")])

    # EffectEngine pre-check should abort the operation and record violation
    sim.run()

    # Event log contains START but no END for the bad op
    types = [r.event_type for r in sim.history_log]
    assert "OPERATION_START" in types
    assert "OPERATION_END" not in types
    assert "OPERATION_ABORT" in types

    assert any(v.origin == "ENGINE" for v in sim.effect_engine._shacl_violations)

    # No effect applied: either the object doesn't exist, or it exists but
    # is not present and does not have the data value we tried to add.
    try:
        obj = KnowledgeGraph.get_object_from_lookup("LabObject_FAKE")
    except Exception:
        # Lookup fails → that's fine (definitely no commit happened)
        return

    # If lookup succeeds, assert no data was written and not present
    # (engine aborted before commit).
    has_vals = getattr(obj, "has_interrupt_events", set())
    assert "boom" not in has_vals
    assert getattr(obj, "is_present", {False}) != {True}
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_mechanical_abort_no_commit.py ###
