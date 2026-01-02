# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_effects_dsl.py ###
from __future__ import annotations

from libsyn_tools.sim.operation.effects_dsl import EffectsBuilder
from libsyn_tools.sim.operation.unitary_edit import (
    AddDataProperty,
    AddObjectProperty,
    Annihilate,
    ChangeDataProperty,
    Create,
    RemoveDataProperty,
    RemoveObjectProperty,
)


def test_effects_builder_orders_edits() -> None:
    edits = (
        EffectsBuilder()
        .create("ex:obj")
        .set_data("ex:obj", "ex:prop", 1.0)
        .add_data("ex:obj", "ex:tag", "a")
        .remove_data("ex:obj", "ex:tag", "a")
        .link("ex:obj", "ex:rel", "ex:other")
        .unlink("ex:obj", "ex:rel", "ex:other")
        .annihilate("ex:obj")
        .build()
    )

    assert [type(edit) for edit in edits] == [
        Create,
        ChangeDataProperty,
        AddDataProperty,
        RemoveDataProperty,
        AddObjectProperty,
        RemoveObjectProperty,
        Annihilate,
    ]
    assert edits[0].instance_1_iri == "ex:obj"
    assert edits[1].property_iri == "ex:prop"
    assert edits[1].data_value == 1.0
    assert edits[2].property_iri == "ex:tag"
    assert edits[2].data_value == "a"
    assert edits[3].property_iri == "ex:tag"
    assert edits[3].data_value == "a"
    assert edits[4].instance_2_iri == "ex:other"
    assert edits[5].instance_2_iri == "ex:other"
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_effects_dsl.py ###
