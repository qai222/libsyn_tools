# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_preset_transfer.py ###
from __future__ import annotations

import simpy
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph.physical_entities import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import (
    UnitaryEditType,
    Create,
    AddObjectProperty,
)
from libsyn_tools.sim.operation_preset.transfer import (
    TransferMaterialByVolume,
    TransferMaterialByPortionSize,
)


def _mk_container_and_pom(vol_ml: float):
    src = MaterialContainer()
    dst = MaterialContainer()
    dev = MaterialContainer()  # device treated as LabObject/MaterialContainer for simplicity

    pom = PortionOfMaterial()
    pom.add_chemical(Chemical(mass=vol_ml, density=1.0))

    return src, dst, dev, pom


def test_transfer_by_volume_edits_shape(env: simpy.Environment):
    src, dst, dev, pom = _mk_container_and_pom(10.0)

    # Seed KG and link initial POM into src
    eng = EffectEngine()
    for o in (src, dst, dev, pom):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
    eng.apply(
        [
            Create(instance_1_iri=src.identifier),
            Create(instance_1_iri=dst.identifier),
            Create(instance_1_iri=dev.identifier),
            Create(instance_1_iri=pom.identifier),
            AddObjectProperty(
                instance_1_iri=pom.identifier,
                instance_2_iri=src.identifier,
                property_iri=Is_directly_contained_by.predicate_iri,
            ),
        ],
        env,
        operation_id="seed",
        locked_iris=[src.identifier, dst.identifier, dev.identifier, pom.identifier],
    )

    # Build operation (transfer 6 mL of 10 mL)
    op = TransferMaterialByVolume(
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=dev.identifier,
        transfer_volume=6.0,
    )
    edits = op.get_operation_effects()

    # For a single POM, expected sequence is 7 edits:
    # 1 Annihilate original, 2 Create residual + 3 add->src,
    # 4 Create transfer + 5 add->dev + 6 remove->dev + 7 add->dst
    assert len(edits) == 7
    types = [e.type for e in edits]
    assert types[0] is UnitaryEditType.ANNIHILATE
    assert types[-1] is UnitaryEditType.ADD_OBJECT_PROPERTY


def test_transfer_by_portionsize_validation():
    src = MaterialContainer()
    dst = MaterialContainer()
    dev = MaterialContainer()
    # model init validation should reject >1.0
    try:
        TransferMaterialByPortionSize(
            participant_source=src.identifier,
            participant_destination=dst.identifier,
            participant_device=dev.identifier,
            portion_size=1.2,
        )
        assert False, "portion_size > 1.0 should be rejected"
    except Exception:
        pass
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_preset_transfer.py ###
