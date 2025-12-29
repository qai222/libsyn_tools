# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_preset_drain.py ###
from __future__ import annotations

import simpy
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph   import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import (
    Create,
    AddObjectProperty,
    UnitaryEditType,
)
from libsyn_tools.sim.operation_preset.drain import DrainExcess


def _mk_src_dst_with_poms(volumes):
    src = MaterialContainer()
    dst = MaterialContainer()
    poms = []
    for v in volumes:
        p = PortionOfMaterial()
        p.add_chemical(Chemical(mass=v, density=1.0))
        poms.append(p)
    return src, dst, poms


def test_drain_excess_edits_shape(env: simpy.Environment):
    # src has 7 + 5 mL; target_volume = 8 → move 4 mL to dst
    src, dst, poms = _mk_src_dst_with_poms([7.0, 5.0])

    eng = EffectEngine()
    for o in (src, dst, *poms):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
    seed_edits = [Create(instance_1_iri=src.identifier), Create(instance_1_iri=dst.identifier)]
    for p in poms:
        seed_edits += [
            Create(instance_1_iri=p.identifier),
            AddObjectProperty(
                instance_1_iri=p.identifier,
                instance_2_iri=src.identifier,
                property_iri=Is_directly_contained_by.predicate_iri,
            ),
        ]
    eng.apply(seed_edits, env, operation_id="seed",
              locked_iris=[src.identifier, dst.identifier, *[p.identifier for p in poms]])

    op = DrainExcess(
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        target_volume=8.0,
    )
    edits = op.get_operation_effects()
    # We don't assert exact count (depends on split across POMs),
    # but we assert the shape uses only valid primitives and non-empty.
    assert edits
    assert all(isinstance(e.type, UnitaryEditType) for e in edits)
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_preset_drain.py ###
