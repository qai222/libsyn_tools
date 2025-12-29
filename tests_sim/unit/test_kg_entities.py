# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_kg_entities.py ###
from __future__ import annotations

from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph   import (
    PortionOfMaterial,
    MaterialContainer,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import (
    Create,
    AddObjectProperty,
)


def _chem(mass: float, density: float = 1.0) -> Chemical:
    return Chemical(mass=mass, density=density)


def test_pom_volume_and_split_mix():
    pom = PortionOfMaterial()
    pom.add_chemical(_chem(5.0, 1.0))  # 5 mL
    pom.add_chemical(_chem(3.0, 1.0))  # +3 mL = 8 mL total

    assert abs(pom.volume - 8.0) < 1e-9

    # Split 25% / 75%
    p25 = pom.get_portion(0.25)
    p75 = pom.get_portion(0.75)
    assert abs(p25.volume + p75.volume - pom.volume) < 1e-9

    # Mix back equals original total (composition-level test is optional here)
    mix = p25.mix_with(p75)
    assert abs(mix.volume - pom.volume) < 1e-9


def test_container_capacity_property_and_direct_volume(env):
    c = MaterialContainer()
    # Accessing capacity with nothing set should raise
    try:
        _ = c.capacity
        assert False, "expected AttributeError when no capacity set"
    except AttributeError:
        pass

    c.has_capacity.add(50.0)
    assert c.capacity == 50.0

    # direct contained volume via overlay helper path
    eng = EffectEngine()
    KnowledgeGraph.get_object_from_lookup(c.identifier)
    p = PortionOfMaterial()
    p.add_chemical(_chem(12.0, 1.0))
    KnowledgeGraph.get_object_from_lookup(p.identifier)

    # Create both & link POM → container
    eng.apply(
        [
            Create(instance_1_iri=c.identifier),
            Create(instance_1_iri=p.identifier),
            AddObjectProperty(
                instance_1_iri=p.identifier,
                instance_2_iri=c.identifier,
                property_iri=Is_directly_contained_by.predicate_iri,
            ),
        ],
        env,
        operation_id="link",
        locked_iris=[c.identifier, p.identifier],
    )

    # `directly_contained_pom_volume` should match 12.0
    assert abs(c.directly_contained_pom_volume - 12.0) < 1e-9
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_kg_entities.py ###
