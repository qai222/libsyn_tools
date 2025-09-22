# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_effect_engine.py ###
from __future__ import annotations

import simpy
from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.effect_engine import EffectEngine
from libsyn_tools.sim.knowledge_graph.physical_entities import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import (
    Create,
    AddObjectProperty,
    AddDataProperty,
)


def _make_container_present() -> MaterialContainer:
    c = MaterialContainer()
    # present only after Create (use engine later) or set directly here for unit isolation
    return c


def _make_pom_with_volume(vol_ml: float) -> PortionOfMaterial:
    pom = PortionOfMaterial()
    # Build a single-ingredient chemical with known volume
    chem = Chemical(mass=vol_ml, density=1.0)  # 1 g/mL → volume == mass
    pom.add_chemical(chem)
    return pom


def _minimal_overflow_shape() -> Graph:
    """
    SHACL shape: any subject with lib:currentVolume > 9.9 triggers a violation.

    We use a SPARQLConstraint over the **overlay** predicate `lib:currentVolume`.
    """

    lib = Namespace("https://libsyn-sim/kg/")
    ttl = f'''
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    [] a sh:NodeShape ;
       sh:targetSubjectsOf lib:currentVolume ;
       sh:sparql [
         a sh:SPARQLConstraint ;
         sh:select """
    SELECT ?this WHERE {{ ?this lib:currentVolume ?v . FILTER(xsd:double(?v) > 9.9) }}
    """ ;
       ] .
    '''
    return Graph().parse(data=ttl, format="turtle")


def test_mechanical_abort_create_on_present(env: simpy.Environment):
    eng = EffectEngine()
    c = _make_container_present()
    # Mark present first via Create, then attempt another Create in a new batch
    # First: make the object exist in KG and present
    KnowledgeGraph.get_object_from_lookup(c.identifier)  # ensure in KG
    eng.apply([Create(instance_1_iri=c.identifier)], env, operation_id="init", locked_iris=[c.identifier])

    # Now: a second Create on a present object must fail (mechanical)
    try:
        eng.apply([Create(instance_1_iri=c.identifier)], env, operation_id="dup", locked_iris=[c.identifier])
        assert False, "Expected mechanical pre-check to abort on CREATE of present"
    except RuntimeError as e:
        assert "CREATE on present object" in str(e)


def test_mechanical_abort_dangling_subject(env: simpy.Environment):
    eng = EffectEngine()
    fake = "LabObject_FAKE"
    # Add data property to non-existent subject (not created in batch) → abort
    from libsyn_tools.sim.knowledge_graph.physical_entities import Has_interrupt_events
    edit = AddDataProperty(
        instance_1_iri=fake,
        property_iri=Has_interrupt_events.predicate_iri,
        data_value="x",
    )
    try:
        eng.apply([edit], env, operation_id="dangling", locked_iris=[])
        assert False, "Expected mechanical pre-check to catch dangling subject"
    except RuntimeError as e:
        assert "dangling subject" in str(e)


def test_apply_addobjectproperty_and_overlay_shacl(env: simpy.Environment):
    # Build container & POM; apply edits to present them and relate via direct containment
    eng = EffectEngine(shapes_graph=_minimal_overflow_shape())

    c = _make_container_present()
    p = _make_pom_with_volume(12.0)  # > 9.9 → should violate overlay shape

    # Ensure both instances are known to the KG
    KnowledgeGraph.get_object_from_lookup(c.identifier)
    KnowledgeGraph.get_object_from_lookup(p.identifier)

    # Prepare edits: create both + link POM → container
    edits = [
        Create(instance_1_iri=c.identifier),
        Create(instance_1_iri=p.identifier),
        AddObjectProperty(
            instance_1_iri=p.identifier,
            instance_2_iri=c.identifier,
            property_iri=Is_directly_contained_by.predicate_iri,
        ),
    ]
    eng.apply(edits, env, operation_id="link", locked_iris=[c.identifier, p.identifier])

    # Overlay must compute currentVolume for container; SHACL should record violation
    # Verify an overlay triple exists via validate_now()
    conforms, report, _ = eng.validate_now()
    assert not conforms  # should still be non-conformant with same state

    # The engine recorded violations during .apply(...)
    viols = eng._shacl_violations  # for unit tests we access internal buffer
    assert len(viols) >= 1
    v = viols[-1]
    assert v.origin == "SHACL" and v.severity == "soft" and v.disposition == "committed"
    assert v.operation_id == "link"
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_effect_engine.py ###
