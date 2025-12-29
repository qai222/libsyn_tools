# === REPLACE FULL FILE tests_sim/unit/test_effect_engine.py ===
from __future__ import annotations

import simpy
from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
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
    AddDataProperty,
)
from libsyn_tools.sim.overlay.current_volume_overlay import CurrentVolumeOverlayProvider


def _make_container_present() -> MaterialContainer:
    c = MaterialContainer()
    return c


def _make_pom_with_volume(vol_ml: float) -> PortionOfMaterial:
    pom = PortionOfMaterial()
    chem = Chemical(mass=vol_ml, density=1.0)
    pom.add_chemical(chem)
    return pom


def _minimal_overflow_shape() -> Graph:
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
    KnowledgeGraph.get_object_from_lookup(c.identifier)
    eng.apply([Create(instance_1_iri=c.identifier)], env, operation_id="init", locked_iris=[c.identifier])

    try:
        eng.apply([Create(instance_1_iri=c.identifier)], env, operation_id="dup", locked_iris=[c.identifier])
        assert False, "Expected mechanical pre-check to abort on CREATE of present"
    except RuntimeError as e:
        assert "CREATE on present object" in str(e)


def test_mechanical_abort_dangling_subject(env: simpy.Environment):
    eng = EffectEngine()
    fake = "LabObject_FAKE"
    from libsyn_tools.sim.knowledge_graph   import Has_interrupt_events
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
    # Engine with shape + explicit currentVolume overlay provider
    eng = EffectEngine(shapes_graph=_minimal_overflow_shape())
    eng.register_overlay_provider(CurrentVolumeOverlayProvider().snapshot)

    c = _make_container_present()
    p = _make_pom_with_volume(12.0)  # > 9.9 → should violate

    KnowledgeGraph.get_object_from_lookup(c.identifier)
    KnowledgeGraph.get_object_from_lookup(p.identifier)

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

    conforms, report, _ = eng.validate_now()
    assert not conforms
