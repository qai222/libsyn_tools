# === REPLACE FULL FILE tests_sim/unit/test_effect_engine.py ===
from __future__ import annotations

import simpy
from rdflib import Graph, Namespace
from rdflib.namespace import SH, XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.effect_engine import EffectEngine, EngineMechanicalError
from libsyn_tools.sim.policy import PolicyBundle, PolicyRule
from libsyn_tools.sim.knowledge_graph   import (
    Has_interrupt_events,
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import (
    Create,
    AddObjectProperty,
    AddDataProperty,
    Annihilate,
    ChangeDataProperty,
)
from libsyn_tools.sim.overlay.current_volume_overlay import CurrentVolumeOverlayProvider
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize


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


def _named_overflow_shape() -> tuple[Graph, str]:
    lib = Namespace("https://libsyn-sim/kg/")
    shape_iri = str(lib.OverflowShape)
    ttl = f'''
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    lib:OverflowShape a sh:NodeShape ;
       sh:targetSubjectsOf lib:currentVolume ;
       sh:sparql [
         a sh:SPARQLConstraint ;
         sh:select """
    SELECT ?this WHERE {{ ?this lib:currentVolume ?v . FILTER(xsd:double(?v) > 9.9) }}
    """ ;
       ] .
    '''
    return Graph().parse(data=ttl, format="turtle"), shape_iri


def _capacity_overflow_shape() -> tuple[Graph, str]:
    lib = Namespace("https://libsyn-sim/kg/")
    shape_iri = str(lib.CapacityOverflowShape)
    ttl = f'''
    PREFIX sh: <{SH}>
    PREFIX xsd: <{XSD}>
    PREFIX lib: <{lib}>
    lib:CapacityOverflowShape a sh:NodeShape ;
       sh:targetSubjectsOf lib:currentVolume ;
       sh:sparql [
         a sh:SPARQLConstraint ;
         sh:select """
    SELECT ?this WHERE {{
      ?this lib:currentVolume ?v ;
            lib:has_capacity ?cap .
      FILTER(xsd:double(?v) > xsd:double(?cap))
    }}
    """ ;
       ] .
    '''
    return Graph().parse(data=ttl, format="turtle"), shape_iri


def _world_reservoir(amount: float):
    base = "https://libsyn-sim/kg/"
    v1 = MaterialContainer(identifier=f"{base}v1")
    res = MaterialContainer(identifier=f"{base}reservoir")
    pip = MaterialContainer(identifier=f"{base}pipette")
    p = PortionOfMaterial(identifier=f"{base}water")
    p.add_chemical(Chemical(mass=amount, density=1.0))
    for o in (v1, res, pip, p):
        KnowledgeGraph.get_object_from_lookup(o.identifier)
        Create(instance_1_iri=o.identifier).apply()
    AddObjectProperty(
        instance_1_iri=p.identifier,
        instance_2_iri=res.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()
    return v1, res, pip, p


def test_mechanical_abort_create_on_present(env: simpy.Environment):
    eng = EffectEngine()
    c = _make_container_present()
    KnowledgeGraph.get_object_from_lookup(c.identifier)
    eng.apply([Create(instance_1_iri=c.identifier)], env, operation_id="init", locked_iris=[c.identifier])

    try:
        eng.apply([Create(instance_1_iri=c.identifier)], env, operation_id="dup", locked_iris=[c.identifier])
        assert False, "Expected mechanical pre-check to abort on CREATE of present"
    except EngineMechanicalError as e:
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
    except EngineMechanicalError as e:
        assert "dangling subject" in str(e)


def test_mechanical_abort_change_on_non_present(env: simpy.Environment):
    eng = EffectEngine()
    obj = MaterialContainer()
    KnowledgeGraph.get_object_from_lookup(obj.identifier)

    eng.apply([Create(instance_1_iri=obj.identifier)], env, operation_id="init", locked_iris=[obj.identifier])
    eng.apply([Annihilate(instance_1_iri=obj.identifier)], env, operation_id="kill", locked_iris=[obj.identifier])

    edit = ChangeDataProperty(
        instance_1_iri=obj.identifier,
        property_iri=Has_interrupt_events.predicate_iri,
        data_value="after",
    )
    try:
        eng.apply([edit], env, operation_id="mutate", locked_iris=[obj.identifier])
        assert False, "Expected mechanical pre-check to abort on non-present object"
    except EngineMechanicalError as e:
        assert "non-present subject" in str(e)


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


def test_shacl_policy_applied_to_violation(env: simpy.Environment):
    shape_graph, shape_iri = _named_overflow_shape()
    policy = PolicyBundle(per_shape={shape_iri: PolicyRule(severity="hard", disposition="aborted")})
    eng = EffectEngine(shapes_graph=shape_graph, policy=policy)
    eng.register_overlay_provider(CurrentVolumeOverlayProvider().snapshot)

    c = _make_container_present()
    p = _make_pom_with_volume(12.0)
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
    result = eng.apply_tx(edits, env, operation_id="link", locked_iris=[c.identifier, p.identifier])

    viols = result.violations
    assert any(
        v.origin == "SHACL"
        and v.shape_iri == shape_iri
        and v.severity == "hard"
        and v.disposition == "aborted"
        for v in viols
    )


def test_apply_tx_returns_result_and_violations(env: simpy.Environment):
    shape_graph, shape_iri = _named_overflow_shape()
    eng = EffectEngine(shapes_graph=shape_graph)
    eng.register_overlay_provider(CurrentVolumeOverlayProvider().snapshot)

    c = _make_container_present()
    p = _make_pom_with_volume(12.0)
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
    result = eng.apply_tx(edits, env, operation_id="link", locked_iris=[c.identifier, p.identifier])

    assert result.batch_id
    assert result.edit_fingerprints
    assert all(v.batch_id == result.batch_id for v in result.violations)
    assert any(v.shape_iri == shape_iri for v in result.violations)
    assert [v.model_dump() for v in result.violations] == [
        v.model_dump() for v in eng._shacl_violations
    ]


def test_apply_tx_rolls_back_on_aborted_violation(env: simpy.Environment):
    shape_graph, shape_iri = _capacity_overflow_shape()
    policy = PolicyBundle(per_shape={shape_iri: PolicyRule(severity="hard", disposition="aborted")})
    eng = EffectEngine(shapes_graph=shape_graph, policy=policy)
    eng.register_overlay_provider(CurrentVolumeOverlayProvider().snapshot)

    v1, res, pip, pom = _world_reservoir(12.0)
    v1.has_capacity.add(5.0)

    op = TransferMaterialByPortionSize(
        identifier="fill",
        participant_source=res.identifier,
        participant_destination=v1.identifier,
        participant_device=pip.identifier,
        portion_size=1.0,
        temporal_cost=0.0,
    )
    edits = op.get_operation_effects()

    before_res_volume = res.directly_contained_pom_volume
    before_v1_volume = v1.directly_contained_pom_volume
    before_pom_present = set(pom.is_present)
    before_pom_containers = set(pom.is_directly_contained_by)

    result = eng.apply_tx(edits, env, operation_id="fill", locked_iris=[res.identifier, v1.identifier, pip.identifier])

    assert result.committed is False
    assert any(v.disposition == "aborted" for v in result.violations)
    assert res.directly_contained_pom_volume == before_res_volume
    assert v1.directly_contained_pom_volume == before_v1_volume
    assert pom.is_present == before_pom_present
    assert pom.is_directly_contained_by == before_pom_containers
