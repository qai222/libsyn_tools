from __future__ import annotations

import pytest
import simpy

from libsyn_tools.sim.effect_engine import EffectEngine, EngineMechanicalError
from libsyn_tools.sim.knowledge_graph import MaterialContainer, PortionOfMaterial, Has_capacity
from libsyn_tools.sim.operation.unitary_edit import ChangeDataProperty


def test_unknown_property_iri_is_mechanical_error() -> None:
    """An unknown property_iri should fail as an EngineMechanicalError.

    This is important because snapshotting uses ontology lookups; without a
    precheck, an unknown property IRI can raise a raw KeyError and bypass the
    normal transaction rollback/error reporting.
    """

    env = simpy.Environment()
    eng = EffectEngine(shapes_graph=None)

    c = MaterialContainer()
    c.has_capacity.add(10.0)
    c.is_present = {True}

    bad = ChangeDataProperty(
        instance_1_iri=c.identifier,
        property_iri="https://libsyn-sim/kg/not_a_real_property",
        data_value=5.0,
    )

    with pytest.raises(EngineMechanicalError):
        eng.apply_tx(edits=[bad], env=env, operation_id="op", locked_iris=[c.identifier])

    assert c.has_capacity == {10.0}


def test_subject_missing_property_field_is_mechanical_error() -> None:
    """A property IRI may exist in the ontology but not be valid on a given
    subject type; catch this mechanically rather than via AttributeError.
    """

    env = simpy.Environment()
    eng = EffectEngine(shapes_graph=None)

    pom = PortionOfMaterial()
    pom.is_present = {True}

    # Has_capacity is a valid data property, but PortionOfMaterial doesn't declare
    # a `has_capacity` field.
    bad = ChangeDataProperty(
        instance_1_iri=pom.identifier,
        property_iri=Has_capacity.predicate_iri,
        data_value=5.0,
    )

    with pytest.raises(EngineMechanicalError):
        eng.apply_tx(edits=[bad], env=env, operation_id="op", locked_iris=[])