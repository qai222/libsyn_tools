from __future__ import annotations

import pytest
import simpy

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim.effect_engine import EffectEngine, EngineMechanicalError
from libsyn_tools.sim.knowledge_graph import Is_directly_contained_by, MaterialContainer, PortionOfMaterial
from libsyn_tools.sim.operation.unitary_edit import ChangeDataProperty, Create, RemoveObjectProperty


def _make_container(*, capacity: float = 10.0) -> MaterialContainer:
    c = MaterialContainer()
    c.has_capacity.add(capacity)
    c.is_present = {True}
    return c


def _make_pom(*, volume: float, container: MaterialContainer) -> PortionOfMaterial:
    chem = Chemical(smiles="O", molecular_weight=18.015, density=1.0, mass=volume)
    pom = PortionOfMaterial()
    pom.add_chemical(chem)
    pom.is_present = {True}
    pom.is_directly_contained_by.add(container)
    return pom


def test_create_unknown_object_is_mechanical_error() -> None:
    env = simpy.Environment()
    engine = EffectEngine(shapes_graph=None)

    with pytest.raises(EngineMechanicalError):
        engine.apply_tx(
            edits=[Create(instance_1_iri="does_not_exist")],
            env=env,
            operation_id="op",
            locked_iris=[],
        )


def test_missing_property_iri_is_mechanical_error() -> None:
    env = simpy.Environment()
    engine = EffectEngine(shapes_graph=None)

    c = _make_container(capacity=10.0)

    # property_iri is required for ChangeDataProperty; ensure we get a friendly mechanical error.
    bad_edit = ChangeDataProperty(instance_1_iri=c.identifier, property_iri=None, data_value=5.0)

    with pytest.raises(EngineMechanicalError):
        engine.apply_tx(edits=[bad_edit], env=env, operation_id="op", locked_iris=[c.identifier])

    assert c.has_capacity == {10.0}


def test_apply_tx_allows_idempotent_remove() -> None:
    """Removing a missing relation should be a no-op while other edits commit."""

    env = simpy.Environment()
    engine = EffectEngine(shapes_graph=None)

    a = _make_container(capacity=10.0)
    b = _make_container(capacity=20.0)
    pom = _make_pom(volume=3.0, container=a)

    edits = [
        # First edit mutates container A.
        ChangeDataProperty(instance_1_iri=a.identifier, property_iri="https://libsyn-sim/kg/has_capacity", data_value=5.0),
        # Second edit fails: attempt to remove containment edge to the wrong container.
        RemoveObjectProperty(
            instance_1_iri=pom.identifier,
            property_iri=Is_directly_contained_by.predicate_iri,
            instance_2_iri=b.identifier,
        ),
    ]

    result = engine.apply_tx(edits=edits, env=env, operation_id="op", locked_iris=[a.identifier, b.identifier])
    assert result.committed is True

    # The first edit should have been committed.
    assert a.has_capacity == {5.0}
    # Containment should remain unchanged.
    assert a in pom.is_directly_contained_by
