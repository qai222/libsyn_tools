from __future__ import annotations

import pytest

from libsyn_tools.sim.knowledge_graph import (
    Has_capacity,
    Is_directly_contained_by,
    MaterialContainer,
    PortionOfMaterial,
)
from libsyn_tools.sim.operation.unitary_edit import (
    AddObjectProperty,
    Annihilate,
    AddDataProperty,
    ChangeDataProperty,
    Create,
    RemoveDataProperty,
    RemoveObjectProperty,
)


@pytest.mark.parametrize(
    "edit, expected",
    [
        (Create(instance_1_iri="missing-subject"), "CREATE: unknown subject object"),
        (Annihilate(instance_1_iri="missing-subject"), "ANNIHILATE: unknown subject object"),
        (
            ChangeDataProperty(
                instance_1_iri="missing-subject",
                property_iri=Has_capacity.predicate_iri,
                data_value=1.0,
            ),
            "CHANGE_DATA_PROPERTY: unknown subject object",
        ),
        (
            AddDataProperty(
                instance_1_iri="missing-subject",
                property_iri=Has_capacity.predicate_iri,
                data_value=1.0,
            ),
            "ADD_DATA_PROPERTY: unknown subject object",
        ),
        (
            RemoveDataProperty(
                instance_1_iri="missing-subject",
                property_iri=Has_capacity.predicate_iri,
                data_value=1.0,
            ),
            "REMOVE_DATA_PROPERTY: unknown subject object",
        ),
    ],
)
def test_direct_unitary_edits_raise_clear_error_for_unknown_subject(edit, expected: str) -> None:
    with pytest.raises(RuntimeError, match=expected):
        edit.apply()


def test_add_object_property_raises_clear_error_for_unknown_object() -> None:
    src = PortionOfMaterial(identifier="https://libsyn-sim/kg/src_for_obj_missing")
    src.is_present = {True}

    with pytest.raises(RuntimeError, match="ADD_OBJECT_PROPERTY: unknown object object"):
        AddObjectProperty(
            instance_1_iri=src.identifier,
            instance_2_iri="https://libsyn-sim/kg/missing_target",
            property_iri=Is_directly_contained_by.predicate_iri,
        ).apply()


def test_remove_object_property_raises_clear_error_for_unknown_subject_and_object() -> None:
    dst = MaterialContainer(identifier="https://libsyn-sim/kg/dst_for_remove_missing")
    dst.is_present = {True}

    with pytest.raises(RuntimeError, match="REMOVE_OBJECT_PROPERTY: unknown subject object"):
        RemoveObjectProperty(
            instance_1_iri="https://libsyn-sim/kg/missing_subject",
            instance_2_iri=dst.identifier,
            property_iri=Is_directly_contained_by.predicate_iri,
        ).apply()

    src = PortionOfMaterial(identifier="https://libsyn-sim/kg/src_for_remove_missing")
    src.is_present = {True}
    with pytest.raises(RuntimeError, match="REMOVE_OBJECT_PROPERTY: unknown object object"):
        RemoveObjectProperty(
            instance_1_iri=src.identifier,
            instance_2_iri="https://libsyn-sim/kg/missing_object",
            property_iri=Is_directly_contained_by.predicate_iri,
        ).apply()
