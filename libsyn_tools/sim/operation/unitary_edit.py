from __future__ import annotations

from abc import abstractmethod
from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph import SimOntology


class UnitaryEditType(str, Enum):
    """ possible types of a unitary edit """

    CREATE = "CREATE"
    """ create a lab object """

    ANNIHILATE = "ANNIHILATE"
    """ annihilate a lab object """

    CHANGE_DATA_PROPERTY = "CHANGE_DATA_PROPERTY"
    """ change a data property of a lab object """

    ADD_DATA_PROPERTY = "ADD_DATA_PROPERTY"
    """ add data to non-functional data property """

    REMOVE_DATA_PROPERTY = "REMOVE_DATA_PROPERTY"
    """ remove data from non-functional data property """

    ADD_OBJECT_PROPERTY = "ADD_OBJECT_PROPERTY"
    """ add an object property between two lab objects """

    REMOVE_OBJECT_PROPERTY = "REMOVE_OBJECT_PROPERTY"
    """ remove an object property between two lab objects """


class UnitaryEdit(BaseModel):
    """A single, atomic mutation of the knowledge graph."""

    type: UnitaryEditType
    instance_1_iri: str
    instance_2_iri: str | None = None
    property_iri: str | None = None
    data_value: Any | None = None

    # Not frozen: tests and debug tooling may monkeypatch `apply` on instances
    # to simulate unexpected failures. The engine treats edits as immutable by
    # convention, but allowing mutation here improves testability.
    # Allow setting extra attributes on instances (e.g. monkeypatching `apply` in tests)
    # while keeping edits conceptually immutable by convention.
    model_config = {"frozen": False, "extra": "allow"}

    @abstractmethod
    def apply(self) -> None: ...

    def compute_inverse(self) -> "UnitaryEdit": ...  # implemented per subclass

    def describe(self) -> str:
        if self.type is UnitaryEditType.CREATE:
            return f"CREATE {self.instance_1_iri}"
        if self.type is UnitaryEditType.ANNIHILATE:
            return f"ANNIHILATE {self.instance_1_iri}"
        if self.type is UnitaryEditType.CHANGE_DATA_PROPERTY:
            return f"SET {self.instance_1_iri} {self.property_iri} = {self.data_value!r}"
        if self.type is UnitaryEditType.ADD_DATA_PROPERTY:
            return f"ADD_DATA {self.instance_1_iri} {self.property_iri} += {self.data_value!r}"
        if self.type is UnitaryEditType.REMOVE_DATA_PROPERTY:
            return f"REMOVE_DATA {self.instance_1_iri} {self.property_iri} -= {self.data_value!r}"
        if self.type is UnitaryEditType.ADD_OBJECT_PROPERTY:
            return f"ADD {self.instance_1_iri} {self.property_iri} {self.instance_2_iri}"
        if self.type is UnitaryEditType.REMOVE_OBJECT_PROPERTY:
            return f"REMOVE {self.instance_1_iri} {self.property_iri} {self.instance_2_iri}"
        return f"{self.type.value} {self.instance_1_iri}"

    # backwards-compat aliases --------------------------------------
    @classmethod
    def create(cls, iri: str):  # noqa: N802
        return Create(instance_1_iri=iri)

    @classmethod
    def annihilate(cls, iri: str):  # noqa: N802
        return Annihilate(instance_1_iri=iri)

    @classmethod
    def change_data_prop(cls, iri: str, data_prop: str, data_value: Any):
        return ChangeDataProperty(instance_1_iri=iri,
                                  property_iri=data_prop,
                                  data_value=data_value)

    @classmethod
    def add_obj_prop(cls, iri1: str, prop: str, iri2: str):
        return AddObjectProperty(instance_1_iri=iri1,
                                 instance_2_iri=iri2,
                                 property_iri=prop)

    @classmethod
    def remove_obj_prop(cls, iri1: str, prop: str, iri2: str):
        return RemoveObjectProperty(instance_1_iri=iri1,
                                    instance_2_iri=iri2,
                                    property_iri=prop)

    def __init_subclass__(cls, **kwargs):
        """
        The five edit primitives declared below (Create, Annihilate, ChangeDataProperty, AddObjectProperty,
        RemoveObjectProperty) form a complete set of atomic mutations for the simulator’s RDF-like knowledge graph.
        Higher-level changes should be expressed as *sequences* of these primitives.
        For that reason further subclassing of `UnitaryEdit` is explicitly blocked.
        """
        if cls.__module__ != __name__:
            raise TypeError(
                "Sub-classing UnitaryEdit outside operation.unitary_edit "
                "is not allowed.  Use a sequence of the existing "
                "edit primitives instead."
            )
        super().__init_subclass__(**kwargs)


def _require_object(iri: str, *, role: str, edit: UnitaryEdit) -> Any:
    obj = KnowledgeGraph.get_object_from_lookup(iri)
    if obj is None:
        raise RuntimeError(
            f"{edit.type.value}: unknown {role} object {iri!r}"
        )
    return obj


def _require_data_property(edit: UnitaryEdit) -> Any:
    if not edit.property_iri:
        raise RuntimeError(f"{edit.type.value}: property_iri is required")
    prop = SimOntology.data_property_lookup.get(edit.property_iri)
    if prop is None:
        raise RuntimeError(
            f"{edit.type.value}: unknown data property {edit.property_iri!r}"
        )
    return prop


def _require_object_property(edit: UnitaryEdit) -> Any:
    if not edit.property_iri:
        raise RuntimeError(f"{edit.type.value}: property_iri is required")
    prop = SimOntology.object_property_lookup.get(edit.property_iri)
    if prop is None:
        raise RuntimeError(
            f"{edit.type.value}: unknown object property {edit.property_iri!r}"
        )
    return prop


class Create(UnitaryEdit):
    type: Literal[UnitaryEditType.CREATE] = UnitaryEditType.CREATE

    def apply(self):
        _require_object(self.instance_1_iri, role="subject", edit=self).is_present = {True}


class Annihilate(UnitaryEdit):
    type: Literal[UnitaryEditType.ANNIHILATE] = UnitaryEditType.ANNIHILATE

    def apply(self):
        _require_object(self.instance_1_iri, role="subject", edit=self).is_present = {False}


class AddDataProperty(UnitaryEdit):
    type: Literal[UnitaryEditType.ADD_DATA_PROPERTY] = UnitaryEditType.ADD_DATA_PROPERTY

    def apply(self):
        subj = _require_object(self.instance_1_iri, role="subject", edit=self)
        data_prop = _require_data_property(self)
        field_name = data_prop.__name__[0].lower() + data_prop.__name__[1:]
        existing_values = getattr(subj, field_name)
        existing_values.add(self.data_value)


class RemoveDataProperty(UnitaryEdit):
    type: Literal[UnitaryEditType.REMOVE_DATA_PROPERTY] = UnitaryEditType.REMOVE_DATA_PROPERTY

    def apply(self):
        subj = _require_object(self.instance_1_iri, role="subject", edit=self)
        data_prop = _require_data_property(self)
        field_name = data_prop.__name__[0].lower() + data_prop.__name__[1:]
        existing_values = getattr(subj, field_name)
        if self.data_value in existing_values:
            existing_values.remove(self.data_value)


class ChangeDataProperty(UnitaryEdit):
    type: Literal[UnitaryEditType.CHANGE_DATA_PROPERTY] = UnitaryEditType.CHANGE_DATA_PROPERTY

    def apply(self):
        subj = _require_object(self.instance_1_iri, role="subject", edit=self)
        data_prop = _require_data_property(self)
        field_name = data_prop.__name__[0].lower() + data_prop.__name__[1:]
        setattr(subj, field_name, {self.data_value})


class AddObjectProperty(UnitaryEdit):
    type: Literal[UnitaryEditType.ADD_OBJECT_PROPERTY] = (
        UnitaryEditType.ADD_OBJECT_PROPERTY
    )

    def apply(self):
        if not self.instance_2_iri:
            raise RuntimeError(f"{self.type.value}: instance_2_iri is required")
        subj = _require_object(self.instance_1_iri, role="subject", edit=self)
        obj = _require_object(self.instance_2_iri, role="object", edit=self)
        obj_prop = _require_object_property(self)
        field_name = obj_prop.__name__[0].lower() + obj_prop.__name__[1:]
        getattr(subj, field_name).add(obj)


class RemoveObjectProperty(UnitaryEdit):
    type: Literal[UnitaryEditType.REMOVE_OBJECT_PROPERTY] = (
        UnitaryEditType.REMOVE_OBJECT_PROPERTY
    )

    def apply(self):
        if not self.instance_2_iri:
            raise RuntimeError(f"{self.type.value}: instance_2_iri is required")
        subj = _require_object(self.instance_1_iri, role="subject", edit=self)
        obj = _require_object(self.instance_2_iri, role="object", edit=self)
        obj_prop = _require_object_property(self)
        field_name = obj_prop.__name__[0].lower() + obj_prop.__name__[1:]
        existing_values = getattr(subj, field_name)
        if obj in existing_values:
            existing_values.remove(obj)


# update exported symbols so `from ... import *` keeps working
__all__ = [
    "UnitaryEditType",
    "UnitaryEdit",
    "Create",
    "Annihilate",
    "ChangeDataProperty",
    "AddObjectProperty",
    "RemoveObjectProperty",
    "AddDataProperty",
    "RemoveDataProperty",
]
