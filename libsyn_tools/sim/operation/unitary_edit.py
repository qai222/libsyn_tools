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

    model_config = {"frozen": True}

    @abstractmethod
    def apply(self) -> None: ...

    def compute_inverse(self) -> "UnitaryEdit": ...  # implemented per subclass

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


class Create(UnitaryEdit):
    type: Literal[UnitaryEditType.CREATE] = UnitaryEditType.CREATE

    def apply(self):
        KnowledgeGraph.get_object_from_lookup(self.instance_1_iri).is_present = {True}


class Annihilate(UnitaryEdit):
    type: Literal[UnitaryEditType.ANNIHILATE] = UnitaryEditType.ANNIHILATE

    def apply(self):
        KnowledgeGraph.get_object_from_lookup(self.instance_1_iri).is_present = {False}


class AddDataProperty(UnitaryEdit):
    type: Literal[UnitaryEditType.ADD_DATA_PROPERTY] = UnitaryEditType.ADD_DATA_PROPERTY

    def apply(self):
        subj = KnowledgeGraph.get_object_from_lookup(self.instance_1_iri)
        data_prop = SimOntology.data_property_lookup[self.property_iri]
        field_name = data_prop.__name__[0].lower() + data_prop.__name__[1:]
        existing_values = getattr(subj, field_name)
        existing_values.add(self.data_value)


class ChangeDataProperty(UnitaryEdit):
    type: Literal[UnitaryEditType.CHANGE_DATA_PROPERTY] = UnitaryEditType.CHANGE_DATA_PROPERTY

    def apply(self):
        subj = KnowledgeGraph.get_object_from_lookup(self.instance_1_iri)
        data_prop = SimOntology.data_property_lookup[self.property_iri]
        field_name = data_prop.__name__[0].lower() + data_prop.__name__[1:]
        setattr(subj, field_name, {self.data_value})


class AddObjectProperty(UnitaryEdit):
    type: Literal[UnitaryEditType.ADD_OBJECT_PROPERTY] = (
        UnitaryEditType.ADD_OBJECT_PROPERTY
    )

    def apply(self):
        subj = KnowledgeGraph.get_object_from_lookup(self.instance_1_iri)
        obj = KnowledgeGraph.get_object_from_lookup(self.instance_2_iri)
        obj_prop = SimOntology.object_property_lookup[self.property_iri]
        field_name = obj_prop.__name__[0].lower() + obj_prop.__name__[1:]
        getattr(subj, field_name).add(obj)


class RemoveObjectProperty(UnitaryEdit):
    type: Literal[UnitaryEditType.REMOVE_OBJECT_PROPERTY] = (
        UnitaryEditType.REMOVE_OBJECT_PROPERTY
    )

    def apply(self):
        subj = KnowledgeGraph.get_object_from_lookup(self.instance_1_iri)
        obj = KnowledgeGraph.get_object_from_lookup(self.instance_2_iri)
        obj_prop = SimOntology.object_property_lookup[self.property_iri]
        field_name = obj_prop.__name__[0].lower() + obj_prop.__name__[1:]
        getattr(subj, field_name).remove(obj)


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
]
