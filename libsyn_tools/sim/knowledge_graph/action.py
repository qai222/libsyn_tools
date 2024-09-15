from __future__ import annotations

from enum import Enum
from typing import Any
from typing import Optional

from loguru import logger
from pydantic import BaseModel
from twa.data_model.base_ontology import KnowledgeGraph

from .base import SimOntology, Field, str_uuid


# TODO it doesn't seem necessary to use ontology classes for instructions
# TODO classes related to instructions are close to `information entity` that may deserve an abs class


class UnitaryEditType(str, Enum):
    """ possible types of a unitary edit """

    CREATE = "CREATE"
    """ create a lab object """

    ANNIHILATE = "ANNIHILATE"
    """ annihilate a lab object """

    CHANGE_DATA_PROPERTY = "CHANGE_DATA_PROPERTY"
    """ change a data property of a lab object """

    ADD_OBJECT_PROPERTY = "ADD_OBJECT_PROPERTY"
    """ add an object property between two lab objects """

    REMOVE_OBJECT_PROPERTY = "REMOVE_OBJECT_PROPERTY"
    """ remove an object property between two lab objects """


class UnitaryEdit(BaseModel):
    """ a unitary edit is a change to the knowledge graph that cannot be further divided """

    type: UnitaryEditType
    """ type of this unitary edit """

    instance_1_iri: Optional[str] = None
    """ the iri of the first instance, usually the subject """

    instance_2_iri: Optional[str] = None
    """ the iri of the second instance, usually the object """

    property_iri: Optional[str] = None
    """ the iri of the predicate """

    data_value: Optional[Any] = None
    """ data property value if this edit is to change a data property """

    def apply(self):
        logger.info(f"applying edit: {self.type}")
        # TODO It is probably better to just use RDFlib
        # TODO type check lab objects
        instance_1 = KnowledgeGraph.get_object_from_lookup(iri=self.instance_1_iri)

        if self.type == UnitaryEditType.CREATE:
            instance_1.is_present = {True, }

        elif self.type == UnitaryEditType.ANNIHILATE:
            instance_1.is_present = {False, }

        elif self.type == UnitaryEditType.CHANGE_DATA_PROPERTY:
            assert instance_1.is_present == {True, }, "changing data property for a lab object that is absent"
            data_property = SimOntology.data_property_lookup[self.property_iri]
            data_property_name = data_property.__class__.__name__
            field_name = data_property_name[0].lower() + data_property_name[1:]
            setattr(instance_1, field_name, {self.data_value, })

        elif self.type in (UnitaryEditType.ADD_OBJECT_PROPERTY, UnitaryEditType.REMOVE_OBJECT_PROPERTY):
            assert instance_1.is_present == {
                True, }, "changing object property of a lab object but the subject is absent"
            instance_2 = KnowledgeGraph.get_object_from_lookup(iri=self.instance_2_iri)
            assert instance_2.is_present == {
                True, }, "changing object property of a lab object but the object is absent"
            object_property = SimOntology.object_property_lookup[self.property_iri]
            object_property_name = object_property.__class__.__name__
            field_name = object_property_name[0].lower() + object_property_name[1:]
            if self.type == UnitaryEditType.ADD_OBJECT_PROPERTY:
                getattr(instance_1, field_name).add(instance_2)
            else:
                getattr(instance_1, field_name).remove(instance_2)

        else:
            raise ValueError(f"unknown unitary edit type: {self.type}")


class Action(BaseModel):
    """
    An `Action` is a process that changes the knowledge graph
    """
    # TODO make this an abstract class

    identifier: str = Field(default_factory=str_uuid)
    """ identifier of this action """

    temporal_cost: Optional[float] = None
    # TODO this may depend on the actual knowledge graph right before its execution
    """ an estimate of how long this action would take """

    presumptions: Optional[Any] = None
    # TODO formalize and implement
    """  
    a set of assumptions of the world that serve as the prerequisites for this action to be executed
    example: the robot arm is not occupied by any other actions
    example: the container should contain at least 10 mL liquid
    """

    action_effects: list[UnitaryEdit] = []
    """ 
    a list of unitary graph edits to the knowledge graph 
    """
    # TODO there may be a hard-to-model difference between "transfer 10 mL from A to B" and
    #  "transfer half of what is inside A to B": the latter depends on the state of A right before execution and it is
    #  almost implied that during this transfer the content of A does not change.

    action_effects_description: Optional[str] = None
    """ free text description for the effects of this action """

    resources: Optional[list[str]] = []
    """
    a list of uuids of the lab objects that will be occupied during the execution of this action,
    used in DES as `resources`
    """

    # TODO It may make sense to also include specific properties, for example
    #  a beaker is being dried in an oven, during this process, the beaker is limited to be inside this oven, i.e.
    #  `beaker.location` is `occupied` but other properties of `beaker` can be changed, such as adding a stirring
    #  bar to this beaker while it's being dried in the oven.
    #  A better way may be formalizing occupations like this as
    #  statements/restrictions that should be active during the execution of this action, and other actions cannot
    #  violate them during their executions.
    #  This could also be realized in `presumptions`, e.g. if the beaker is being dried
    #  any other action should check if it will change `beaker.location`

    def get_action_effects(self) -> list[UnitaryEdit]:
        pass

    def get_resources(self) -> list[str]:
        pass

    def execute(self):
        """ applying action effects """
        logger.info(f"execute action: {self.identifier}")
        for edit in self.action_effects:
            edit.apply()
