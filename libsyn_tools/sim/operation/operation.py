from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import Any, Optional, Union

import simpy
from loguru import logger
from pydantic import BaseModel
from simpy.resources.resource import Request
from twa.data_model.base_ontology import BaseClass

from libsyn_tools.sim.knowledge_graph import SimOntology, Field, str_uuid
from libsyn_tools.sim.operation.selector import Selector, LiteralSelector


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

    model_config = {"frozen": True}

    @classmethod
    def create(cls, iri: str):
        return cls(type=UnitaryEditType.CREATE, instance_1_iri=iri)

    @classmethod
    def annihilate(cls, iri: str):
        return cls(type=UnitaryEditType.ANNIHILATE, instance_1_iri=iri)

    @classmethod
    def change_data_prop(cls, iri: str, data_prop: str, data_value: str):
        return cls(type=UnitaryEditType.CHANGE_DATA_PROPERTY, instance_1_iri=iri, property_iri=data_prop,
                   data_value=data_value)

    @classmethod
    def add_obj_prop(cls, iri1: str, property_iri: str, iri2: str):
        return cls(type=UnitaryEditType.ADD_OBJECT_PROPERTY, instance_1_iri=iri1, instance_2_iri=iri2,
                   property_iri=property_iri)

    @classmethod
    def remove_obj_prop(cls, iri1: str, property_iri: str, iri2: str):
        return cls(type=UnitaryEditType.REMOVE_OBJECT_PROPERTY, instance_1_iri=iri1, instance_2_iri=iri2,
                   property_iri=property_iri)

    def apply(self):
        logger.debug(f"applying edit: {self.type}")
        # TODO It is probably better to just use RDFlib
        # TODO type check lab objects
        instance_1 = BaseClass.object_lookup[self.instance_1_iri]

        if self.type == UnitaryEditType.CREATE:
            instance_1.is_present = {True, }

        elif self.type == UnitaryEditType.ANNIHILATE:
            instance_1.is_present = {False, }

        elif self.type == UnitaryEditType.CHANGE_DATA_PROPERTY:
            assert instance_1.is_present == {True, }, (f"changing data property='{self.property_iri}' for a lab "
                                                       f"object='{self.instance_1_iri}' that is absent")
            data_property = SimOntology.data_property_lookup[self.property_iri]
            # TODO this only applies to functional data property
            data_property_name = data_property.__class__.__name__
            field_name = data_property_name[0].lower() + data_property_name[1:]
            setattr(instance_1, field_name, {self.data_value, })

        elif self.type in (UnitaryEditType.ADD_OBJECT_PROPERTY, UnitaryEditType.REMOVE_OBJECT_PROPERTY):
            assert instance_1.is_present == {
                True, }, "changing object property of a lab object but the subject is absent"
            instance_2 = BaseClass.object_lookup[self.instance_2_iri]
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


class Presumption(BaseModel):
    """
    A generic presumption (precondition) that must be satisfied
    before an Action can be executed safely .
    """

    @abstractmethod
    def validate(self, knowledge_graph: Any) -> bool:
        """
        Checks whether this presumption is satisfied by the current state
        of the knowledge graph. Returns True if satisfied, or False otherwise.

        In your use-case, you might raise a ValidationError instead of
        returning False.
        """
        pass


def _collect_participant_specs(operation: "Operation") -> dict[str, StrOrSelector]:
    """
    Scan all fields that start with 'participant_' and build a mapping
    role → spec (str | Selector).
    """
    return {
        fname[len("participant_"):]: getattr(operation, fname)
        for fname in operation.model_fields
        if fname.startswith("participant_")
    }


def _write_participant_iris(operation: "Operation", resolved: dict[str, str]):
    """Overwrite participant_* fields with literal IRIs (in-place)."""
    for role, iri in resolved.items():
        setattr(operation, f"participant_{role}", iri)


StrOrSelector = Union[str, Selector]


class Operation(ABC, BaseModel):
    """
    A `Operation` is a process that changes the knowledge graph.

    Concrete operation define `participant_<role>` attrs and implement `get_action_effects()` to return a list of
    UnitaryEdits built *after* all roles are resolved.
    """

    identifier: str = Field(default_factory=str_uuid)
    """ identifier of this operation """

    temporal_cost: Optional[float] = None
    # TODO this may depend on the actual knowledge graph right before its execution
    """ an estimate of how long this operation would take """

    scheduled_start_time: Optional[float] = None
    """ 
    the scheduled start time, the actual start time in a simulation of this operation cannot be earlier than the 
    scheduled start time 
    """

    required_precedents: list[str] = Field(default_factory=list)
    """ the uuids of the required precedent operations that must precede this operation """

    presumptions: list[Presumption] = Field(default_factory=list)
    # TODO formalize and implement
    # TODO we could define functions to validate presumptions in subclasses,
    #  or we can use SHACL like in https://github.com/RDFLib/pySHACL
    """
    a set of assumptions of the world that serve as the prerequisites for this operation to be executed
    example: the robot arm is not occupied by any other operations
    example: the container should contain at least 10 mL liquid
    """

    operation_effects: list[UnitaryEdit] = Field(default_factory=list)
    """ 
    a list of unitary graph edits to the knowledge graph 
    """
    # TODO there may be a hard-to-model difference between "transfer 10 mL from A to B" and
    #  "transfer half of what is inside A to B": the latter depends on the state of A right before execution and it is
    #  almost implied that during this transfer the content of A does not change.

    operation_effects_description: Optional[str] = None
    """ free text description for the effects of this action """

    resources: list[str] = Field(default_factory=list)
    """
    a list of uuids of the lab objects that will be occupied during the execution of this operation, used in DES as 
    `resources`.
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
    #
    # Cullen: at the high-level wait for smth is always better than causing problems. disable general and enable specific.

    _locks: list[Request] = Field(default_factory=list, exclude=True)
    """
    runtime-only attributes (excluded from serialisation)
    """

    def execute(self):
        """ applying operation effects """
        logger.info(f"execute action: {self.identifier}")
        for edit in self.operation_effects:
            edit.apply()

    @abstractmethod
    def get_operation_effects(self) -> list[UnitaryEdit]:
        pass

    def pre_act(self, env: simpy.Environment):
        """
        * resolve selectors  → IRIs
        * acquire locks      (capacity-1 per LabObject)
        * populate resources (resolved IRIs)
        * compute action_effects
        Returns a SimPy Event so the scheduler can `yield` on it.
        """
        return env.process(self._pre_act_implementation(env))

    def _pre_act_implementation(self, env: simpy.Environment) -> None:
        participant_specs = _collect_participant_specs(self)

        # deterministic ordering prevents dead-locks -------------------
        ordered_specs = sorted(participant_specs.items(), key=lambda kv: str(kv[1]))

        resolved: dict[str, str] = {}

        for role, spec in ordered_specs:
            if isinstance(spec, Selector):
                iri, req = yield env.process(spec.resolve(env))
            elif isinstance(spec, str):
                # lock via LiteralSelector to keep path uniform
                iri, req = yield env.process(
                    LiteralSelector(spec).resolve(env)
                )
            else:
                raise TypeError(
                    f"Participant '{role}' has unsupported type {type(spec)}"
                )
            resolved[role] = iri
            self._locks.append(req)

        # overwrite participant_* fields with pure strings -------------
        _write_participant_iris(self, resolved)

        # also expose them via resources[] for backward compatibility
        self.resources = list(resolved.values())

        # build list of graph edits now that everything is bound -------
        self.operation_effects = self.get_operation_effects()

    def post_act(self):
        for req in self._locks:
            req.resource.release(req)
        self._locks.clear()

    class Config:
        arbitrary_types_allowed = True
        validate_assignment = False
        extra = "forbid"


Operation.model_rebuild()
UnitaryEdit.model_rebuild()
Presumption.model_rebuild()
