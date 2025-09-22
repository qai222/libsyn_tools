from __future__ import annotations

from abc import ABC, abstractmethod
from enum import StrEnum, auto
from typing import Any, Optional, Union

import simpy
from pydantic import BaseModel, Field
from simpy.resources.resource import Request

from libsyn_tools.sim.operation.runtime import get_object_for_resource
from libsyn_tools.sim.operation.selector import Selector, LiteralSelector, FilterStoreRegistry
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit, UnitaryEditType
from libsyn_tools.utils import str_uuid


class _OpState(StrEnum):
    NEW = auto()  # never prepared
    PREPARED = auto()  # selectors resolved, locks held
    RUNNING = auto()  # currently inside OperationProcess.run()
    FINISHED = auto()


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

    sim_state: _OpState = Field(default=_OpState.NEW, exclude=True)
    """ life cycle tag """

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

    resolved_resources: dict[str, str] = Field(default_factory=dict)
    """
    resolved participant-iri pairs
    """

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

    locks: list[Request] = Field(default_factory=list, exclude=True)
    """
    runtime-only attributes (excluded from serialisation)
    """

    def execute(self):
        """ applying operation effects """
        raise RuntimeError(
            "Direct execute() is removed. Run the operation inside `Simulation` or call EffectEngine.apply(...) explicitly."
        )
        # logger.info(f"execute action: {self.identifier}")
        # for edit in self.operation_effects:
        #     edit.apply()

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
        if self.sim_state is not _OpState.NEW:
            raise RuntimeError(f"{self.identifier}: pre_act called in state {self.sim_state}")
        self.sim_state = _OpState.PREPARED
        return env.process(self._pre_act_implementation(env))

    def _mark_running(self):
        if self.sim_state is not _OpState.PREPARED:
            raise RuntimeError(f"{self.identifier}: run() without successful pre_act")
        self.sim_state = _OpState.RUNNING

    def _pre_act_implementation(self, env: simpy.Environment) -> None:
        participant_specs = _collect_participant_specs(self)

        # deterministic ordering prevents dead-locks -------------------
        ordered_specs = sorted(participant_specs.items(), key=lambda kv: str(kv[1]))

        resolved: dict[str, str] = {}
        acquired: dict[str, simpy.events.Event] = {}  # iri → lock (for dedup)

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

            if iri in acquired:
                req.resource.release(req)  # we already hold the lock
                req = acquired[iri]
            else:
                acquired[iri] = req

            resolved[role] = iri
            self.locks.append(req)

        # overwrite participant_* fields with pure strings -------------
        self.resolved_resources = resolved  # remember bindings
        _write_participant_iris(self, resolved)

        # also expose them via resources[] for backward compatibility
        self.resources = list(resolved.values())

        # build list of graph edits now that everything is bound -------
        self.operation_effects = self.get_operation_effects()

        created_iris = [e.instance_1_iri for e in self.operation_effects
                        if e.type is UnitaryEditType.CREATE]
        if len(created_iris) != len(set(created_iris)):
            raise RuntimeError(
                f"Duplicate CREATE IRIs detected in {self.identifier}: {created_iris}"
            )

    def post_act(self, env: simpy.Environment):
        """
        Release all held locks and reinsert surviving objects into their FilterStores.

        Robust to ANNIHILATE:
        - If an object was annihilated during `apply(...)`, its Resource will have
          been removed from `_RESOURCE_MAP`. Reverse lookups for such resources
          will fail; we still release the lock (on the Resource instance we hold),
          but skip reinsertion (object is not present).
        """
        if self.sim_state is not _OpState.RUNNING:
            raise RuntimeError(f"{self.identifier}: post_act called in state {self.sim_state}")

        for req in self.locks:
            obj = None
            # Try to map Resource -> LabObject; this may fail if the object was annihilated.
            try:
                obj = get_object_for_resource(req.resource)
            except KeyError:
                # Resource no longer registered (likely ANNIHILATE). We can still release the lock below.
                pass

            # Always release the SimPy lock we hold
            req.resource.release(req)

            # Reinsert only if we successfully mapped and the object is still present
            if obj is not None and getattr(obj, "is_present", {False}) == {True}:
                FilterStoreRegistry.put_obj_into_filter_store(obj, env)

        self.locks.clear()
        self.sim_state = _OpState.FINISHED

    class Config:
        arbitrary_types_allowed = True
        validate_assignment = False
        extra = "forbid"


Operation.model_rebuild()
UnitaryEdit.model_rebuild()
Presumption.model_rebuild()
