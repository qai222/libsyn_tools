from __future__ import annotations

from abc import ABC, abstractmethod
from enum import StrEnum, auto
import math
from typing import Any, Optional, Union

import simpy
from pydantic import BaseModel, ConfigDict, Field, field_validator
from simpy.resources.resource import Request

from libsyn_tools.sim.operation.runtime import get_object_for_resource
from libsyn_tools.sim.knowledge_graph import identifier_from_iri
from twa.data_model.base_ontology import KnowledgeGraph
from libsyn_tools.sim.operation.selector import (
    Selector,
    LiteralSelector,
    FilterStoreRegistry,
    SelectorCancelled,
)
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
    # participant_* fields may be optional (None). Skip None so that
    # presets can accept optional devices/resources without breaking pre_act.
    specs: dict[str, StrOrSelector] = {}
    for fname in operation.model_fields:
        if not fname.startswith("participant_"):
            continue
        value = getattr(operation, fname)
        if value is None:
            continue
        specs[fname[len("participant_"):]] = value
    return specs


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
    The scheduled start time expressed in base time. The actual start time in a simulation of this operation
    cannot be earlier than the scheduled start time; the simulation_speed_factor multiplies the base-time delay
    to produce the sim-time timeout.
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

    remediation: bool = Field(default=False)
    """ whether this operation was spawned for remediation """

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
    lock_acquired_times: dict[str, float] = Field(default_factory=dict, exclude=True)

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

    def describe_effects(self) -> list[str]:
        return [edit.describe() for edit in self.operation_effects]

    def pre_act(self, env: simpy.Environment):
        """
        * resolve selectors  → IRIs
        * acquire locks      (capacity-1 per LabObject)
        * populate resources (resolved IRIs)
        * compute action_effects
        Returns a SimPy Event so the scheduler can `yield` on it.
        """
        if self.temporal_cost is not None:
            try:
                is_finite = math.isfinite(self.temporal_cost)
            except TypeError as exc:
                raise ValueError(f"{self.identifier}: temporal_cost must be finite and >= 0") from exc
            if not is_finite or self.temporal_cost < 0:
                raise ValueError(f"{self.identifier}: temporal_cost must be finite and >= 0")
        if self.scheduled_start_time is not None:
            try:
                is_finite = math.isfinite(self.scheduled_start_time)
            except TypeError as exc:
                raise ValueError(f"{self.identifier}: scheduled_start_time must be finite and >= 0") from exc
            if not is_finite or self.scheduled_start_time < 0:
                raise ValueError(f"{self.identifier}: scheduled_start_time must be finite and >= 0")
        if self.sim_state is not _OpState.NEW:
            raise RuntimeError(f"{self.identifier}: pre_act called in state {self.sim_state}")
        self.sim_state = _OpState.PREPARED
        self.lock_acquired_times = {}
        return env.process(self._pre_act_implementation(env))

    def _mark_running(self):
        if self.sim_state is not _OpState.PREPARED:
            raise RuntimeError(f"{self.identifier}: run() without successful pre_act")
        self.sim_state = _OpState.RUNNING

    def _release_all_locks(self, env: simpy.Environment) -> None:
        reinserts: dict[str, Any] = {}
        for req in self.locks:
            obj = None
            try:
                obj = get_object_for_resource(req.resource)
            except KeyError:
                pass
            except Exception:
                obj = None

            resource = getattr(req, "resource", None)
            if resource is not None:
                try:
                    users = getattr(resource, "users", None)
                    queue = getattr(resource, "queue", None)
                    if users is not None and req in users:
                        resource.release(req)
                    elif hasattr(req, "cancel"):
                        req.cancel()
                    elif queue is not None and req in queue:
                        try:
                            queue.remove(req)
                        except ValueError:
                            pass
                except Exception:
                    pass

            if obj is not None and getattr(obj, "is_present", {False}) == {True}:
                reinserts[obj.identifier] = obj
        for obj in reinserts.values():
            FilterStoreRegistry.safe_put_obj_into_filter_store(
                obj,
                env,
                context=f"operation.release_all_locks:{self.identifier}",
            )
        self.locks.clear()

    def _pre_act_implementation(self, env: simpy.Environment) -> None:
        try:
            participant_specs = _collect_participant_specs(self)
    
            # deterministic ordering prevents dead-locks -------------------
            def _ordering_key(role: str, spec: StrOrSelector):
                if isinstance(spec, LiteralSelector):
                    return ("literal", role, "LiteralSelector", identifier_from_iri(spec._iri))
                if isinstance(spec, str):
                    return ("literal", role, "LiteralSelector", identifier_from_iri(spec))
                pool_type = getattr(spec, "pool_type", "")
                return (pool_type, role, spec.__class__.__name__)

            ordered_specs = sorted(participant_specs.items(), key=lambda kv: _ordering_key(kv[0], kv[1]))
    
            resolved: dict[str, str] = {}
            acquired: dict[str, simpy.events.Event] = {}  # iri → lock (for dedup)
            resolve_proc: simpy.events.Process | None = None

    
            for role, spec in ordered_specs:
                if isinstance(spec, (str, LiteralSelector)):
                    if isinstance(spec, LiteralSelector):
                        literal_iri = spec._iri
                    else:
                        literal_iri = spec
                    obj = KnowledgeGraph.get_object_from_lookup(iri=literal_iri)
                    if obj is None:
                        obj = KnowledgeGraph.get_object_from_lookup(iri=identifier_from_iri(literal_iri))
                    if obj is not None and obj.identifier in acquired:
                        resolved[role] = obj.identifier
                        continue
                if isinstance(spec, Selector):
                    resolve_proc = env.process(spec.resolve(env))
                    result = yield resolve_proc
                    resolve_proc = None
                elif isinstance(spec, str):
                    # lock via LiteralSelector to keep path uniform
                    resolve_proc = env.process(LiteralSelector(spec).resolve(env))
                    result = yield resolve_proc
                    resolve_proc = None
                else:
                    raise TypeError(
                        f"Participant '{role}' has unsupported type {type(spec)}"
                    )

                if result is None:
                    raise SelectorCancelled(f"selector-cancelled:{role}")
                iri, req = result

                if iri in acquired:
                    req.resource.release(req)  # we already hold the lock
                    req = acquired[iri]
                else:
                    acquired[iri] = req
                    self.locks.append(req)
                    if iri not in self.lock_acquired_times:
                        self.lock_acquired_times[iri] = env.now
    
                resolved[role] = iri
    
            # overwrite participant_* fields with pure strings -------------
            self.resolved_resources = resolved  # remember bindings
            _write_participant_iris(self, resolved)
    
            # also expose them via resources[] for backward compatibility
            unique_resources: list[str] = []
            seen_resources: set[str] = set()
            for iri in resolved.values():
                if iri in seen_resources:
                    continue
                seen_resources.add(iri)
                unique_resources.append(iri)
            self.resources = unique_resources
    
            # build list of graph edits now that everything is bound -------
            self.operation_effects = self.get_operation_effects()
    
            created_iris = [e.instance_1_iri for e in self.operation_effects
                            if e.type is UnitaryEditType.CREATE]
            if len(created_iris) != len(set(created_iris)):
                raise RuntimeError(
                    f"Duplicate CREATE IRIs detected in {self.identifier}: {created_iris}"
                )
    
        except simpy.Interrupt as intr:
            # pre_act() can be interrupted if the parent operation is cancelled while
            # blocked on selector resolution / lock acquisition. Clean up and
            # propagate so callers can abort safely.
            # Ensure any in-flight selector resolution process is also cancelled;
            # otherwise it may later acquire a lock/store item and strand it.
            if resolve_proc is not None:
                try:
                    # SimPy Process has `.triggered` once it is done; no `is_alive`.
                    if not getattr(resolve_proc, "triggered", True):
                        resolve_proc.interrupt(getattr(intr, "cause", None) or "cancel")
                except Exception:
                    pass

            self._release_all_locks(env)
            raise
        except Exception:
            self._release_all_locks(env)
            raise
    def post_act(self, env: simpy.Environment):
        """
        Release all held locks and reinsert surviving objects into their FilterStores.

        Robust to ANNIHILATE:
        - If an object was annihilated during `apply(...)`, its Resource will have
          been removed from the runtime context. Reverse lookups for such resources
          will fail; we still release the lock (on the Resource instance we hold),
          but skip reinsertion (object is not present).
        """
        if self.sim_state is not _OpState.RUNNING:
            raise RuntimeError(f"{self.identifier}: post_act called in state {self.sim_state}")

        reinserts: dict[str, Any] = {}
        for req in self.locks:
            obj = None
            # Try to map Resource -> LabObject; this may fail if the object was annihilated.
            try:
                obj = get_object_for_resource(req.resource)
            except KeyError:
                # Resource no longer registered (likely ANNIHILATE). We can still release the lock below.
                pass
            except Exception:
                obj = None

            resource = getattr(req, "resource", None)
            if resource is not None:
                try:
                    users = getattr(resource, "users", None)
                    queue = getattr(resource, "queue", None)
                    if users is not None and req in users:
                        resource.release(req)
                    elif hasattr(req, "cancel"):
                        req.cancel()
                    elif queue is not None and req in queue:
                        try:
                            queue.remove(req)
                        except ValueError:
                            pass
                except Exception:
                    pass

            # Reinsert only if we successfully mapped and the object is still present
            if obj is not None and getattr(obj, "is_present", {False}) == {True}:
                reinserts[obj.identifier] = obj

        # Fallback reinsertion for cases where resource reverse lookup failed.
        for iri in self.resources:
            obj = KnowledgeGraph.get_object_from_lookup(iri)
            if obj is not None and getattr(obj, "is_present", {False}) == {True}:
                reinserts[obj.identifier] = obj

        for obj in reinserts.values():
            FilterStoreRegistry.safe_put_obj_into_filter_store(
                obj,
                env,
                context=f"operation.post_act:{self.identifier}",
            )

        self.locks.clear()
        self.sim_state = _OpState.FINISHED

    def cleanup(self, env: simpy.Environment) -> None:
        """Best-effort cleanup for abort/interrupt paths.

        Normal executions should call :meth:`post_act` (which enforces that the
        operation reached the RUNNING state). However, interruptions can occur
        while the operation is still waiting to start (scheduled start time,
        precedents) or while it is preparing (PREPARED) and holding some locks.

        This helper:
        - If RUNNING: delegates to ``post_act``.
        - If PREPARED/NEW: releases any acquired locks and attempts to reinsert
          surviving objects into their FilterStores.
        - If FINISHED: no-op.
        """

        if self.sim_state is _OpState.FINISHED:
            return
        if self.sim_state is _OpState.RUNNING:
            # Use the strict/normal post_act path.
            self.post_act(env)
            return

        # NEW or PREPARED: release any locks we may have acquired so far.
        self._release_all_locks(env)
        self.sim_state = _OpState.FINISHED

    @field_validator("temporal_cost")
    @classmethod
    def _validate_temporal_cost(cls, value: Optional[float]) -> Optional[float]:
        if value is None:
            return value
        try:
            is_finite = math.isfinite(value)
        except TypeError as exc:
            raise ValueError("temporal_cost must be finite and >= 0") from exc
        if not is_finite or value < 0:
            raise ValueError("temporal_cost must be finite and >= 0")
        return value

    @field_validator("scheduled_start_time")
    @classmethod
    def _validate_scheduled_start_time(cls, value: Optional[float]) -> Optional[float]:
        if value is None:
            return value
        try:
            is_finite = math.isfinite(value)
        except TypeError as exc:
            raise ValueError("scheduled_start_time must be finite and >= 0") from exc
        if not is_finite or value < 0:
            raise ValueError("scheduled_start_time must be finite and >= 0")
        return value

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        validate_assignment=False,
        extra="forbid",
    )


Operation.model_rebuild()
UnitaryEdit.model_rebuild()
Presumption.model_rebuild()
