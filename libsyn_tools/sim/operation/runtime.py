from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Dict, TYPE_CHECKING

import simpy
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph import LabObject, BaseClass
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit

if TYPE_CHECKING:
    pass

class RuntimeContext:
    """
    Per-simulation runtime storage for artifacts that must not be global.

    Attached to each simpy.Environment as ``env._libsyn_runtime_ctx``.
    """

    def __init__(self) -> None:
        self.resource_map: Dict[str, simpy.Resource] = {}
        self.runtime_cache: Dict[str, _RuntimeState] = {}
        self.filter_stores: Dict[str, simpy.FilterStore] = {}


def get_runtime_context(env: simpy.Environment, *, create: bool = True) -> RuntimeContext:
    ctx = getattr(env, "_libsyn_runtime_ctx", None)
    if ctx is None:
        if not create:
            raise RuntimeError("RuntimeContext not attached to environment")
        ctx = RuntimeContext()
        setattr(env, "_libsyn_runtime_ctx", ctx)
    return ctx


# Deprecated module globals (kept for compatibility only; do not use internally).
_RESOURCE_MAP: Dict[str, simpy.Resource] = {}


def register_object_as_resource(obj: LabObject, env: simpy.Environment):
    ctx = get_runtime_context(env)
    if obj.instance_iri in ctx.resource_map:
        raise RuntimeError(f"Resource for {obj.instance_iri} already registered")
    ctx.resource_map[obj.instance_iri] = simpy.Resource(env, capacity=1)


def get_resource_for_object(obj: LabObject, env: simpy.Environment | None = None) -> simpy.Resource:
    if env is None:
        runtime = getattr(obj, "_runtime", None)
        env = getattr(runtime, "env", None) if runtime else None
    if env is None:
        raise RuntimeError(f"No environment available for {obj.instance_iri}")
    try:
        return get_runtime_context(env, create=False).resource_map[obj.instance_iri]
    except KeyError:
        raise RuntimeError(
            f"No simpy.Resource registered for {obj.instance_iri}. "
            "Call ActionSimulation.build_resources() before env.run()."
        )


@dataclass(frozen=True)
class OperationUsageRecord:
    operation_id: str
    operation_type: str
    pool_type: str | None
    lock_acquired_sim_time: float | None
    pool_type_error: str | None = None


class _RuntimeState:
    def __init__(self, env: simpy.Environment, obj: LabObject):
        self.env = env
        self.obj = obj
        self.lock = get_resource_for_object(obj, env)
        self.recent_edits: deque[UnitaryEdit] = deque(maxlen=128)
        self.recent_operations: deque["Operation"] = deque(maxlen=128)
        self.recent_operation_records: deque[OperationUsageRecord] = deque(maxlen=128)
        # TODO we could use weakref but is it necessary? or maybe just use (timestamp, action id)?
        # from weakref import ref
        # self.recent_actions: deque[ref[Action]] = deque(maxlen=128)


# Deprecated module global (kept for compatibility only; do not use internally).
_RUNTIME_CACHE: Dict[str, _RuntimeState] = {}


def get_runtime_state(obj: LabObject, env: simpy.Environment) -> _RuntimeState:
    """
    Ensure a LabObject has an attached runtime state for *this* env and
    return it.  Idempotent & fast (dict lookup).
    """
    if not _needs_runtime_tracking(obj):
        raise RuntimeError(
            f"{obj.__class__.__name__} objects do not participate in "
            "locking / runtime state."
        )

    ctx = get_runtime_context(env)
    rs = ctx.runtime_cache.get(obj.identifier)
    if rs is None:
        rs = _RuntimeState(env, obj)
        ctx.runtime_cache[obj.identifier] = rs
        # allow convenient access: obj._runtime  (purely in-memory)
        setattr(obj, "_runtime", rs)
    return rs


def _needs_runtime_tracking(obj: BaseClass) -> bool:
    """
    Return True for objects that should be locked, pooled and keep
    a _RuntimeState entry (i.e. real LabObjects – glassware, pumps,
    robots …).  PortionOfMaterial and other data-only nodes return False.
    """

    return isinstance(obj, LabObject)


def get_object_for_resource(res: simpy.Resource) -> LabObject:
    """
    Reverse lookup: simpy.Resource → LabObject.

    Used when releasing locks so we can put the object back into the
    correct FilterStore (FIX 3.3).
    """
    env = getattr(res, "env", None) or getattr(res, "_env", None)
    if env is None:
        raise RuntimeError("Resource has no associated environment")
    ctx = get_runtime_context(env, create=False)
    for iri, r in ctx.resource_map.items():
        if r is res:
            return KnowledgeGraph.get_object_from_lookup(iri)
    raise KeyError("Resource not registered in runtime context")
