from __future__ import annotations

from collections import deque
from typing import Dict, TYPE_CHECKING

import simpy
from twa.data_model.base_ontology import KnowledgeGraph
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit
from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject, BaseClass

if TYPE_CHECKING:
    from libsyn_tools.sim.operation.operation import Operation

# each lab object is mapped to a capacity=1 resource
_RESOURCE_MAP: Dict[str, simpy.Resource] = {}


def register_object_as_resource(obj: LabObject, env: simpy.Environment):
    if obj.instance_iri in _RESOURCE_MAP:
        raise RuntimeError(f"Resource for {obj.instance_iri} already registered")
    _RESOURCE_MAP[obj.instance_iri] = simpy.Resource(env, capacity=1)


def get_resource_for_object(obj: LabObject) -> simpy.Resource:
    try:
        return _RESOURCE_MAP[obj.instance_iri]
    except KeyError:
        raise RuntimeError(
            f"No simpy.Resource registered for {obj.instance_iri}. "
            "Call ActionSimulation.build_resources() before env.run()."
        )


class _RuntimeState:
    def __init__(self, env: simpy.Environment, obj: LabObject):
        self.env = env
        self.obj = obj
        self.lock = get_resource_for_object(obj)
        self.recent_edits: deque[UnitaryEdit] = deque(maxlen=128)
        self.recent_operations: deque["Operation"] = deque(maxlen=128)
        # TODO we could use weakref but is it necessary? or maybe just use (timestamp, action id)?
        # from weakref import ref
        # self.recent_actions: deque[ref[Action]] = deque(maxlen=128)


# global cache: iri ➜ runtime state
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

    rs = _RUNTIME_CACHE.get(obj.identifier)
    if rs is None:
        rs = _RuntimeState(env, obj)
        _RUNTIME_CACHE[obj.identifier] = rs
        # allow convenient access: obj._runtime  (purely in-memory)
        setattr(obj, "_runtime", rs)
    return rs

def _needs_runtime_tracking(obj: BaseClass) -> bool:
    """
    Return True for objects that should be locked, pooled and keep
    a _RuntimeState entry (i.e. real LabObjects – glassware, pumps,
    robots …).  PortionOfMaterial and other data-only nodes return False.
    """
    from libsyn_tools.sim.knowledge_graph.physical_entities import LabObject

    return isinstance(obj, LabObject)

def get_object_for_resource(res: simpy.Resource) -> LabObject:
    """
    Reverse lookup: simpy.Resource → LabObject.

    Used when releasing locks so we can put the object back into the
    correct FilterStore (FIX 3.3).
    """
    for iri, r in _RESOURCE_MAP.items():
        if r is res:
            return KnowledgeGraph.get_object_from_lookup(iri)
    raise KeyError("Resource not registered in _RESOURCE_MAP")