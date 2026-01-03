# ### THIS IS THE START OF CONTENT OF tests_sim/conftest.py ###
from __future__ import annotations

import os
import random
import sys
from typing import Iterator

import pytest
import simpy
from loguru import logger
from rdflib import Graph
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import OperationProcess as _OP
from libsyn_tools.sim.knowledge_graph import (
    LabObject,
    MaterialContainer,
    PortionOfMaterial,
)
from libsyn_tools.sim.operation.runtime import _RESOURCE_MAP, _RUNTIME_CACHE

# Capture the original add_event_log once (before any spawner patches)
_ORIG_ADD_EVENT_LOG = _OP.add_event_log


@pytest.fixture(autouse=True)
def _configure_logging() -> Iterator[None]:
    """
    Configure loguru for test runs (stderr sink).
    Keep quiet by default; override with TEST_LOG_LEVEL=INFO/DEBUG as needed.
    """
    logger.remove()
    logger.add(sys.stderr, level=os.getenv("TEST_LOG_LEVEL", "WARNING"))
    yield
    logger.remove()


@pytest.fixture
def env() -> Iterator[simpy.Environment]:
    """Fresh SimPy Environment per test (for component tests that need it)."""
    yield simpy.Environment()


@pytest.fixture
def seed() -> int:
    """Deterministic seed for tests that need it."""
    return 12345


@pytest.fixture(autouse=True)
def fresh_kg() -> Iterator[None]:
    """
    Fresh KG & runtime per test:

    - clear RDF graph
    - reset per-class object registries to **empty dicts**
    - clear runtime caches and filter stores
    """
    # 1) Clear KG triples
    g: Graph = KnowledgeGraph.graph()
    g.remove((None, None, None))

    # 2) Ensure object_lookup is a dict for each concrete class
    for cls in (PortionOfMaterial, MaterialContainer, LabObject):
        if not hasattr(cls, "object_lookup") or cls.object_lookup is None:
            setattr(cls, "object_lookup", {})
        else:
            try:
                cls.object_lookup.clear()
            except Exception:
                setattr(cls, "object_lookup", {})

    # 3) Clear runtime artefacts
    _RESOURCE_MAP.clear()
    _RUNTIME_CACHE.clear()

    # 4) Deterministic ambient RNG if anyone relies on random()
    random.seed(0)

    yield


@pytest.fixture(autouse=True)
def restore_opprocess_add_event_log() -> Iterator[None]:
    """
    Some spawners patch OperationProcess.add_event_log at the **class** level.
    Restore the original method & clear guard flags after each test.
    """
    yield
    _OP.add_event_log = _ORIG_ADD_EVENT_LOG
    for flag in ("_kginsp_patched", "_interrupt_wrapper_installed"):
        if hasattr(_OP, flag):
            try:
                delattr(_OP, flag)
            except Exception:
                pass

def _clear_all_object_lookups() -> None:
    """Clear global TWA object lookups to keep tests isolated.

    The TWA data model stores instances in per-class `object_lookup` dicts.
    Many sim utilities iterate these lookups, so cross-test pollution can
    create flaky tests.
    """
    cls_lookup = getattr(KnowledgeGraph, "class_lookup", None) or {}
    for cls in cls_lookup.values():
        obj_lookup = getattr(cls, "object_lookup", None)
        if isinstance(obj_lookup, dict):
            obj_lookup.clear()


@pytest.fixture(autouse=True)
def clean_kg_between_tests():
    _clear_all_object_lookups()
    yield
    _clear_all_object_lookups()
# ### THIS IS THE END OF CONTENT OF tests_sim/conftest.py ###
