# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_sppt_overlay_incremental.py ###
from __future__ import annotations

from dataclasses import dataclass

from rdflib import Namespace, Literal, URIRef
from rdflib.namespace import RDF, XSD

from libsyn_tools.sim.lifecycle import LifecycleCallbacks
from libsyn_tools.sim.knowledge_graph.ontology import (
    Has_begin_time,
    Has_end_time,
    Has_interval,
    Has_participant,
)
from libsyn_tools.sim.overlay.sppt_overlay import SPPTOverlayProvider

LIB = Namespace("https://libsyn-sim/kg/")


@dataclass
class _FakeEnv:
    now: float


@dataclass
class _FakeOp:
    identifier: str
    resources: list[str]


@dataclass
class _FakeProc:
    env: _FakeEnv
    operation: _FakeOp


def test_sppt_overlay_snapshot_is_incremental() -> None:
    callbacks = LifecycleCallbacks()
    provider = SPPTOverlayProvider(callbacks)

    op = _FakeOp(identifier="OP1", resources=["resA", "resB"])
    proc = _FakeProc(env=_FakeEnv(now=1.5), operation=op)

    callbacks.emit_operation_start(proc)
    proc.env.now = 3.5
    callbacks.emit_operation_end(proc)

    g1 = provider.snapshot()
    g2 = provider.snapshot()

    assert g1 is g2

    proc_iri = LIB[f"Process/{op.identifier}"]
    int_iri = LIB[f"Interval/{op.identifier}"]

    assert (proc_iri, RDF.type, LIB.Process) in g1
    assert (int_iri, RDF.type, LIB.TimeInterval) in g1
    assert (int_iri, URIRef(Has_begin_time.predicate_iri), Literal(1.5, datatype=XSD.double)) in g1
    assert (int_iri, URIRef(Has_end_time.predicate_iri), Literal(3.5, datatype=XSD.double)) in g1
    assert (proc_iri, URIRef(Has_interval.predicate_iri), int_iri) in g1
    assert (proc_iri, URIRef(Has_participant.predicate_iri), LIB["resA"]) in g1
    assert (proc_iri, URIRef(Has_participant.predicate_iri), LIB["resB"]) in g1

    proc.env.now = 4.25
    callbacks.emit_operation_end(proc)

    assert (int_iri, URIRef(Has_end_time.predicate_iri), Literal(3.5, datatype=XSD.double)) not in g1
    assert (int_iri, URIRef(Has_end_time.predicate_iri), Literal(4.25, datatype=XSD.double)) in g1


def test_sppt_overlay_accepts_canonical_resource_iris() -> None:
    """Resources may already be canonical IRIs (e.g. objects created with
    identifier="https://libsyn-sim/kg/..."), and the overlay should not double-prefix.
    """

    callbacks = LifecycleCallbacks()
    provider = SPPTOverlayProvider(callbacks)

    op = _FakeOp(
        identifier="OP2",
        resources=["https://libsyn-sim/kg/resA", "https://libsyn-sim/kg/resB"],
    )
    proc = _FakeProc(env=_FakeEnv(now=0.0), operation=op)

    callbacks.emit_operation_start(proc)
    proc.env.now = 1.0
    callbacks.emit_operation_end(proc)

    g = provider.snapshot()
    proc_iri = LIB[f"Process/{op.identifier}"]
    assert (proc_iri, URIRef(Has_participant.predicate_iri), LIB["resA"]) in g
    assert (proc_iri, URIRef(Has_participant.predicate_iri), LIB["resB"]) in g
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_sppt_overlay_incremental.py ###
