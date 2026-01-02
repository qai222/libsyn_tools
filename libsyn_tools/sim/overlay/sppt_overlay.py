from __future__ import annotations

"""
SPPT overlay provider: materialize time-boxed event skeletons from lifecycle callbacks.

What this provider emits
------------------------
For each completed Operation (via lifecycle callbacks):
  • a Process node:        <lib:Process/{op_id}> a lib:Process .
  • a TimeInterval node:   <lib:Interval/{op_id}> a lib:TimeInterval ;
                              lib:has_begin_time t0 ; lib:has_end_time t1 .
  • a Process→Interval link:    lib:has_interval
  • Process→Substance participants: lib:has_participant

Deliberately *not* emitted here (leave to domain-specific providers):
  • lib:occurs_in (Place), lib:precedes/causes
  • role reification (instrument/source/destination)
  • derived shortcuts (e.g., in_contact_with, exposed_to)
  • currentVolume (there is a separate provider)

Design notes
------------
• No import of Simulation/OperationProcess to avoid circular deps; we depend only on
  attributes observed on the proc object: proc.env.now, proc.operation.identifier,
  proc.operation.resources.
• Participants are emitted with the canonical lib: namespace, i.e., LIB[iri_string]
  so they match test expectations and any code that builds LIB[...] IRIs.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, List, Any

from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, XSD

from libsyn_tools.sim.knowledge_graph.ontology import (
    Has_participant, Has_interval, Has_begin_time, Has_end_time
)

LIB = Namespace("https://libsyn-sim/kg/")


@dataclass
class _OpSpan:
    t0: Optional[float] = None
    t1: Optional[float] = None
    participants: List[str] = field(default_factory=list)


class SPPTOverlayProvider:
    """
    Collect operation spans via lifecycle callbacks and materialize them as SPPT events.
    Register this provider with the EffectEngine to include its triples in the overlay union.
    """

    def __init__(self, callbacks) -> None:
        self._spans: Dict[str, _OpSpan] = {}
        self._g = Graph()
        callbacks.on_operation_start.append(self._on_start)
        callbacks.on_operation_end.append(self._on_end)

    # ---- lifecycle subscribers ----
    def _on_start(self, proc: Any) -> None:
        op = proc.operation
        span = self._spans.setdefault(op.identifier, _OpSpan())
        span.t0 = float(proc.env.now)
        # After pre_act(), resources hold resolved participant IRIs (strings)
        if op.resources:
            span.participants = list(op.resources)

    def _on_end(self, proc: Any) -> None:
        op = proc.operation
        span = self._spans.setdefault(op.identifier, _OpSpan())
        span.t1 = float(proc.env.now)
        if not span.participants and op.resources:
            span.participants = list(op.resources)
        self._materialize_span(op.identifier, span)

    def _materialize_span(self, op_id: str, span: _OpSpan) -> None:
        proc_iri = LIB[f"Process/{op_id}"]
        int_iri = LIB[f"Interval/{op_id}"]
        self._g.remove((proc_iri, None, None))
        self._g.remove((int_iri, None, None))

        if span.t0 is None or span.t1 is None:
            return

        self._g.add((proc_iri, RDF.type, LIB.Process))
        self._g.add((int_iri, RDF.type, LIB.TimeInterval))
        self._g.add(
            (
                int_iri,
                URIRef(Has_begin_time.predicate_iri),
                Literal(span.t0, datatype=XSD.double),
            )
        )
        self._g.add(
            (
                int_iri,
                URIRef(Has_end_time.predicate_iri),
                Literal(span.t1, datatype=XSD.double),
            )
        )
        self._g.add((proc_iri, URIRef(Has_interval.predicate_iri), int_iri))
        for p in span.participants:
            self._g.add((proc_iri, URIRef(Has_participant.predicate_iri), LIB[p]))

    # ---- overlay snapshot ----
    def snapshot(self) -> Graph:
        """
        Return the current rdflib.Graph snapshot of SPPT events.
        This is a persistent graph; callers must treat it as read-only.
        """
        return self._g
