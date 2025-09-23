from __future__ import annotations

"""
SPPT overlay provider: builds time-boxed Process nodes from lifecycle callbacks.

For each operation:
- On start: record t0
- On end:   record t1; emit a Process instance with:
    lib:has_participant (for each resolved participant IRI)
    lib:has_interval    → TimeInterval node with begin/end

Note: No import of Simulation/OperationProcess (avoid circular). We only
rely on the attributes actually used on the 'proc' argument (env.now,
operation.identifier, operation.resources).
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, List, Any

from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, XSD

from libsyn_tools.sim.knowledge_graph.sppt import (
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
    """

    def __init__(self, callbacks) -> None:
        self._spans: Dict[str, _OpSpan] = {}
        callbacks.on_operation_start.append(self._on_start)
        callbacks.on_operation_end.append(self._on_end)

    # ---- lifecycle subscribers ----
    def _on_start(self, proc: Any) -> None:
        op = proc.operation
        span = self._spans.setdefault(op.identifier, _OpSpan())
        span.t0 = float(proc.env.now)
        if op.resources:
            span.participants = list(op.resources)

    def _on_end(self, proc: Any) -> None:
        op = proc.operation
        span = self._spans.setdefault(op.identifier, _OpSpan())
        span.t1 = float(proc.env.now)
        if not span.participants and op.resources:
            span.participants = list(op.resources)

    # ---- overlay snapshot ----
    def snapshot(self) -> Graph:
        g = Graph()
        for op_id, span in self._spans.items():
            if span.t0 is None or span.t1 is None:
                continue

            proc_iri = LIB[f"Process/{op_id}"]
            int_iri = LIB[f"Interval/{op_id}"]

            # rdf:type assertions for Process and TimeInterval
            g.add((proc_iri, RDF.type, LIB.Process))
            g.add((int_iri, RDF.type, LIB.TimeInterval))

            # interval begin/end (use proper predicate IRIs)
            g.add((int_iri, URIRef(Has_begin_time.predicate_iri),
                   Literal(span.t0, datatype=XSD.double)))
            g.add((int_iri, URIRef(Has_end_time.predicate_iri),
                   Literal(span.t1, datatype=XSD.double)))

            # link process -> interval
            g.add((proc_iri, URIRef(Has_interval.predicate_iri), int_iri))

            # participants (use proper predicate IRI)
            for p in span.participants:
                # Use namespaced IRIs so triples match tests (LIB.has_participant, LIB[p])
                g.add((proc_iri, URIRef(Has_participant.predicate_iri), LIB[p]))

        return g
