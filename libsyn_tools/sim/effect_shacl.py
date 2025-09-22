from __future__ import annotations

from typing import Optional, Any, List
from uuid import uuid4

from pydantic import BaseModel, Field
from rdflib import Graph, Literal, URIRef, BNode
from rdflib.namespace import SH, RDF

"""
Uniform SHACL/engine violation record and small helpers.

A single, typed record schema for *both* SHACL (policy/soft) and engine-mechanical (hard) violations so post-hoc
triage and replay are uniform.

Engine emits `origin="ENGINE"` with `severity="hard"`; 

SHACL validation emits `origin="SHACL"` with `severity="soft"`.

`violation_id` for dedupe/trace, 
`batch_id` ties a group of edits (one operation application) together,
`edit_fingerprints` captures the edit batch content without coupling to internal objects.
"""


class SHACLViolationRecord(BaseModel):
    # universal fields (engine + SHACL)
    violation_id: str = Field(default_factory=lambda: str(uuid4()))
    sim_time: float
    operation_id: str
    origin: str = Field(..., description='One of {"ENGINE","SHACL"}')
    severity: str = Field(..., description='One of {"hard","soft"}')
    disposition: str = Field(
        ...,
        description='One of {"aborted","committed","ignored"} '
                    '(engine hard errors abort; SHACL commits)'
    )

    # optional run metadata
    batch_id: Optional[str] = None
    edit_fingerprints: Optional[List[str]] = None
    seed: Optional[int] = None

    # SHACL-specific decoration (optional for ENGINE)
    shape_iri: Optional[str] = Field(default=None)
    focus_iri: Optional[str] = Field(default=None)
    message: Optional[str] = Field(default=None)
    report_graph_ttl: Optional[str] = Field(default=None)


def _node_to_str(node: Optional[Any]) -> Optional[str]:
    """
    Convert an rdflib node (URIRef, Literal, or BNode) to a human-readable string.
    Blank nodes get a skolemised “_:<id>” form; None stays None.
    """
    if node is None:
        return None
    if isinstance(node, (URIRef, Literal)):
        return str(node)
    if isinstance(node, BNode):
        return f"_:{node}"
    return str(node)


def _first(graph: Graph, subject: Optional[Any], predicate: URIRef) -> Optional[str]:
    """
    Return the *first* object matching (?subject, predicate, ?o).
    If subject is None, the graph is scanned for ANY subject ?s.
    """
    if subject is not None:
        obj = graph.value(subject, predicate, any=False)
        return _node_to_str(obj)
    # no subject given → find the first triple with that predicate
    for _, _, o in graph.triples((None, predicate, None)):
        return _node_to_str(o)
    return None


def _iter_validation_results(r_graph: Graph):
    """Yield every node in r_graph that is rdf:type sh:ValidationResult."""
    for vr in r_graph.subjects(RDF.type, SH.ValidationResult):
        yield vr
