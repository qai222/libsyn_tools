# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_effect_shacl.py ###
from __future__ import annotations

from rdflib import Graph, Namespace, BNode, Literal
from rdflib.namespace import RDF, SH

from libsyn_tools.sim.effect_shacl import (
    SHACLViolationRecord,
    _first,
    _iter_validation_results,
)


def test_violation_record_schema_defaults():
    rec = SHACLViolationRecord(
        sim_time=1.23,
        operation_id="op-1",
        origin="SHACL",
        severity="soft",
        disposition="committed",
    )
    d = rec.model_dump()
    assert d["violation_id"]
    assert d["operation_id"] == "op-1"
    assert d["origin"] == "SHACL"
    assert d["severity"] == "soft"
    assert d["disposition"] == "committed"
    # SHACL decorations are optional
    assert "shape_iri" in d and d["shape_iri"] is None


def test_iter_helpers_over_synthetic_report():
    g = Graph()
    lib = Namespace("https://libsyn-sim/kg/")

    vr = BNode()
    g.add((vr, RDF.type, SH.ValidationResult))
    g.add((vr, SH.sourceShape, lib.SomeShape))
    g.add((vr, SH.focusNode, lib.SomeFocus))
    g.add((vr, SH.resultMessage, Literal("Oops")))

    vrs = list(_iter_validation_results(g))
    assert len(vrs) == 1

    s = _first(g, vrs[0], SH.sourceShape)
    f = _first(g, vrs[0], SH.focusNode)
    m = _first(g, vrs[0], SH.resultMessage)
    assert s.endswith("SomeShape")
    assert f.endswith("SomeFocus")
    assert m == "Oops"
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_effect_shacl.py ###
