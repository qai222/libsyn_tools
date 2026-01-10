from __future__ import annotations

from rdflib import Graph, Namespace, URIRef

from libsyn_tools.sim.graph_utils import union_view_for_shacl


def test_union_graph_view_behaves_like_graph() -> None:
    lib = Namespace("https://libsyn-sim/kg/")
    g1 = Graph()
    g2 = Graph()
    g1.add((lib.a, lib.p, lib.b))
    g2.add((lib.a, lib.p, lib.c))

    union = union_view_for_shacl([g1, g2])

    objs = {obj for _, _, obj in union.triples((lib.a, lib.p, None))}
    assert objs == {lib.b, lib.c}
    assert union.value(subject=lib.a, predicate=lib.p) in {lib.b, lib.c}

    query_res = union.query(
        "SELECT ?o WHERE { ?s ?p ?o }",
        initBindings={"s": URIRef(lib.a), "p": URIRef(lib.p)},
    )
    found = {row.o for row in query_res}
    assert found == {lib.b, lib.c}
