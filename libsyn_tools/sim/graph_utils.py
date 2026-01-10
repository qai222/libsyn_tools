from __future__ import annotations

from collections.abc import Sequence

from rdflib import Graph
from rdflib.graph import ReadOnlyGraphAggregate


class HashableReadOnlyGraphAggregate(ReadOnlyGraphAggregate):
    """
    Read-only graph aggregate that can be hashed (for PySHACL internals).
    """

    def __hash__(self) -> int:  # pragma: no cover - trivial
        return id(self)


class UnionGraphView(Graph):
    """
    Read-only Graph facade over a ReadOnlyGraphAggregate.

    This avoids treating the aggregate as a ConjunctiveGraph so consumers that
    clone Graphs (e.g., PySHACL) iterate over triples instead of contexts.
    """

    def __init__(self, aggregate: ReadOnlyGraphAggregate):
        super().__init__()
        self._aggregate = aggregate
        if aggregate.namespace_manager is not None:
            self.namespace_manager = aggregate.namespace_manager

    def __iter__(self):
        return iter(self._aggregate)

    def __len__(self) -> int:  # pragma: no cover - trivial delegation
        return len(self._aggregate)

    def __contains__(self, triple) -> bool:  # pragma: no cover - trivial delegation
        return triple in self._aggregate

    def triples(self, triple):
        return self._aggregate.triples(triple)

    def value(self, subject=None, predicate=None, object=None, default=None, any=True):
        return self._aggregate.value(
            subject=subject,
            predicate=predicate,
            object=object,
            default=default,
            any=any,
        )

    def query(self, *args, **kwargs):
        return self._aggregate.query(*args, **kwargs)


def union_view(*graphs: Graph | None) -> Graph:
    """
    Return a read-only union view over the provided graphs.

    This is a lightweight view (no triple copying) and should not be mutated.
    """
    return union_view_many(graphs)


def union_view_many(graphs: Sequence[Graph | None]) -> Graph:
    """
    Return a read-only union view over the provided graphs.

    This is a lightweight view (no triple copying) and should not be mutated.
    """
    usable_graphs = [graph for graph in graphs if graph is not None]
    return HashableReadOnlyGraphAggregate(usable_graphs)


def union_view_for_shacl(graphs: Sequence[Graph | None]) -> Graph:
    """
    Return a read-only Graph view suitable for SHACL validation.

    This wraps the aggregate in a Graph facade so PySHACL clones by iterating
    over triples instead of using empty context lists.
    """
    return UnionGraphView(union_view_many(graphs))
