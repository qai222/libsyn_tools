from __future__ import annotations

from pathlib import Path

from rdflib import Graph


def load_core_shapes() -> Graph:
    shapes_path = Path(__file__).with_name("core_shapes.ttl")
    return Graph().parse(shapes_path, format="turtle")
