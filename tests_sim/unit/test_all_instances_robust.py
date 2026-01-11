from __future__ import annotations

from libsyn_tools.sim.knowledge_graph.ontology import LabObject, Substance


def test_all_instances_handles_missing_object_lookup(monkeypatch):
    monkeypatch.setattr(Substance, "object_lookup", None, raising=False)
    monkeypatch.setattr(LabObject, "object_lookup", None, raising=False)

    substance_instances = list(Substance.all_instances())
    lab_instances = list(LabObject.all_instances())

    assert all(hasattr(obj, "instance_iri") for obj in substance_instances)
    assert all(hasattr(obj, "instance_iri") for obj in lab_instances)
