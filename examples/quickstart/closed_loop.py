from __future__ import annotations

from pathlib import Path
from uuid import uuid4

from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import Simulation
from libsyn_tools.sim.knowledge_graph import (
    MaterialContainer,
    PortionOfMaterial,
    Is_directly_contained_by,
)
from libsyn_tools.sim.operation.unitary_edit import Create, AddObjectProperty
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByVolume
from libsyn_tools.sim.operation_preset.wait import Wait
from libsyn_tools.sim.remediation import make_drain_to_capacity
from libsyn_tools.sim.shapes import load_core_shapes
from libsyn_tools.sim.spawner import PolicyEnforcerSpawner, ValidationAuditSpawner


def run(out_dir: str | Path | None = None) -> Simulation:
    base = "https://libsyn-sim/kg/"
    shapes = load_core_shapes()

    source = MaterialContainer(identifier=f"{base}source")
    destination = MaterialContainer(identifier=f"{base}destination")
    destination.has_capacity.add(1.0)
    waste = MaterialContainer(identifier=f"{base}waste")
    device = MaterialContainer(identifier=f"{base}device")
    pom = PortionOfMaterial(identifier=f"{base}pom-source")
    pom.add_chemical(Chemical(mass=1.0, density=1.0))

    for obj in (source, destination, waste, device, pom):
        Create(instance_1_iri=obj.identifier).apply()

    AddObjectProperty(
        instance_1_iri=pom.identifier,
        instance_2_iri=source.identifier,
        property_iri=Is_directly_contained_by.predicate_iri,
    ).apply()

    transfer = TransferMaterialByVolume(
        identifier="fill-dest",
        participant_source=source.identifier,
        participant_destination=destination.identifier,
        participant_device=device.identifier,
        transfer_volume=1.0,
        temporal_cost=0.0,
    )

    audit_anchor = Wait(identifier="VALIDATION_AUDIT", temporal_cost=0.0)
    sim = Simulation([transfer, audit_anchor], shacl_shapes=shapes)

    def _drain_factory(record):
        if record.focus_iri is None:
            return None
        op = make_drain_to_capacity(record.focus_iri, waste.identifier)
        op.identifier = f"drain-{uuid4().hex[:6]}"
        op.temporal_cost = 0.0
        return op

    ValidationAuditSpawner(inspect_interval=1.0).attach(sim)
    PolicyEnforcerSpawner(shape_dispatch={f"{base}CapacityShape": _drain_factory}).attach(sim)

    if out_dir:
        sim.run_and_report(out_dir, until=2.1)
    else:
        sim.run(until=2.1)
    return sim


if __name__ == "__main__":
    run()
