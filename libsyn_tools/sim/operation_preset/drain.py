from __future__ import annotations

from typing import List

from loguru import logger
from pydantic import Field
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import MaterialContainer, PortionOfMaterial
from libsyn_tools.sim.knowledge_graph import Is_directly_contained_by, LabObject
from libsyn_tools.sim.operation.effects_dsl import EffectsBuilder
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit

_EPS = 1e-6


class DrainExcess(Operation):
    """Trim a container down to `target_volume` by moving the surplus into *destination*.

    Parameters
    ----------
    participant_source:
        Container that is too full.
    participant_destination:
        Container that receives the surplus.
    target_volume:
        Volume **after** the operation (same units as Chemical.volume).

    Notes
    -----
    The operation drains one PortionOfMaterial (POM) at a time until the target
    volume is reached.

    Important correctness detail:
        If an entire POM is drained, we do *not* create a zero-volume residual POM
        (``PortionOfMaterial.get_portion_by_volume(0)`` would raise).
    """

    participant_source: StrOrSelector
    participant_destination: StrOrSelector
    participant_device: StrOrSelector | None = None  # optional pipette, pump, …

    target_volume: float = Field(..., gt=0)

    def get_operation_effects(self) -> List[UnitaryEdit]:
        src = KnowledgeGraph.get_object_from_lookup(self.participant_source)
        dst = KnowledgeGraph.get_object_from_lookup(self.participant_destination)

        if not isinstance(src, MaterialContainer):
            raise RuntimeError("DrainExcess: source must be a MaterialContainer")
        if not isinstance(dst, MaterialContainer):
            raise RuntimeError("DrainExcess: destination must be a MaterialContainer")

        builder = EffectsBuilder()
        prop_iri = Is_directly_contained_by.predicate_iri

        current_v = src.directly_contained_pom_volume
        surplus = current_v - self.target_volume
        if surplus <= _EPS:
            logger.info(f"{src.identifier} already ≤ target volume")
            return []

        poms = LabObject.get_directly_contained_individuals(src, PortionOfMaterial, only_present=True)
        for pom in poms:
            if surplus <= _EPS:
                break

            take = min(pom.volume, surplus)
            if take <= _EPS:
                continue

            portion = pom.get_portion_by_volume(take)
            residual_vol = pom.volume - take
            residual = pom.get_portion_by_volume(residual_vol) if residual_vol > _EPS else None
            surplus -= take

            # remove original POM from the source
            builder.unlink(pom.identifier, prop_iri, src.identifier)
            builder.annihilate(pom.identifier)

            # residual stays in source (skip if drained fully)
            if residual is not None:
                builder.create(residual.identifier)
                builder.link(residual.identifier, prop_iri, src.identifier)

            # moved portion goes to destination
            builder.create(portion.identifier)
            builder.link(portion.identifier, prop_iri, dst.identifier)

        logger.info(
            f"DrainExcess: moved {current_v - self.target_volume:.3g} "
            f"from {src.identifier} to {dst.identifier}"
        )
        return builder.build()
