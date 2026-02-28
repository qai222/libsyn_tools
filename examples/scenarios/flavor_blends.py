from __future__ import annotations

import sys
from pathlib import Path
from pprint import pprint, pformat
from typing import Tuple, List, Dict, Any
import pandas as pd
from scipy.stats import bootstrap
import numpy as np
from loguru import logger
from pydantic import BaseModel, Field
from twa.data_model.base_ontology import KnowledgeGraph
import matplotlib.pyplot as plt
import seaborn as sns
from libsyn_tools.chem_schema import Chemical
from libsyn_tools.sim import (
    Create,
    LabObject,
    MaterialContainer,
    PortionOfMaterial,
    Simulation, AddObjectProperty, Is_directly_contained_by, SimFunctionalDataProperty, ChangeDataProperty, Operation
)
from libsyn_tools.sim.operation.operation import Operation
from libsyn_tools.sim.operation_preset import TransferMaterialByVolume

sns.set_theme(style="whitegrid")
logger.remove()
# logger.add(sys.stderr, level="INFO")
logger.add(sys.stderr, level="WARNING")


# TODO right now we do not consider rack swap, the only operations are liquid transfer and refill off-bed reservoirs

class FlavorBlendSim(BaseModel):
    t_single_liquid_transfer_off_bed: float = Field(default=10.0)
    """ single liquid transfer from off bed reservoir temporal cost in s """

    t_wash: float = Field(default=45.0)
    """ time used to wash the flow path (needle + syringe) """

    # t_rack_swap: float = Field(default=120.0)
    # """ time used to swap one rack """

    t_reservoir_refill: float = Field(default=300.0)
    """ time used to top off a reservoir container """

    v_reservoir_capacity: float = Field(default=100.0)
    """ reservoir capacity in mL """

    v_sample_vial: float = Field(default=7.0)
    """ sample (destination) vial capacity in mL """

    v_sample_target: float = Field(default=3.0)
    """ sample final volume in mL """

    p_reservoir_refill_threshold: float = Field(default=10.0)
    """ percentage below which reservoir is refilled """

    # n_flavor_vials: int = Field(default=210)
    # """ number of vials for flavor blends """
    #
    # n_terpene_vials: int = Field(default=210)
    # """ number of vials for terpene blends """
    #
    # n_rack_slots: int = Field(default=6)
    # """ number of rack slots on the platform """

    n_vials_per_rack: int = Field(default=70)
    """ number of vials per rack, default from Gilson rack code 22 """

    n_terpene_rack: int = Field(default=3)
    """ number of racks for terpene samples """

    n_flavor_rack: int = Field(default=3)
    """ number of racks for flavor samples """

    reservoir_map: dict[str, str] = Field(default_factory=dict)
    """ ingredient name to reservoir iri """

    formulations: dict[str, dict[str, float]] = Field(default_factory=dict)
    """ d[formulation_name][ingredient_name] -> ingredient volume """

    vial_formulation_map: dict[str, str] = Field(default_factory=dict)
    """ d[vial_iri] -> formulation name """

    seed: int = Field(default=42)
    """ random seed """

    def random_formulation(self, ingredients: list[str], rng: np.random.RandomState) -> dict[str, float]:
        k = len(ingredients)
        vols = rng.dirichlet(np.ones(k), size=1)[0] * self.v_sample_target
        return dict(zip(ingredients, vols.tolist()))

    @classmethod
    def random_init(
            cls, n_terpene_formulations: int = 7, n_flavor_formulations: int = 7,
            n_terpene_ingredients: int = 30, n_flavor_ingredients: int = 20,
    ):
        sim = cls()
        rng = np.random.RandomState(seed=sim.seed)
        terpene_ingredients = [f"ING_TERPENE_{i}" for i in range(n_terpene_ingredients)]
        flavor_ingredients = [f"ING_FLAVOR_{i}" for i in range(n_flavor_ingredients)]
        sim.reservoir_map = {ing_name: f"reservoir-{ing_name}" for ing_name in terpene_ingredients + flavor_ingredients}

        terpene_formulations = []
        flavor_formulations = []

        for i in range(n_terpene_formulations):
            formulation_name = f"TERPENE_FORMULATION_{i}"
            sim.formulations[formulation_name] = sim.random_formulation(ingredients=terpene_ingredients, rng=rng)
            terpene_formulations.append(formulation_name)
        for i in range(n_flavor_formulations):
            formulation_name = f"FLAVOR_FORMULATION_{i}"
            sim.formulations[formulation_name] = sim.random_formulation(ingredients=flavor_ingredients, rng=rng)
            flavor_formulations.append(formulation_name)

        for i in range(sim.n_flavor_rack * sim.n_vials_per_rack):
            f = str(rng.choice(flavor_formulations, 1)[0])
            vial_iri = f"vial_{str(i).zfill(4)}"
            sim.vial_formulation_map[vial_iri] = f

        for i in range(sim.n_flavor_rack * sim.n_vials_per_rack,
                       sim.n_flavor_rack * sim.n_vials_per_rack + sim.n_terpene_rack * sim.n_vials_per_rack):
            f = str(rng.choice(terpene_formulations, 1)[0])
            vial_iri = f"vial_{str(i).zfill(4)}"
            sim.vial_formulation_map[vial_iri] = f
        return sim


class Is_clean(SimFunctionalDataProperty): pass


class LiquidTransferDevice(LabObject):
    is_clean: Is_clean[bool] = Field(default={True, })


class Wash(Operation):
    operation_effects_description: str = "wash the transfer device between transferring two different ingredients"
    participant_wash_station: str
    participant_transfer_device: str

    def get_operation_effects(self):
        return [
            ChangeDataProperty(instance_1_iri=self.participant_transfer_device, property_iri=Is_clean.predicate_iri,
                               data_value=True)
        ]


class ReservoirRefill(Operation):
    operation_effects_description: str = "top off an ingredient reservoir"
    participant_platform: str
    participant_transfer_device: str  # add to block
    participant_reservoir: str

    refill_volume: float  # mL to add
    ingredient_name: str

    def get_operation_effects(self):
        chem = Chemical.make_up_from_smiles("O")  # placeholder
        chem.quantify_by_volume(self.refill_volume)
        chem.smiles = self.ingredient_name
        pom = PortionOfMaterial()
        pom.add_chemical(chem)

        return [
            Create(instance_1_iri=pom.identifier),
            AddObjectProperty(
                instance_1_iri=pom.identifier,
                instance_2_iri=self.participant_reservoir,
                property_iri=Is_directly_contained_by.predicate_iri
            )
        ]


def init_world(fbs: FlavorBlendSim):
    # create platform components
    platform = LabObject(identifier="gx-281", is_present={True, })

    liquid_transfer_device = LiquidTransferDevice(identifier="liquid_transfer_device", is_present={True, }, )
    liquid_transfer_device.is_directly_contained_by.add(platform)

    wash_station = LabObject(identifier="wash_station", is_present={True, }, )
    wash_station.is_directly_contained_by.add(platform)

    # create reservoirs
    reservoirs = []
    for ingredient_name, reservoir_iri in fbs.reservoir_map.items():
        reservoir = MaterialContainer(identifier=reservoir_iri, is_present={True, }, )
        chem = Chemical.make_up_from_smiles("O")
        chem.quantify_by_volume(fbs.v_reservoir_capacity)
        chem.smiles = ingredient_name
        pom = PortionOfMaterial(is_present={True, })
        pom.add_chemical(chem)
        pom.is_directly_contained_by.add(reservoir)
        reservoirs.append(reservoir)

    # create vials
    vials = []
    terpene_vials = []
    flavor_vials = []
    for vial_iri, formulation_name in fbs.vial_formulation_map.items():
        vial = MaterialContainer(identifier=vial_iri, is_present={True, }, )
        vials.append(vial)
        if formulation_name.startswith("TERPENE_FORMULATION"):
            terpene_vials.append(vial)
        elif formulation_name.startswith("FLAVOR_FORMULATION"):
            flavor_vials.append(vial)
        else:
            raise ValueError("formulation name must start with TERPENE_FORMULATION or FLAVOR_FORMULATION")

    # create racks
    n_placed_vial_total = 0
    terpene_racks = []
    for i in range(fbs.n_terpene_rack):
        rack = LabObject(identifier=f"rack_terpene_{i}", is_present={True, }, )
        rack.is_directly_contained_by.add(platform)
        n_placed_vial = 0
        while n_placed_vial < fbs.n_vials_per_rack:
            v = terpene_vials.pop()  # raise on empty
            v.is_directly_contained_by.add(rack)
            n_placed_vial += 1
        n_placed_vial_total += n_placed_vial
        terpene_racks.append(rack)
    assert n_placed_vial_total == fbs.n_terpene_rack * fbs.n_vials_per_rack

    n_placed_vial_total = 0
    flavor_racks = []
    for i in range(fbs.n_flavor_rack):
        rack = LabObject(identifier=f"rack_flavor_{i}", is_present={True, }, )
        rack.is_directly_contained_by.add(platform)
        n_placed_vial = 0
        while n_placed_vial < fbs.n_vials_per_rack:
            v = flavor_vials.pop()
            v.is_directly_contained_by.add(rack)
            n_placed_vial += 1
        n_placed_vial_total += n_placed_vial
        flavor_racks.append(rack)
    assert n_placed_vial_total == fbs.n_flavor_rack * fbs.n_vials_per_rack

    return platform, liquid_transfer_device, wash_station, reservoirs, vials, flavor_racks, terpene_racks


def build_ops(
        fbs: FlavorBlendSim,
        platform: LabObject,
        transfer_device: LabObject,
        wash_station: LabObject,
        vials: list[MaterialContainer],
) -> tuple[list[Operation], dict[str, int | Any]]:
    """
    Serial schedule: for ing₁ → all vials, WASH, ing₂ → …, WASH …
    with refills inserted.
    """
    ops: list[Operation] = []
    res_level = {i: fbs.v_reservoir_capacity for i in fbs.reservoir_map}
    last_op: Operation | None = None

    t_xfer_terpene = 0
    t_xfer_flavor = 0
    t_refill = 0
    t_wash = 0

    rng = np.random.RandomState(fbs.seed)
    def sample_liquid_transfer_time() -> float:
        """
        Returns a positive transfer duration.

        * 95 % of the time:  ~N(BASE_MU, BASE_SD)
        *  5 % of the time:  that base  +  Exp(mean=C_DELAY/P_GLITCH)
          (≈ occasional jam / re-try etc.)
        """
        base_mu = fbs.t_single_liquid_transfer_off_bed  # e.g. 9 s
        base_sd = 0.2 * base_mu  # mild normal jitter
        p_glitch = 0.05
        t = max(rng.normal(base_mu, base_sd), 0.1)

        if rng.random() < p_glitch:
            max_extra = 30.0
            lam = 1 / ((4 * base_mu) / p_glitch)  # example scale for truncated exp
            cdf_max = 1 - np.exp(-lam * max_extra)
            u = rng.random()
            extra = -np.log(1 - u * cdf_max) / lam
            t += extra
        return t


    # per ingredient transfer, saves on washing cost
    for ingredient_name, reservoir_iri in fbs.reservoir_map.items():
        for vial in vials:
            formulation_name = fbs.vial_formulation_map[vial.instance_iri]
            assert formulation_name in fbs.formulations, f"formulation not defined: {formulation_name}"
            try:
                vol = fbs.formulations[formulation_name][ingredient_name]
            except KeyError:
                continue

            # add refill
            if res_level[ingredient_name] - vol < fbs.v_reservoir_capacity * fbs.p_reservoir_refill_threshold / 100:
                refill = ReservoirRefill(
                    participant_transfer_device=transfer_device.instance_iri,
                    participant_reservoir=reservoir_iri,
                    participant_platform=platform.instance_iri,
                    refill_volume=fbs.v_reservoir_capacity - res_level[ingredient_name],
                    temporal_cost=fbs.t_reservoir_refill,
                    ingredient_name=ingredient_name
                )
                if last_op:
                    refill.required_precedents = [last_op.identifier]
                t_refill += refill.temporal_cost
                ops.append(refill)
                last_op = refill
                res_level[ingredient_name] = fbs.v_reservoir_capacity

            # actual transfer
            # TODO set "is_clean" to {False,} in transfer
            xfer = TransferMaterialByVolume(
                participant_source=reservoir_iri,
                participant_destination=vial.identifier,
                participant_device=transfer_device.instance_iri,
                transfer_volume=vol,
                # temporal_cost=fbs.t_single_liquid_transfer_off_bed,
                temporal_cost=sample_liquid_transfer_time(),
            )
            if ingredient_name.startswith("ING_TERPENE"):
                t_xfer_terpene += xfer.temporal_cost
            elif ingredient_name.startswith("ING_FLAVOR"):
                t_xfer_flavor += xfer.temporal_cost
            if last_op:
                xfer.required_precedents = [last_op.identifier]
            ops.append(xfer)
            last_op = xfer
            res_level[ingredient_name] -= vol

        if last_op:
            wash = Wash(
                participant_transfer_device=transfer_device.instance_iri,
                participant_wash_station=wash_station.identifier,
                temporal_cost=fbs.t_wash,
            )
            wash.required_precedents = [last_op.identifier]
            t_wash += wash.temporal_cost
            ops.append(wash)
            last_op = wash
    time_estimates = {
        "t_xfer_terpene": t_xfer_terpene,
        "t_xfer_flavor": t_xfer_flavor,
        "t_refill": t_refill,
        "t_wash": t_wash,
    }
    logger.warning(pformat(time_estimates))
    return ops, time_estimates

def throughput_estimate():
    data = []
    n_runs = 10
    for t_single_mean in range(7, 14):
        for i_run in range(1, n_runs+1):
            fbs = FlavorBlendSim.random_init()
            fbs.seed = i_run
            fbs.t_single_liquid_transfer_off_bed = t_single_mean
            platform, liquid_transfer_device, wash_station, reservoirs, vials, flavor_racks, terpene_racks = init_world(fbs)
            ops, estimates = build_ops(fbs, platform, liquid_transfer_device, wash_station, vials)
            KnowledgeGraph.clear_object_lookup()
            makespan = sum(estimates.values()) / 3600
            makespan_flavor = makespan - estimates["t_xfer_terpene"] / 3600  # this is an overestimate as washing includes the other type of ingredient
            makespan_terpene = makespan - estimates["t_xfer_flavor"] / 3600
            throughput_terpene = fbs.n_terpene_rack * fbs.n_vials_per_rack / makespan_terpene
            throughput_flavor = fbs.n_flavor_rack * fbs.n_vials_per_rack / makespan_flavor
            throughput = len(fbs.vial_formulation_map) / makespan
            data.append(
                {
                    "t_single_mean": t_single_mean,
                    "throughput_terpene": throughput_terpene,
                    "throughput_flavor": throughput_flavor,
                    "throughput": throughput,
                    "makespan": makespan
                }
            )

    df = pd.DataFrame(data)
    groups = df.groupby("t_single_mean")

    records = []
    for t, sub in groups:
        v_terp = sub["throughput_terpene"].to_numpy()
        v_flav = sub["throughput_flavor"].to_numpy()
        v_all = sub["throughput"].to_numpy()

        ci_terp = bootstrap((v_terp,), np.mean, confidence_level=0.95,
                            n_resamples=10_000, vectorized=False).confidence_interval
        ci_flav = bootstrap((v_flav,), np.mean, confidence_level=0.95,
                            n_resamples=10_000, vectorized=False).confidence_interval
        ci_all = bootstrap((v_all,), np.mean, confidence_level=0.95, n_resamples=10_000, vectorized=False).confidence_interval

        records.append(dict(
            t_single_mean=t,
            mean_terpene=v_terp.mean(),
            low_terpene=ci_terp.low,
            up_terpene=ci_terp.high,
            mean_flavor=v_flav.mean(),
            low_flavor=ci_flav.low,
            up_flavor=ci_flav.high,
            mean_all=v_all.mean(),
            low_all=ci_all.low,
            up_all=ci_all.high,
        ))

    agg = pd.DataFrame(records).sort_values("t_single_mean")
    agg.to_csv("throughput_estimate.csv", index=False)

    # Build asymmetric yerr arrays
    yerr_terp = np.vstack([agg["mean_terpene"] - agg["low_terpene"],
                           agg["up_terpene"] - agg["mean_terpene"]])
    yerr_flav = np.vstack([agg["mean_flavor"] - agg["low_flavor"],
                           agg["up_flavor"] - agg["mean_flavor"]])
    yerr_all = np.vstack([agg["mean_all"] - agg["low_all"],
                          agg["up_all"] - agg["mean_all"]])

    fig, ax = plt.subplots(figsize=(6, 4))

    ax.errorbar(agg["t_single_mean"], agg["mean_terpene"],
                yerr=yerr_terp, fmt="-o", color="red", label=f"30-component terpene blends")

    ax.errorbar(agg["t_single_mean"], agg["mean_flavor"],
                yerr=yerr_flav, fmt="--o", color="blue", label="20-component flavor blends")

    # ax.errorbar(agg["t_single_mean"], agg["mean_all"], yerr=yerr_all, fmt=":o", color="black", label="All blends")
    ax.set_xlabel("Temporal cost per single dispense (seconds)", fontsize=14)
    ax.set_ylabel("Throughput (blends per hour)", fontsize=14)
    ax.legend(fontsize=14)
    plt.tight_layout()
    fig.savefig("throughput_estimate.png", dpi=600)
    return df






def main():
    fbs = FlavorBlendSim.random_init()
    platform, liquid_transfer_device, wash_station, reservoirs, vials, flavor_racks, terpene_racks = init_world(fbs)
    ops, _ = build_ops(fbs, platform, liquid_transfer_device, wash_station, vials)
    sim = Simulation(operations=ops)
    sim.run()

    out_dir = Path("sim_outputs")
    out_dir.mkdir(exist_ok=True)
    sim.export_event_log(out_dir / "event_log_granular.csv")
    KnowledgeGraph.graph().serialize(out_dir / "state_after.ttl", format="turtle")

    logger.info(f"Finished @ t = {sim.env.now:.1f} s – event log & KG written to {out_dir}/")


if __name__ == "__main__":
    throughput_estimate()
    # main()
