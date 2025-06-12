import math
import os.path

import numpy as np
import pandas as pd
from loguru import logger

from libsyn_tools.chem_schema import Chemical, OperationType, ReactionNetwork
from libsyn_tools.opt.schema import Solver, OperationNetwork
from libsyn_tools.sim import LabObject, PortionOfMaterial, UnitaryEdit, UnitaryEditType, \
    TransferMaterialByPortionSize, ActionSimulation, ConvertByHeating, Evaporate, KnowledgeGraph
from libsyn_tools.utils import json_load, FilePath

"""
Load a scheduler instance and simulate it with or without operation time sampling.
"""
np.random.seed(42)


def load_schedule_instance(scheduler_instance_folder: FilePath = "../../workplace_opt/RUNS/FDA-03-09-0-0"):
    solver_file = os.path.join(scheduler_instance_folder, "solver_milp.json")
    operation_network_file = os.path.join(scheduler_instance_folder, "operation_network.json")
    reaction_network_file = os.path.join(scheduler_instance_folder, "reaction_network.json")

    solver = Solver(**json_load(solver_file))
    operation_network = OperationNetwork(**json_load(operation_network_file))
    reaction_network = ReactionNetwork(**json_load(reaction_network_file))
    return solver, operation_network, reaction_network, os.path.basename(scheduler_instance_folder)


def init_world():
    solver, operation_network, reaction_network, schedule_instance_name = load_schedule_instance()
    scheduler_makespan = max(solver.output.end_times.values())
    logger.info(f"makespan from the scheduler: {scheduler_makespan}")
    reaction_vessels = dict()
    inventory_smi2moles = dict()
    for r in reaction_network.chemical_reactions:
        vessel = LabObject(identifier=f"ReactionVessel-{r.identifier}")
        UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=vessel.instance_iri).apply()
        reaction_vessels[r.identifier] = vessel
        for c in r.reactants + r.reagents:
            if c.smiles not in inventory_smi2moles:
                inventory_smi2moles[c.smiles] = c.moles
            else:
                inventory_smi2moles[c.smiles] += c.moles

    chemical_storage = LabObject(identifier="SimplifiedChemicalStorage")
    chemical = Chemical.make_up_from_smiles("O")
    chemical.quantify_by_moles(1)
    pom = PortionOfMaterial(identifier=f"InitialStorage")
    pom.add_chemical(chemical)
    pom.is_directly_contained_by.add(chemical_storage)
    UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=pom.instance_iri).apply()

    UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=chemical_storage.instance_iri).apply()
    # TODO an all-encompassing storage, criminally simplified
    # chemical_storages = dict()
    # for smi, moles in inventory_smi2moles.items():
    #     storage = LabObject(identifier=f"ChemicalStorage:{smi}")
    #     UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=storage.instance_iri).apply()
    #     chemical_storages[smi] = storage
    #     chemical = Chemical.make_up_from_smiles(smi)
    #     chemical.quantify_by_moles(moles)
    #     pom = PortionOfMaterial(identifier=f"InitialStorage:{smi}")
    #     pom.add_chemical(chemical)
    #     pom.is_directly_contained_by.add(storage)
    #     UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=pom.instance_iri).apply()

    fms = dict()
    for fm in solver.input.functional_modules:
        fm_lab_object = LabObject(identifier=fm.identifier)
        fms[fm.identifier] = fm_lab_object
        UnitaryEdit(type=UnitaryEditType.CREATE, instance_1_iri=fm_lab_object.instance_iri).apply()
    return operation_network, solver, reaction_network, fms, reaction_vessels, chemical_storage, schedule_instance_name


def weibull_shape_from_cv(target_cv, tol=1e-6, low=0.1, high=100):
    """
    Numerically invert the CV relation for the Weibull distribution:
      CV(k) = sqrt(Γ(1+2/k)/(Γ(1+1/k)^2) - 1)
    to find k such that CV(k) ≈ target_cv.
    """
    # Binary search for k in the interval [low, high]
    while high - low > tol:
        mid = (low + high) / 2
        cv_mid = math.sqrt(math.gamma(1 + 2 / mid) / (math.gamma(1 + 1 / mid) ** 2) - 1)
        if cv_mid > target_cv:
            # If cv_mid is too high, increase k to lower the CV.
            low = mid
        else:
            high = mid
    return (low + high) / 2


def sample_operation_time(dist_type, x, dist_param):
    """
    Returns a random sample from the specified distribution with mean x and
    coefficient of variation (CV) equal to dist_param.

    Parameters:
      dist_type (str): The type of distribution. One of:
          'normal'    - Uses sigma = dist_param * x.
          'lognormal' - Uses sigma_ln computed so that CV = sqrt(exp(sigma_ln^2)-1) = dist_param.
          'exponential'- The exponential distribution has CV = 1; dist_param must equal 1.
          'gamma'     - Uses shape k = 1/(dist_param^2) and theta = x/k.
          'weibull'   - Numerically finds shape k such that the Weibull CV equals dist_param.
      x (float): The desired mean operation time.
      dist_param (float): The target coefficient of variation (CV = std/mean).

    Returns:
      float: A random sample from the specified distribution.
    """
    dist_type = dist_type.lower()

    if dist_type == 'normal':
        # For the normal distribution, simply set sigma = dist_param * x.
        sigma = dist_param * x
        sample = np.random.normal(loc=x, scale=sigma)
        # Ensure nonnegative operation time.
        while sample < 0:
            sample = np.random.normal(loc=x, scale=sigma)
        return sample

    elif dist_type == 'lognormal':
        # For a lognormal, we have:
        #   CV = sqrt(exp(sigma_ln^2)-1)
        # Solve for sigma_ln: sigma_ln = sqrt(log(dist_param^2 + 1))
        sigma_ln = math.sqrt(math.log(dist_param ** 2 + 1))
        # To have mean = x, we set:
        #   x = exp(mu_ln + sigma_ln^2/2)  -->  mu_ln = log(x) - sigma_ln^2/2
        mu_ln = math.log(x) - (sigma_ln ** 2) / 2
        return np.random.lognormal(mean=mu_ln, sigma=sigma_ln)

    elif dist_type == 'exponential':
        # The exponential distribution has a fixed CV = 1.
        if abs(dist_param - 1) > 1e-6:
            raise ValueError("Exponential distribution always has a CV of 1. Use dist_param=1.")
        return np.random.exponential(scale=x)

    elif dist_type == 'gamma':
        # For the gamma distribution, if T ~ Gamma(k, theta) then:
        #   mean = k * theta and variance = k * theta^2, so CV = 1/sqrt(k).
        # To have CV = dist_param, set: 1/sqrt(k) = dist_param  -->  k = 1/(dist_param^2)
        k = 1 / (dist_param ** 2)
        theta = x / k  # so that mean = k * theta = x.
        return np.random.gamma(shape=k, scale=theta)

    elif dist_type == 'weibull':
        # For the Weibull distribution, let the shape be k and the scale be λ.
        # The mean is λ * Gamma(1 + 1/k) and the CV is given by:
        #   CV = sqrt(Γ(1+2/k)/(Γ(1+1/k)^2) - 1)
        # We numerically solve for k such that the CV equals dist_param.
        k_weibull = weibull_shape_from_cv(dist_param)
        # Then, to have mean = x, set λ = x / Gamma(1 + 1/k).
        lambda_weibull = x / math.gamma(1 + 1 / k_weibull)
        return np.random.weibull(a=k_weibull) * lambda_weibull

    elif dist_type == 'exact':
        return x

    else:
        raise ValueError(
            "Unsupported distribution type. Choose from 'exact', 'normal', 'lognormal', 'exponential', 'gamma', or 'weibull'.")


def run_sim_once(
        operation_network, solver, reaction_network,
        fms, reaction_vessels, chemical_storage,
        dist_type, dist_cv
):
    """
    Simulate the library synthesis using a given scheduling instance.

    :param operation_network: the network built from operations and the precedence relations among them
    :param solver: the solver with solved schedule
    :param reaction_network: the corresponding reaction network
    :param fms: the functional modules
    :param reaction_vessels: a dictionary of {<reaction id>: <LabObject>}
    :param chemical_storage: the simplified chemical storage as a `LabObject`
    :param dist_type: the type of distribution, see `sample_operation_time`
    :param dist_cv: the desired coefficient of variation for the selected distribution
    """
    use_scheduler = True
    actions = []
    for o in operation_network.operations:
        oid = o.identifier
        fmid = solver.output.assignments[oid]
        time_cost = solver.output.end_times[oid] - solver.output.start_times[oid]
        time_cost = sample_operation_time(dist_type=dist_type, x=time_cost, dist_param=dist_cv)
        reaction_id = o.from_reaction
        reaction = reaction_network.entity_dictionary[reaction_id]
        reaction_vessel = reaction_vessels[reaction_id]
        if o.type in (OperationType.TransferLiquid, OperationType.TransferSolid, OperationType.MakeSolution):
            action = TransferMaterialByPortionSize(
                identifier=o.identifier,
                temporal_cost=time_cost,
                scheduled_start_time=solver.output.start_times[o.identifier] if use_scheduler else 0,
                required_precedents=o.precedents,
                action_effects_description=str(o.annotations),
                source_iri=chemical_storage.instance_iri,
                destination_iri=reaction_vessel.instance_iri,
                transfer_device_iri=fms[fmid].instance_iri,
                portion_size=0.5  # TODO this is a dummy value for now
            )
        elif o.type == OperationType.Heating:
            action = ConvertByHeating(
                identifier=o.identifier,
                temporal_cost=time_cost,
                scheduled_start_time=solver.output.start_times[o.identifier] if use_scheduler else 0,
                required_precedents=o.precedents,
                action_effects_description=str(o.annotations),
                container_iri=reaction_vessel.instance_iri,
                heating_module_iri=fms[fmid].instance_iri,
                conversion_ratio=1.0,
                reaction=reaction,
            )
        elif o.type == OperationType.ConcentrationAndPurification:
            action = Evaporate(
                identifier=o.identifier,
                temporal_cost=time_cost,
                scheduled_start_time=solver.output.start_times[o.identifier] if use_scheduler else 0,
                required_precedents=o.precedents,
                action_effects_description=str(o.annotations),
                container_iri=reaction_vessel.instance_iri,
                evaporation_module_iri=fms[fmid].instance_iri,
                reaction=reaction,
            )
        else:
            raise ValueError(f"Operation type {o.type} not supported")
        actions.append(action)

    sim = ActionSimulation(actions=actions)

    sim.run()

    return sim, sim.history_log[-1].timestamp


def run_all(export_sim=False, test=False):
    if test:
        export_sim = True
    operation_network, solver, reaction_network, fms, reaction_vessels, chemical_storage, schedule_instance_name = init_world()
    results = []
    for dist_type in [
        'normal',
        'gamma',
        'weibull'
    ]:
        for dist_cv in [0.1, 0.2, 0.3, 0.4, 0.5]:
            for i_repeat in range(10):
                sim, makespan = run_sim_once(operation_network, solver, reaction_network, fms, reaction_vessels,
                                        chemical_storage, dist_type, dist_cv)
                results.append(
                    {
                        "dist_type": dist_type,
                        "dist_cv": dist_cv,
                        "i_repeat": i_repeat,
                        "makespan": makespan,
                    }
                )

                if export_sim:
                    out_file = f"{schedule_instance_name}_{dist_type}_{dist_cv}_{i_repeat}.ttl"
                    KnowledgeGraph.graph().serialize(destination=out_file, format="turtle")
                    logger.info(f"Exported updated knowledge graph to {out_file}")
                    sim.export_event_log(out_file.replace(".ttl", ".csv"))
                if test:
                    break
            if test:
                break

    pd.DataFrame.from_records(results).to_csv(__file__.replace(".py", ".csv"), index=False)


if __name__ == '__main__':
    run_all(test=True)
