from __future__ import annotations

import os.path

from loguru import logger

from libsyn_tools.chem_schema import FunctionalModule, OperationType
from libsyn_tools.chem_schema.network_operation import OperationNetwork
from libsyn_tools.chem_schema.network_reaction import ReactionNetwork
from libsyn_tools.chem_schema.reaction import ChemicalReaction
from libsyn_tools.opt import random_functional_modules, SchedulerInput, SolverMILP, SolverBaseline
from libsyn_tools.utils import *


class Workflow(BaseModel):
    routes_file: FilePath

    work_folder: FilePath

    scraper_output: FilePath | None = None

    scraper_input: FilePath = "scraper_input.json"

    functional_modules_json: FilePath = "functional_modules.json"

    reaction_network_json: FilePath = "reaction_network.json"

    operation_network_json: FilePath = "operation_network.json"

    scheduler_input_json: FilePath = "scheduler_input.json"

    solver_milp_json: FilePath = "solver_milp.json"

    solver_baseline_json: FilePath = "solver_baseline.json"

    def create_scraper_input(self):
        """
        collect more info based on sparrow's output. this produces the input for `ChemScraper`.
        skip this if you are fine with made up physical properties
        """
        reaction_smiles = parse_sparrow_routes(routes_file=self.routes_file)
        reactions = [ChemicalReaction.from_reaction_smiles(smi) for smi in reaction_smiles]
        reaction_network = ReactionNetwork.from_reactions(reactions)
        json_dump(
            sorted(reaction_network.molecular_smiles_table.keys()),
            os.path.join(self.work_folder, self.scraper_input)
        )

    def export_reaction_network(self, target_sample_seed: int = None, query_askcos: bool = False,
                                sample_n_target: int = None, dummy_quantify=False, specified_targets=None):
        """
        export the reaction network
        """
        reaction_smiles = parse_sparrow_routes(self.routes_file, target_sample_seed, sample_n_target, specified_targets)
        reactions = [
            ChemicalReaction.from_reaction_smiles(smi, query_askcos, self.scraper_output)
            for smi in reaction_smiles
        ]
        reaction_network = ReactionNetwork.from_reactions(reactions, ignore_multi_prod=True)
        if sample_n_target is not None:
            assert sample_n_target <= len(reaction_network.target_smiles)
        if dummy_quantify:
            reaction_network.dummy_quantify(by_mass=False)
        json_dump(reaction_network.model_dump(), os.path.join(self.work_folder, self.reaction_network_json))

    def export_operation_network(self, rng: random.Random):
        reaction_network = json_load(os.path.join(self.work_folder, self.reaction_network_json))
        reaction_network = ReactionNetwork(**reaction_network)
        reaction_network.dummy_quantify(by_mass=False)
        for r in reaction_network.chemical_reactions:
            if len([c for c in r.reactants + r.reagents if c.state_of_matter == StateOfMatter.LIQUID]) == 0:
                r.reactants[0].state_of_matter = StateOfMatter.LIQUID
        operation_network = OperationNetwork.from_reaction_network(reaction_network, rng)
        logger.info(f"# of operations: {len(operation_network.operations)}")
        json_dump(operation_network.model_dump(), os.path.join(self.work_folder, self.operation_network_json))

    def export_functional_modules(
            self, rng: random.Random = None,
            max_capacity=1, max_module_number_per_type=1, random=True,
            module_settings: dict[OperationType, tuple[int, int]] = None,
    ):
        if random:
            assert module_settings is None
            fms = random_functional_modules(
                rng,
                max_capacity=max_capacity,
                max_module_number_per_type=max_module_number_per_type,
            )
        else:
            assert module_settings
            assert not random
            fms = []
            for t in OperationType:
                try:
                    capacity, module_number = module_settings[t]
                except KeyError:
                    capacity = 1
                    module_number = 1
                for i in range(module_number):
                    m = FunctionalModule(
                        name=f"{t}-{i}",
                        can_process=[t, ],
                        capacity=capacity,
                    )
                    fms.append(m)
        json_dump([fm.model_dump() for fm in fms], os.path.join(self.work_folder, self.functional_modules_json))

    def export_scheduler_input(self, rng: random.Random,
                               # temperature_threshold,
                               dummy_work_shifts=False):

        operation_network = json_load(os.path.join(self.work_folder, self.operation_network_json))
        operation_network = OperationNetwork(**operation_network)

        functional_modules = json_load(os.path.join(self.work_folder, self.functional_modules_json))
        functional_modules = [FunctionalModule(**fm) for fm in functional_modules]

        required_operation_types = set([o.type for o in operation_network.operations])
        functional_modules = [
            fm for fm in functional_modules if set(fm.can_process).issubset(required_operation_types)
        ]
        si = SchedulerInput.build_input(
            rng,
            operation_network, functional_modules,
            # temperature_threshold=temperature_threshold
        )
        # logger.info(f"{operation_network.max_processing_times}")
        if dummy_work_shifts:
            work_shifts = si.get_dummy_work_shifts()
        else:
            work_shifts = None
        si.frak_W = work_shifts

        json_dump(si.model_dump(), os.path.join(self.work_folder, self.scheduler_input_json))

    def export_solver(self, baseline, time_limit, gb_threads):
        si = json_load(os.path.join(self.work_folder, self.scheduler_input_json))
        si = SchedulerInput(**si)

        if baseline:
            solver = SolverBaseline(
                input=si,
                operation_network=OperationNetwork(
                    **json_load(os.path.join(self.work_folder, self.operation_network_json))),
                reaction_network=ReactionNetwork(
                    **json_load(os.path.join(self.work_folder, self.reaction_network_json)))
            )
            solver.solve()
            solver.output.notes['validation'] = solver.output.validate_schedule(solver.input)
        else:
            solver = SolverMILP(input=si, time_limit=time_limit)
            solver.solve(logfile=os.path.join(self.work_folder, "gurobi.log"), threads=gb_threads)
            if solver.opt_log['gurobi status'] == 3 or (solver.opt_log['gurobi status'] == 9 and solver.opt_log[
                'gurobi solution count'] == 0):  # infeasible model or zero solution
                solver.output.notes['validation'] = {}
            else:
                solver.output.notes['validation'] = solver.output.validate_schedule(solver.input)

        if baseline:
            json_dump(solver.model_dump(), os.path.join(self.work_folder, self.solver_baseline_json))
        else:
            json_dump(solver.model_dump(), os.path.join(self.work_folder, self.solver_milp_json))
