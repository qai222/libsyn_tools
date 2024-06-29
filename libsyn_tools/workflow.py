from __future__ import annotations

import os.path

from loguru import logger

from libsyn_tools.chem_schema.network_operation import OperationNetwork
from libsyn_tools.chem_schema.network_reaction import ReactionNetwork
from libsyn_tools.chem_schema.reaction import ChemicalReaction
from libsyn_tools.opt import random_functional_modules, SchedulerInput, SolverMILP
from libsyn_tools.utils import *


class Workflow(BaseModel):
    routes_file: FilePath

    work_folder: FilePath

    scraper_output: FilePath | None = None

    query_askcos: bool = False

    sample_seed: int | None = None

    sample_n_target: int | None = None

    reaction_network_json: FilePath = "reaction_network.json"

    operation_network_json: FilePath = "operation_network.json"

    scheduler_input_json: FilePath = "scheduler_input.json"

    scheduler_output_json: FilePath = "scheduler_output.json"

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
            os.path.join(self.work_folder, "scraper_input.json")
        )

    def export_reaction_network(self, dummy_quantify=False):
        """
        export the reaction network
        """
        reaction_smiles = parse_sparrow_routes(self.routes_file, self.sample_seed, self.sample_n_target)
        reactions = [
            ChemicalReaction.from_reaction_smiles(smi, self.query_askcos, self.scraper_output)
            for smi in reaction_smiles
        ]
        reaction_network = ReactionNetwork.from_reactions(reactions)
        if dummy_quantify:
            reaction_network.dummy_quantify(by_mass=False)
        json_dump(reaction_network.model_dump(), os.path.join(self.work_folder, self.reaction_network_json))

    def export_operation_network(self):
        reaction_network = json_load(os.path.join(self.work_folder, self.reaction_network_json))
        reaction_network = ReactionNetwork(**reaction_network)
        for r in reaction_network.chemical_reactions:
            if len([c for c in r.reactants + r.reagents if c.state_of_matter == StateOfMatter.LIQUID]) == 0:
                r.reactants[0].state_of_matter = StateOfMatter.LIQUID
        operation_network = OperationNetwork.from_reaction_network(reaction_network)
        logger.info(f"# of operations: {len(operation_network.operations)}")
        json_dump(operation_network.model_dump(), os.path.join(self.work_folder, self.operation_network_json))

    def export_scheduler(self, max_capacity=1, max_module_number_per_type=1, temperature_threshold=50,
                         dummy_work_shifts=False, ):

        operation_network = json_load(os.path.join(self.work_folder, self.operation_network_json))
        operation_network = OperationNetwork(**operation_network)

        functional_modules = random_functional_modules(
            max_capacity=max_capacity,
            max_module_number_per_type=max_module_number_per_type,
        )
        required_operation_types = set([o.type for o in operation_network.operations])
        functional_modules = [fm for fm in functional_modules if set(fm.can_process).issubset(required_operation_types)]
        si = SchedulerInput.build_input(
            operation_network, functional_modules, temperature_threshold=temperature_threshold
        )
        if dummy_work_shifts:
            work_shifts = si.get_dummy_work_shifts()
        else:
            work_shifts = None
        si.frak_W = work_shifts
        solver = SolverMILP(input=si)
        solver.solve()
        json_dump(solver.input.model_dump(), os.path.join(self.work_folder, self.scheduler_input_json))
        json_dump(solver.output.model_dump(), os.path.join(self.work_folder, self.scheduler_output_json))
