from __future__ import annotations

import os.path

from libsyn_tools.utils import *
from .network_reaction import ReactionNetwork
from .reaction import ChemicalReaction


class NetworkWorkflow(BaseModel):
    routes_file: FilePath

    work_folder: FilePath

    scraper_output: FilePath | None = None

    query_askcos: bool = False

    sample_seed: int | None = None

    sample_n_target: int | None = None

    network_json: FilePath = "reaction_network.json"

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
        json_dump(reaction_network.model_dump(), os.path.join(self.work_folder, self.network_json))
