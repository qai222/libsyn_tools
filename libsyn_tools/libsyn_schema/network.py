from __future__ import annotations

import functools
from collections import defaultdict

import networkx as nx

from .base import Entity
from .chemical import Chemical
from .reaction import ChemicalReaction


class ReactionNetwork(Entity):
    """
    the bipartite graph of SMILES and `ChemicalReaction`.
    note the chemicals in this graph are indexed by SMILES,
    two chemicals of the same composition but different forms will be grouped to the same node
    """

    chemical_reactions: list[ChemicalReaction]
    """ all chemical reactions in this network """

    molecular_nodes: list[str] = []
    """ a list of molecular SMILES """

    edges: list[tuple[str, str]] = []
    """ directed edges, can be either
    [molecular SMILES, reaction identifier] or [reaction identifier, molecular SMILES]  """

    @property
    def nx_digraph(self) -> nx.DiGraph:
        """
        convert to a networkx directed graph, nodes are **identifiers**
        """
        g = nx.DiGraph()
        for node in self.molecular_nodes:
            g.add_node(node)
        for cr in self.chemical_reactions:
            g.add_node(cr.identifier)
        for u, v in self.edges:
            g.add_edge(u, v)
        return g

    @functools.cached_property
    def entity_dictionary(self) -> dict[str, ChemicalReaction | Chemical]:
        """ a lookup table for all chemicals and reactions in this network """
        d = dict()
        for reaction in self.chemical_reactions:
            d[reaction.identifier] = reaction
            for chemical in reaction.chemicals:
                d[chemical.identifier] = chemical
        return d

    @property
    def molecular_smiles_table(self) -> dict[str, list[Chemical]]:
        """ SMILES -> a list of chemicals that have the SMILES"""
        d = defaultdict(list)
        for eid, e in self.entity_dictionary.items():
            if isinstance(e, Chemical):
                d[e.smiles].append(e)
        return d

    @classmethod
    def from_reactions(cls, reactions: list[ChemicalReaction]) -> ReactionNetwork:
        molecular_smiles = []
        edges = []
        for reaction in reactions:
            for chemical in reaction.reactants:
                molecular_smiles.append(chemical.smiles)
                assert len(chemical.is_consumed_by) + len(chemical.is_produced_by) == 1
                if chemical.is_produced_by:
                    edge = (reaction.identifier, chemical.smiles)
                else:
                    edge = (chemical.smiles, reaction.identifier)
                edges.append(edge)
        molecular_smiles = sorted(set(molecular_smiles))
        edges = sorted(set(edges))
        return cls(chemical_reactions=reactions, molecular_nodes=molecular_smiles, edges=edges)

    def quantify(self, targets: dict[str, float], intended_ratios: dict[str, dict[str, float] | None],
                 expected_yields: dict[str, float | None]) -> None:
        """
        quantify the reactions in this network. all reactions have to be uni-product.
        (for multi-product reactions, one could simply modify the reaction SMILES if there is only one desired product.)
        for every molecular SMILES, there is at most one reaction has it as the product.

        :param expected_yields: the yield table for uni-product reactions, using reaction identifiers as the domain.
        :param intended_ratios: the intended ratios for uni-product reactions, using reaction identifiers as the domain
        :param targets: a map from molecular SMILES to its expected mass in grams, every one of these molecules can only
        be made from one reaction
        :return:
        """
        # lookup table for product -> reaction
        product_smi_to_reaction = dict()
        for cr in self.chemical_reactions:
            assert cr.is_uniproduct, "network can only be quantified when all reactions are uni-product"
            product_smi_to_reaction[cr.product_smiles] = cr

        # create a tree from the reversed network to run bfs
        g = self.nx_digraph.reverse(copy=True)
        dummy_root = "DUMMY_ROOT"
        for target_smi in targets:
            g.add_edge(dummy_root, target_smi)
        bfs_edges = nx.bfs_edges(g, dummy_root)
        bfs_nodes = [dummy_root] + [v for u, v in bfs_edges]

        # (intermediate) product SMILES -> mass
        required_product_quantities = {k: v for k, v in targets.items()}  # make a copy

        # quantify the reaction network
        for n in bfs_nodes:
            # skip root and molecular smiles
            if n == dummy_root:
                continue
            if n in self.molecular_smiles_table:
                continue
            # quantify the reaction
            reaction = self.entity_dictionary[n]
            reaction: ChemicalReaction
            reaction.quantify(
                product_smiles=reaction.product_smiles,
                product_mass=required_product_quantities[reaction.product_smiles],
                intended_ratios=intended_ratios[reaction.identifier],
                expected_yield=expected_yields[reaction.identifier]
            )
            for chemical in reaction.reagents + reaction.reactants:
                required_product_quantities[chemical.smiles] = chemical.mass
