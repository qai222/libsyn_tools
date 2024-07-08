from __future__ import annotations

import functools
from collections import defaultdict

import networkx as nx

from .base import Entity
from .chemical import Chemical
from .reaction import ChemicalReaction
from ..utils import is_uuid, calculate_mw, CytoEdge, CytoEdgeData, CytoNodeData, CytoNode, drawing_url, StateOfMatter


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
            assert u in g.nodes and v in g.nodes
            g.add_edge(u, v)
        return g

    @property
    def target_smiles(self) -> list[str]:
        """
        a list of target SMILES, obtained from the directed graph
        """
        smis = []
        g = self.nx_digraph
        for node in g.nodes:
            if g.out_degree(node) == 0 and g.in_degree(node) >= 1:
                assert not is_uuid(node)
                smis.append(node)
        return smis

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
    def from_reactions(cls, reactions: list[ChemicalReaction], ignore_multi_prod:bool=False) -> ReactionNetwork:
        molecular_smiles = []
        edges = []
        network_reactions = []
        for reaction in reactions:
            if ignore_multi_prod and not reaction.is_uniproduct:
                continue
            for chemical in reaction.chemicals:
                molecular_smiles.append(chemical.smiles)
                assert len(chemical.is_consumed_by) + len(chemical.is_produced_by) == 1
                if chemical.is_produced_by:
                    edge = (reaction.identifier, chemical.smiles)
                else:
                    edge = (chemical.smiles, reaction.identifier)
                edges.append(edge)
            network_reactions.append(reaction)
        molecular_smiles = sorted(set(molecular_smiles))
        edges = sorted(set(edges))
        return cls(chemical_reactions=network_reactions, molecular_nodes=molecular_smiles, edges=edges)

    def dummy_quantify(self, by_mass=False):
        """
        dummy case for quantifying a network

        :return:
        """
        for cr in self.chemical_reactions:
            cr.intended_ratios = {c.identifier: 1.0 for c in cr.reactants + cr.reagents}
            cr.reaction_extent = 1.0
            cr.expected_yields = {cr.products[0].identifier: 1.0}
        self.quantify(
            targets={smi: 1.0 for smi in self.target_smiles},
            intended_ratios={cr.identifier: None for cr in self.chemical_reactions},
            expected_yields={cr.identifier: None for cr in self.chemical_reactions},
            by_mass=by_mass
        )

    def quantify(self, targets: dict[str, float], intended_ratios: dict[str, dict[str, float] | None],
                 expected_yields: dict[str, float | None], by_mass=True) -> None:
        """
        quantify the reactions in this network. all reactions have to be uni-product.
        (for multi-product reactions, one could simply modify the reaction SMILES if there is only one desired product.)
        for every molecular SMILES, there is at most one reaction has it as the product.

        :param by_mass: if False quantify by moles
        :param expected_yields: the yield table for uni-product reactions, using reaction identifiers as the domain.
        :param intended_ratios: the intended ratios for uni-product reactions, using reaction identifiers as the domain
        :param targets: a map from molecular SMILES to its expected mass in grams, every one of these molecules can only
        be made from one reaction
        :return:
        """
        if not by_mass:
            targets = {k: calculate_mw(k) * v for k, v in targets.items()}  # convert moles to grams

        # lookup table for product -> reaction
        product_smi_to_reaction = dict()
        for cr in self.chemical_reactions:
            assert cr.is_uniproduct, f"network can only be quantified when all reactions are uni-product, but you have: {cr.reaction_smiles}"
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

    def to_cytoscape_elements(self, fig_size=80) -> list[CytoEdge | CytoNode]:
        """
        convert to cytoscape elements where nodes are SMILES and edges are production/consumption
        """
        nodes = []
        edges = []
        g = self.nx_digraph

        # add molecular SMILES nodes
        for molecular_smiles in self.molecular_nodes:
            node_data = CytoNodeData(id=molecular_smiles, label=molecular_smiles, )
            node_data['url'] = drawing_url(molecular_smiles, size=fig_size)
            chemicals = self.molecular_smiles_table[molecular_smiles]
            node_data['list_of_chemicals'] = [c.model_dump() for c in chemicals]
            node = CytoNode(data=node_data, classes="RN__MolecularSMILES")
            if g.in_degree(molecular_smiles) == 0 and g.out_degree(molecular_smiles) > 0:
                node['classes'] += " RN__MolecularSMILES_starting_material"
            elif g.in_degree(molecular_smiles) > 0 and g.out_degree(molecular_smiles) == 0:
                node['classes'] += " RN__MolecularSMILES_targeting_material"
            else:
                node['classes'] += " RN__MolecularSMILES_intermediate_material"
            nodes.append(node)

        # add chemical reaction nodes
        for cr in self.chemical_reactions:
            node_data = CytoNodeData(id=cr.identifier, label=cr.reaction_smiles, )
            node_data['url'] = drawing_url(cr.reaction_smiles, size=fig_size)
            node_data['chemcail_reaction'] = cr.model_dump()
            node = CytoNode(data=node_data, classes="RN__ChemicalReaction")
            nodes.append(node)

        for u, v in self.edges:
            edge_data = CytoEdgeData(id=str((u, v)), source=u, target=v)
            reaction_identifier, molecular_smiles, relation = (u, v, "production") if is_uuid(u) else (
                v, u, "consumption")
            reaction = self.entity_dictionary[reaction_identifier]
            chemical_in_this_reaction = reaction.smiles2chemical[molecular_smiles]
            edge_data['reaction_to_chemical_relation'] = relation
            edge_data['chemical_in_this_reaction'] = chemical_in_this_reaction.model_dump()
            edge = CytoEdge(data=edge_data, classes="RN__Edge")
            edges.append(edge)
        return nodes + edges

    @property
    def summary(self):
        solids = []
        liquids = []
        for c in self.entity_dictionary.values():
            if isinstance(c, Chemical):
                if c.state_of_matter == StateOfMatter.LIQUID:
                    liquids.append(c)
                elif c.state_of_matter == StateOfMatter.SOLID:
                    solids.append(c)

        return {
            "# reactions": len(self.chemical_reactions),
            "# targets": len([n for n in self.molecular_nodes if self.nx_digraph.out_degree(n) == 0]),
            "# molecular SMILES": len(self.molecular_smiles_table),
            "# solids": len(solids),
            "# liquids": len(liquids),
            "# starting materials": len([n for n in self.molecular_nodes if self.nx_digraph.in_degree(n) == 0]),
        }
