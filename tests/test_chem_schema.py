import networkx as nx
import pytest

from libsyn_tools.chem_schema import *
from libsyn_tools.utils import is_uuid, parse_sparrow_routes, StateOfMatter


@pytest.fixture(scope="session")
def test_reaction_network1():
    reactions = [
        ChemicalReaction.from_reaction_smiles("CC>>c1ccccc1"),
        ChemicalReaction.from_reaction_smiles("CC>COC>c1cScc1"),
        ChemicalReaction.from_reaction_smiles("c1cScc1>>c1cNcc1"),
    ]
    return ReactionNetwork.from_reactions(reactions)


@pytest.fixture(scope="session")
def test_reaction_network2():
    reaction_smiles = parse_sparrow_routes(routes_file="data/routes.json")
    reactions = [ChemicalReaction.from_reaction_smiles(smi) for smi in reaction_smiles]
    return ReactionNetwork.from_reactions(reactions)


@pytest.fixture(scope="session")
def quantified_reaction_1():
    s1 = "CCCCCCCCCCCCCCCCCCCC"
    s2 = "SiSiSiSiSiSiSiSiSiSi"
    l1 = "CC"
    l2 = "CCOCC"
    reactions_smiles = ".".join([s1, s2, l1, l2]) + ">>CCSiCC"
    reaction = ChemicalReaction.from_reaction_smiles(reactions_smiles)
    reaction.quantify(
        product_mass=2.3, product_smiles="CCSiCC",
        intended_ratios={l1: 1, l2: 10.2, s1: 1.4, s2: 1.3},
        expected_yield=0.8
    )
    return reaction


@pytest.mark.order(2)
class TestReactionNetwork:

    def test_molecular_smiles_table(self, test_reaction_network1):
        assert len(test_reaction_network1.molecular_smiles_table) == 5

    def test_nx_digraph(self, test_reaction_network1, test_reaction_network2):
        for test_network in [test_reaction_network1, test_reaction_network2]:
            g = test_network.nx_digraph
            assert len(g.nodes) in [8, 80]
            for u, v in g.edges:
                molecule_smiles = None
                reaction_instance = None
                if is_uuid(u):
                    reaction_instance = test_network.entity_dictionary[u]
                else:
                    molecule_smiles = u
                if is_uuid(v):
                    reaction_instance = test_network.entity_dictionary[v]
                else:
                    molecule_smiles = v
                assert molecule_smiles is not None
                assert reaction_instance is not None
                assert molecule_smiles in reaction_instance.reaction_smiles

    def test_entity_dictionary(self, test_reaction_network1, test_reaction_network2):
        for test_network in [test_reaction_network1, test_reaction_network2]:
            for k, v in test_network.entity_dictionary.items():
                assert is_uuid(k)
                assert isinstance(v, (ChemicalReaction, Chemical))

    def test_quantify(self, test_reaction_network1, test_reaction_network2):
        for test_network in [test_reaction_network1, test_reaction_network2]:
            test_network.dummy_quantify()
            for cr in test_network.chemical_reactions:
                for c in cr.chemicals:
                    assert c.moles == 1


@pytest.mark.order(3)
class TestOperationGraph:

    def test_from_reaction(self, quantified_reaction_1):
        og = OperationGraph.from_reaction(quantified_reaction_1)
        assert len(og.starting_operations) == 2
        assert og.ending_operation
        assert og.lmin
        assert og.lmax


@pytest.mark.order(4)
class TestOperationNetwork:

    def test_from_reaction_network(self, test_reaction_network1, test_reaction_network2):
        test_reaction_network1.dummy_quantify()
        network1 = OperationNetwork.from_reaction_network(test_reaction_network1)
        assert len([*nx.connected_components(network1.nx_digraph.to_undirected())]) == 2
        # TODO: the operation graph has 2 connected components because while CC>>c1ccccc1 and CC>COC>c1cScc1
        #  share the same starting CC, transferring CC does not depend on any operations.
        #  So long there is a SMILES in the reaction network has in_degree == 0 and out_degree > 1
        #  then this would happen.

        test_reaction_network2.dummy_quantify()
        # the following reaction does not have liquid
        for r in test_reaction_network2.chemical_reactions:
            if len(r.reactants + r.reagents) == 1:
                r.reactants[0].state_of_matter = StateOfMatter.LIQUID
        network2 = OperationNetwork.from_reaction_network(test_reaction_network2)
        assert len([*nx.connected_components(network2.nx_digraph.to_undirected())]) > 3
