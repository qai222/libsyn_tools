import pytest

from libsyn_tools.chem_schema import *
from libsyn_tools.utils import is_uuid, parse_sparrow_routes


@pytest.mark.order(2)
class TestReactionNetwork:

    @pytest.fixture(scope="session")
    def test_network1(self):
        reactions = [
            ChemicalReaction.from_reaction_smiles("CC>>c1ccccc1"),
            ChemicalReaction.from_reaction_smiles("CC>COC>c1cScc1"),
            ChemicalReaction.from_reaction_smiles("c1cScc1>>c1cNcc1"),
        ]
        return ReactionNetwork.from_reactions(reactions)

    @pytest.fixture(scope="session")
    def test_network2(self):
        reaction_smiles = parse_sparrow_routes(routes_file="data/routes.json")
        reactions = [ChemicalReaction.from_reaction_smiles(smi) for smi in reaction_smiles]
        return ReactionNetwork.from_reactions(reactions)

    def test_molecular_smiles_table(self, test_network1):
        assert len(test_network1.molecular_smiles_table) == 5

    def test_nx_digraph(self, test_network1, test_network2):
        for test_network in [test_network1, test_network2]:
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

    def test_entity_dictionary(self, test_network1, test_network2):
        for test_network in [test_network1, test_network2]:
            for k, v in test_network.entity_dictionary.items():
                assert is_uuid(k)
                assert isinstance(v, (ChemicalReaction, Chemical))

    def test_quantify(self, test_network1, test_network2):
        for test_network in [test_network1]:
            test_network.dummy_quantify()
            for cr in test_network.chemical_reactions:
                for c in cr.chemicals:
                    assert c.moles == 1
