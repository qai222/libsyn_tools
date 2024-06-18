import pytest

from libsyn_tools.utils.io import parse_sparrow_routes, parse_chemscraper_output


@pytest.mark.order(1)
class TestIO:

    def test_parse_sparrow_routes(self):
        reaction_smiles = parse_sparrow_routes("data/routes.json")
        assert len(reaction_smiles) == 30
        for smi in reaction_smiles:
            assert ">>" in smi

    def test_parse_chemscraper_output(self):
        output = parse_chemscraper_output("data/scraper_output.json")
        for smi, props in output.items():
            state_of_matter, density = props
            assert isinstance(density, float) or density is None
            assert isinstance(state_of_matter, str) or state_of_matter is None
