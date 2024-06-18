from __future__ import annotations

from loguru import logger

from libsyn_tools.utils import get_provenance_model, StateOfMatter, estimate_property_dummy, FilePath, \
    parse_chemscraper_output
from .base import Entity


class ChemicalBase(Entity):
    """
    http://purl.allotrope.org/ontologies/material#AFM_0001097
    A chemical substance is a portion of material that is matter of constant composition best characterized
    by the entities (molecules, formula units, atoms) it is composed of. [IUPAC]
    """

    smiles: str | None = None
    """ the SMILES representation for this chemical """

    state_of_matter: StateOfMatter | None = None
    """ the state of matter of this chemical under 1 atm and room temperature """

    molecular_weight: float | None = None
    """ (average) molecular weight of this chemical in g/mol """

    density: float | None = None
    """ density of this chemical in g/mL """

    mass: float | None = None
    """ mass of this chemical in g """

    is_produced_by: list[str] = []
    """ object property to reaction """

    is_consumed_by: list[str] = []
    """ object property to reaction """

    def __repr__(self):
        return (f"{self.__class__.__name__}: ({self.identifier})\n"
                f"{self.smiles} @ {self.mass} (g)")

    @property
    def volume(self) -> float:
        return self.mass / self.density

    @property
    def moles(self):
        return self.mass / self.molecular_weight

    def quantify_by_moles(self, moles: float) -> None:
        self.mass = self.molecular_weight * moles

    def quantify_by_volume(self, volume: float) -> None:
        self.mass = self.density * volume


Chemical_ = get_provenance_model(ChemicalBase, "Chemical_")

Chemical_: type[ChemicalBase]


class Chemical(Chemical_):

    @property
    def has_quantity(self) -> bool:
        """
        :return: if the chemical has an associated quantity
        """
        return self.mass is not None

    @classmethod
    def make_up_from_smiles(cls, smiles: str, quantity_unit: str = None, quantity_value: float = None):
        """
        make up a chemical based on a SMILES

        :param smiles: the molecular SMILES
        :param quantity_unit: the unit of the associated quantity
        :param quantity_value: the value of the associated quantity
        :return: a Chemical instance
        """
        logger.info(f"making up a chemical based on: {smiles}")
        properties = estimate_property_dummy(smiles)
        return cls(
            smiles=smiles,
            state_of_matter=properties['state_of_matter'],
            molecular_weight=properties['mw'],
            density=properties['density'],
            quantity_unit=quantity_unit,
            qunatity_value=quantity_value,
            provenance__state_of_matter='dummy estimator based on mw',
            provenance__molecular_weight='periodic table',
            provenance__density='dummy estimator based on mw',
        )

    @classmethod
    def from_smiles(cls, smiles: str, scraper_output: FilePath = None, quantity_unit: str = None,
                    quantity_value: float = None):
        """
        create an unquantified chemical instance from a SMILES string

        :param quantity_unit: the unit of the associated quantity
        :param quantity_value: the value of the associated quantity
        :param smiles: molecular SMILES
        :param scraper_output: the output from `ChemScraper`
        :return:
        """
        if scraper_output is None:
            return Chemical.make_up_from_smiles(smiles=smiles, quantity_unit=quantity_unit,
                                                quantity_value=quantity_value)

        output_data = parse_chemscraper_output(scraper_output)
        if smiles not in output_data:
            return Chemical.make_up_from_smiles(smiles=smiles, quantity_unit=quantity_unit,
                                                quantity_value=quantity_value)
        else:
            state_of_matter, density = output_data[smiles]
            made_up_properties = estimate_property_dummy(smiles)
            form_p = f"from scraper output: {scraper_output}"
            density_p = f"from scraper output: {scraper_output}"
            mw_p = 'periodic table'
            if state_of_matter is None:
                state_of_matter = made_up_properties['state_of_matter']
                form_p = 'dummy estimator based on mw'
            if density is None:
                density = made_up_properties['density']
                density_p = 'dummy estimator based on mw'
        return cls(
            smiles=smiles,
            state_of_matter=state_of_matter,
            molecular_weight=made_up_properties['mw'],
            density=density,
            quantity_unit=quantity_unit,
            qunatity_value=quantity_value,
            provenance__state_of_matter=form_p,
            provenance__molecular_weight=mw_p,
            provenance__density=density_p,
        )


Chemical: type[ChemicalBase]
# a hack to add provenance, still couldn't find a way to type hint without explicitly typing down the provenance fields
