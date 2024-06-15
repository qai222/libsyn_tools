from __future__ import annotations

from loguru import logger
from pydantic import BaseModel, Field

from libsyn_tools.utils import str_uuid, UNKNOWN, get_provenance_model, StateOfMatter, estimate_property_dummy


class ChemicalBase(BaseModel):
    """
    http://purl.allotrope.org/ontologies/material#AFM_0001097
    A chemical substance is a portion of material that is matter of constant composition best characterized
    by the entities (molecules, formula units, atoms) it is composed of. [IUPAC]
    """

    identifier: str = Field(default_factory=str_uuid)
    """ a unique identifier for this chemical """

    smiles: str = UNKNOWN
    """ the SMILES representation for this chemical """

    state_of_matter: StateOfMatter = StateOfMatter.UNKNOWN
    """ the state of matter of this chemical under 1 atm and room temperature """

    molecular_weight: float = 0
    """ (average) molecular weight of this chemical in g/mol, 0 if unknown """

    density: float = 0
    """ density of this chemical in g/mL, 0 if unknown """

    quantity_value: float = 0
    """ the value of the associated quantity, 0 if unknown """

    quantity_unit: str = UNKNOWN
    """ the unit of the associated quantity """


Chemical_ = get_provenance_model(ChemicalBase, "Chemical_")

Chemical_: type[ChemicalBase]


class Chemical(Chemical_):

    @property
    def has_quantity(self) -> bool:
        """
        :return: if the chemical has an associated quantity
        """
        return self.quantity_value > 0 and self.quantity_unit != UNKNOWN

    @classmethod
    def make_up_from_smiles(cls, smiles: str, quantity_unit: str = UNKNOWN, qunatity_value: float = 0):
        """
        make up a chemical based on a SMILES

        :param smiles: the molecular SMILES
        :return: a Chemical
        """
        logger.warning(f"making up a chemical based on: {smiles}")
        properties = estimate_property_dummy(smiles)
        return cls(
            smiles=smiles,
            state_of_matter=properties['state_of_matter'],
            molecular_weight=properties['mw'],
            density=properties['density'],
            quantity_unit=quantity_unit,
            qunatity_value=qunatity_value,
            provenance__state_of_matter='dummy estimator based on mw',
            provenance__molecular_weight='periodic table',
            provenance__density='dummy estimator based on mw',

        )


Chemical: type[ChemicalBase]
# a hack to add provenance, still couldn't find a way to type hint without explicitly type down the provenance fields
