from __future__ import annotations

from enum import Enum
from typing import Any

from loguru import logger

from libsyn_tools.utils import FilePath, query_askcos_condition_rec
from .chemical import Chemical, Entity


class ChemicalReactionSpecificationLevel(str, Enum):
    """
    A chemical reaction can be specified at different levels

    - the transformation level
        - stoichiometry
        - chemical identities
    - the quantified level
        - intended molar ratios (different from stoichiometry)
        - reaction time
        - expected yields
        - reaction extent (the actual consumption/production of chemicals),
            for a product its production is expected_yield * reaction_extent * stoichiometry coefficient

    depends on the level of specification of chemical reactions, the reaction network can be quantified or unquantified.

    A chemical reaction specified at the transformation level cannot have quantified chemicals.

    A chemical reaction specified at the quantified level can have the quantities of the chemicals
    inferred from reaction extent and intended molar ratios.
    """
    TRANSFORMATION = "TRANSFORMATION"
    QUANTIFIED = "QUANTIFIED"


class ChemicalReaction(Entity):
    """
    The class for both transformation and quantified chemical reactions.
    """

    reactants: list[Chemical]
    """ reaction_smiles.split('>')[0] """

    reagents: list[Chemical]
    """ reaction_smiles.split('>')[1] """

    products: list[Chemical]
    """ reaction_smiles.split('>')[2] """

    stoichiometry: dict[str, float]
    """ a mapping from `chemical.identifier` to `stoichiometric coefficient`, reagents do not appear here """

    temperature: float | None = None
    """ in C """

    intended_ratios: dict[str, float] | None = None
    """ the intended molar ratios for all reactants + reagents involved in the reaction. similar to 'equivalence' or 
    the keys are `chemical.identifier` -- they can be different from stoichiometry as excess amounts can be used """

    expected_yields: dict[str, float] | None = None
    """ the expected yields for all products """

    reaction_extent: float | None = None
    """ the reaction extent factor, the production/consumption of a chemical is calculated by 
    stoichiometry * reaction_extent * (1 if reactant, -1 * expected_yield if product) 
    reagents will always be fully consumed. this concept is similar to 'batch size'. """

    def model_post_init(self, __context: Any) -> None:
        # warn about multi-product reactions
        if not self.is_uniproduct:
            logger.critical(f"current predictive models only work for uniproduct, but this reaction is not, "
                            f"please be cautious about the routes generated for: {self.products}")

        # check specifications for intended ratios
        if self.intended_ratios is not None:
            for c in self.reactants + self.reagents:
                if c.identifier not in self.intended_ratios:
                    logger.critical(f"no intended ratio specified for: {c}")

        # check yield specifications
        if self.expected_yields is not None:
            for p in self.products:
                if p.identifier not in self.expected_yields:
                    logger.critical("no yield specified for: {p}")

        # set object properties
        for c in self.products:
            assert not c.is_produced_by and not c.is_consumed_by
            c.is_produced_by = [self.identifier]
        for c in self.reactants + self.reagents:
            assert not c.is_produced_by and not c.is_consumed_by
            c.is_consumed_by = [self.identifier]

    @property
    def is_uniproduct(self) -> bool:
        return len(self.products) == 1

    @property
    def specification_level(self) -> ChemicalReactionSpecificationLevel:
        """ the specification level, see the doc string of `ChemicalReactionSpecificationLevel` for more info """
        if self.intended_ratios and self.expected_yields and self.reaction_extent:
            return ChemicalReactionSpecificationLevel.QUANTIFIED
        else:
            return ChemicalReactionSpecificationLevel.TRANSFORMATION

    @property
    def reaction_smiles(self) -> str:
        """ the reaction smiles reconstructed """
        assert all(c.smiles for c in self.reactants + self.reagents + self.products)
        smi_reactants = ".".join([r.smiles for r in self.reactants])
        smi_reagents = ".".join([r.smiles for r in self.reagents])
        smi_products = ".".join([r.smiles for r in self.products])
        return ">".join([smi_reactants, smi_reagents, smi_products])

    @property
    def product_smiles(self) -> str:
        """ molecular smiles of the only product """
        assert self.is_uniproduct
        return self.products[0].smiles

    @property
    def unique_molecular_smiles(self) -> list[str]:
        """ a sorted list of molecular smiles for all chemicals """
        smis = []
        for c in self.reactants + self.reagents + self.products:
            smis.append(c.smiles)
        return sorted(set(smis))

    @classmethod
    def from_reaction_smiles(cls, reaction_smiles: str, query_askcos: bool = False,
                             scraper_output: FilePath = None) -> ChemicalReaction:
        """
        create a ChemicalReaction from reagent free reaction smiles,
        this method allows sending queries to ASKCOS to obtain reagents and quantities

        :param scraper_output: output from ChemScraper for physical properties
        :param reaction_smiles: the reaction smiles, if sending queries to ASKCOS then it should be *reagent-free*
        :param query_askcos: if sending queries to ASKCOS to obtain reagents and quantities
        :return: a ChemicalReaction instance
        """
        assert reaction_smiles.count(">") == 2, f"invalid reaction smiles: {reaction_smiles}"

        stoichiometry = dict()
        identifier2smiles = dict()  # within one reaction this should be bijective
        r_smiles1, r_smiles2, r_smiles3 = reaction_smiles.split(">")

        # reactants
        reactants = []
        for smi in r_smiles1.split('.'):
            chemical = Chemical.from_smiles(smi, scraper_output)
            reactants.append(chemical)
            identifier2smiles[chemical.identifier] = smi
            stoichiometry[chemical.identifier] = 1.0

        # products
        products = []
        for smi in r_smiles3.split('.'):
            chemical = Chemical.from_smiles(smi, scraper_output)
            products.append(chemical)
            identifier2smiles[chemical.identifier] = smi
            stoichiometry[chemical.identifier] = -1.0

        # if query ASKCOS, get reagents and conditions from the query response
        if query_askcos:
            intended_molar_ratios = dict()
            assert ">>" in reaction_smiles, (f"requiring ASKCOS query but "
                                             f"the reaction smiles is not reagent-free: {reaction_smiles}")
            response, query = query_askcos_condition_rec(reaction_smiles, return_query=True)
            result = response['result'][0]
            temperature = result['temperature']
            smi2ratio_reactants = result['reactants']
            smi2ratio_reagents = result['reagents']

            # use askcos reagents
            reagents = []
            for smi, ratio in smi2ratio_reagents.items():
                chemical = Chemical.from_smiles(smi, scraper_output)
                reagents.append(chemical)
                identifier2smiles[chemical.identifier] = smi
                intended_molar_ratios[chemical.identifier] = ratio  # assuming unique

            # update reactant ratios
            smiles2identifier = {v: k for k, v in identifier2smiles.items()}
            for smi, ratio in smi2ratio_reactants.items():
                intended_molar_ratios[smiles2identifier[smi]] = ratio
        # otherwise use what the reaction smiles gives
        else:
            intended_molar_ratios = None
            temperature = None
            reagents = []
            if r_smiles2:
                for smi in r_smiles2.split('.'):
                    chemical = Chemical.from_smiles(smi, scraper_output)
                    reagents.append(chemical)
            else:
                logger.warning(f"no reagents found in the SMILES of: {reaction_smiles}")

        return cls(reactants=reactants, products=products, reagents=reagents, stoichiometry=stoichiometry,
                   intended_ratios=intended_molar_ratios, temperature=temperature, )

    @property
    def smiles2chemical(self) -> dict[str, Chemical]:
        """ lookup table from smiles to `Chemical` instances of this reaction,
        this assumes smiles are unique (this is not a very robust assumption) """
        return {c.smiles: c for c in self.reactants + self.products + self.reagents}

    def quantify(self, product_smiles: str, product_mass: float, intended_ratios: dict[str, float] | None,
                 expected_yield: float | None) -> None:
        """
        quantify a uni-product transformation
        (technically the 'uni-product' here can be relaxed to 'of only one desired product')

        :param intended_ratios: keyed by molecular SMILES
        :param expected_yield: expected yield for this product
        :param product_smiles: the product smiles
        :param product_mass: the product mass in gram
        :return:
        """
        assert self.is_uniproduct, \
            "cannot quantify multi-product reaction"
        assert len(set(self.smiles2chemical.keys())) == len(self.smiles2chemical), \
            "cannot quantify a reaction with non-unique molecular smiles"

        # add specifications
        product = self.smiles2chemical[product_smiles]
        product.mass = product_mass
        logger.warning(product)

        if expected_yield is None:
            assert self.expected_yields is not None
            expected_yield = self.expected_yields[self.products[0].identifier]
        else:
            self.expected_yields = {product.identifier: expected_yield}

        if intended_ratios is None:
            assert self.intended_ratios is not None
        else:
            self.intended_ratios = {self.smiles2chemical[smi].identifier: val for smi, val in intended_ratios.items()}

        # calculate reaction extent.
        # reaction extent = moles / (the stoichiometry coefficient of the product * its expected yield)
        # for uni-product reactions the stoichiometry coefficient is often -1
        self.reaction_extent = abs(product.moles / expected_yield / self.stoichiometry[product.identifier])

        # set quantities
        self.reset_chemical_quantities(chemical_identifier=None)
        product.mass = product_mass  # redo this as the last line resets product mass
        for rid, ratio in self.intended_ratios.items():  # reagents and reactants
            self.chemical_dictionary[rid].quantify_by_moles(ratio * self.reaction_extent)

    @property
    def chemicals(self) -> list[Chemical]:
        """ the list of all chemicals involved in the reaction """
        return self.reactants + self.products + self.reagents

    @property
    def chemical_dictionary(self) -> dict[str, Chemical]:
        """ map from identifiers to chemical instances of this reaction """
        return {c.identifier: c for c in self.chemicals}

    def reset_chemical_quantities(self, chemical_identifier=None):
        if chemical_identifier is None:
            for c in self.chemicals:
                c.mass = None
        else:
            self.chemical_dictionary[chemical_identifier].mass = None
