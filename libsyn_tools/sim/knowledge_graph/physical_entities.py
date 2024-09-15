from __future__ import annotations

import json
from typing import Any, Type

from twa.data_model.base_ontology import DatatypeProperty
from twa.data_model.base_ontology import ObjectProperty

from .base import Individual, SimOntology, BaseClass
from ...chem_schema import Chemical


class Has_ingredient(DatatypeProperty):
    rdfs_isDefinedBy = SimOntology


class PortionOfMaterial(Individual):
    """ a portion of some material that has a chemical composition """

    has_ingredient: Has_ingredient[str] = set()

    is_present: Is_present[bool] = {False, }

    is_directly_contained_by: Is_directly_contained_by[LabObject] = set()

    def add_chemical(self, chemical: Chemical):
        self.has_ingredient.add(json.dumps(chemical.model_dump()))

    def get_portion(self, portion_size: float) -> PortionOfMaterial:
        """
        create a new individual that has the given portion size
        all of its ingredients are portioned respectively

        :param portion_size:
        :return:
        """
        assert portion_size <= 1, "cannot expand a portion of material"
        assert len(self.has_ingredient), "the portion of material has no ingredient"

        new_pom = PortionOfMaterial()
        chemicals = self.get_ingredients()
        for chemical in chemicals:
            new_chemical = chemical.split([portion_size, 1 - portion_size])[0]
            new_pom.has_ingredient.add(json.dumps(new_chemical.model_dump()))
        return new_pom

    def get_ingredients(self) -> list[Chemical]:
        return [Chemical(**json.loads(chemical_json_str)) for chemical_json_str in self.has_ingredient]

    def mix_with(self, other: PortionOfMaterial) -> PortionOfMaterial:
        new_pom = PortionOfMaterial()
        chemicals = self.get_ingredients() + other.get_ingredients()
        for chemical in chemicals:
            new_pom.has_ingredient.add(json.dumps(chemical.model_dump()))
        return new_pom


class Is_made_of(DatatypeProperty):
    rdfs_isDefinedBy = SimOntology


class Is_present(DatatypeProperty):
    """ if a lab object is present or has been annihilated """
    rdfs_isDefinedBy = SimOntology
    owl_maxQualifiedCardinality = 1


class Is_directly_contained_by(ObjectProperty):
    rdfs_isDefinedBy = SimOntology


class Is_part_of(ObjectProperty):
    """
    can be proper or improper, transitive
    """
    # TODO transitive
    rdfs_isDefinedBy = SimOntology


class Is_immediate_part_of(ObjectProperty):
    """
    a relation between A and B where A is a part of B and there is no intermediate layer between A and B,
    e.g. a box in a room in a house, the box is a part of the house but there is an intermediate layer (the room)
    between the box and the house
    """
    # TODO this should be a sub property of `Is_part_of` and it is not transitive
    rdfs_isDefinedBy = SimOntology


class LabObject(Individual):

    def model_post_init(self, __context: Any) -> None:
        # TODO: do we want to put this in `Individual`?
        # NOTE adding this as it seems to be necessary for other actually overwritten methods to
        # work when multi-inheritance is used
        # i.e. JuniorLabObject and JuniorInstruction
        return super().model_post_init(__context)

    is_made_of: Is_made_of[str] = set()

    is_present: Is_present[bool] = {False, }

    is_directly_contained_by: Is_directly_contained_by[LabObject] = set()

    is_part_of: Is_part_of[LabObject] = set()

    # TODO location?

    @staticmethod
    def get_directly_contained_individuals(container: LabObject, instance_class: Type[BaseClass]):
        # TODO pls tell me there is a faster way...
        # TODO can we have a function returns SPARQL results as BaseClass instances?
        directly_contains = []
        for instance in instance_class.object_lookup.values():
            if container in instance.is_directly_contained_by:
                directly_contains.append(instance)
        return directly_contains

    @staticmethod
    def get_parts(lab_object: LabObject):
        parts = []
        for instance in LabObject.object_lookup.values():
            if lab_object in instance.is_part_of:
                parts.append(instance)
        return parts
