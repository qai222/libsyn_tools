from __future__ import annotations

import copy
import json
from typing import Type, TypeVar

from pydantic import Field
from twa.data_model.base_ontology import DatatypeProperty, ObjectProperty

from libsyn_tools.chem_schema import Chemical
from .base import Individual, SimOntology, BaseClass

T = TypeVar("T", bound=BaseClass)

class Has_ingredient(DatatypeProperty):
    rdfs_isDefinedBy = SimOntology


class PortionOfMaterial(Individual):
    """ a portion of some material that has a chemical composition """

    has_ingredient: Has_ingredient[str] = Field(default_factory=set)
    """
    JSON strings of libsyn_tools.chem_schema.chemical.Chemical objects
    """

    is_present: Is_present[bool] = {False, }

    is_directly_contained_by: Is_directly_contained_by[LabObject] = Field(default_factory=set)

    @property
    def volume(self):
        v = 0
        for c in self.get_ingredients():
            v += c.volume
        return v

    def add_chemical(self, chemical: Chemical):
        """
        directly add chemicals to this pom, usually not used in action effect directly
        """

        self.has_ingredient.add(json.dumps(chemical.model_dump()))

    def get_portion_by_volume(self, volume: float) -> PortionOfMaterial:
        assert 0 < volume <= self.volume + 1e-9
        portion_size = volume / self.volume
        return self.get_portion(portion_size)


    def get_portion(self, portion_size: float) -> PortionOfMaterial:
        """
        create a new individual that has the given portion size
        all of its ingredients are portioned respectively

        :param portion_size:
        :return:
        """
        assert 0 < portion_size <= 1, "cannot expand/zero a portion of material"
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
        """
        mixing with another pom to yield a new pow
        """
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
    """ not transitive """
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


class Has_capacity(DatatypeProperty):
    """
    a relation between A and B where A is a part of B and there is no intermediate layer between A and B,
    e.g. a box in a room in a house, the box is a part of the house but there is an intermediate layer (the room)
    between the box and the house
    """
    # TODO this should be a sub property of `Is_part_of` and it is not transitive
    rdfs_isDefinedBy = SimOntology


class LabObject(Individual):
    is_made_of: Is_made_of[str] = Field(default_factory=set)

    is_present: Is_present[bool] = Field(default={False, })

    is_directly_contained_by: Is_directly_contained_by[LabObject] = Field(default_factory=set)

    # is_part_of: Is_part_of[LabObject] = set()
    is_immediate_part_of: Is_immediate_part_of[LabObject] = Field(default_factory=set)

    has_capacity: Has_capacity[float] = Field(default_factory=set)

    # TODO location?
    # TODO capacity?

    @property
    def capacity(self):
        assert len(self.has_capacity) == 1, f"calling capacity for an ill-specified lab object: {self.has_capacity}"
        return next(self.has_capacity)

    def would_overfill(self, delta: float) -> bool:
        return self.current_volume + delta > self.capacity

    @property
    def current_volume(self) -> float:
        vol = 0
        for pom in PortionOfMaterial.object_lookup.values():
            pom: PortionOfMaterial
            if self in pom.is_directly_contained_by and pom.is_present == {True}:
                vol += pom.volume
        return vol

    @staticmethod
    def get_directly_contained_individuals(
            container: LabObject,
            instance_class: Type[T],
            only_present: bool =True
    ) -> list[T]:
        # TODO pls tell me there is a faster way...
        # TODO can we have a function returns SPARQL results as BaseClass instances?
        directly_contains = []
        for instance in instance_class.object_lookup.values():
            if only_present and instance.is_present != {True}:
                continue
            if container in instance.is_directly_contained_by:
                directly_contains.append(instance)
        return directly_contains

    @staticmethod
    def get_immediate_parts(lab_object: LabObject):
        parts = []
        for instance in LabObject.object_lookup.values():
            if lab_object in instance.is_immediate_part_of:
                parts.append(instance)
        return parts

    @staticmethod
    def get_all_parts(lab_object: LabObject, visited=None):
        if visited is None:
            visited = set()
        # Get the direct parts of the current lab_object
        direct_parts = LabObject.get_immediate_parts(lab_object)
        for part in direct_parts:
            if part not in visited:
                visited.add(part)
                # Recursively collect parts of the current part
                LabObject.get_all_parts(part, visited)
        return visited
