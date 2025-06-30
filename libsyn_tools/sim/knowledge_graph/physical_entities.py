from __future__ import annotations

import json
from typing import Type, TypeVar

from pydantic import Field

from libsyn_tools.chem_schema import Chemical
from .base import Individual, BaseClass, SimFunctionalDataProperty, SimDataProperty, SimObjectProperty

T = TypeVar("T", bound=BaseClass)


class Has_ingredient(SimDataProperty):
    """ data property to JSON string of libsyn_tools.chem_schema.chemical.Chemical """
    pass


class Is_contained_by(SimObjectProperty):
    """ transitive """
    # TODO transitive
    pass


class Is_directly_contained_by(SimObjectProperty):
    """ not transitive """
    # TODO subproperty of Is_contained_by
    pass


class PortionOfMaterial(Individual):
    """ a portion of some material that has a chemical composition """

    has_ingredient: Has_ingredient[str] = Field(default_factory=set)
    """
    JSON strings of libsyn_tools.chem_schema.chemical.Chemical objects
    """

    is_directly_contained_by: Is_directly_contained_by[LabObject] = Field(default_factory=set)

    @property
    def volume(self):
        v = 0
        for c in self.get_ingredients():
            if c.volume is None:
                raise ValueError(f"Chemical {c!r} lacks a `volume` value")
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
        assert 0 <= portion_size <= 1, f"cannot expand a portion of material: {portion_size}"
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


class Is_made_of(SimDataProperty):
    pass


class Is_part_of(SimObjectProperty):
    """ can be proper or improper, transitive """
    # TODO transitive
    pass


class Is_immediate_part_of(SimObjectProperty):
    """ not transitive """
    # TODO this should be a sub property of `Is_part_of` and it is not transitive
    pass


class Has_pool_type(SimFunctionalDataProperty):
    # TODO we only allow one pool type for now
    pass


class LabObject(Individual):
    """
    Pure-data representation of anything that can appear in the lab KG.
    """

    is_made_of: Is_made_of[str] = Field(default_factory=set)

    is_directly_contained_by: Is_directly_contained_by[LabObject] = Field(default_factory=set)

    is_immediate_part_of: Is_immediate_part_of[LabObject] = Field(default_factory=set)
    # TODO location?

    has_pool_type: Has_pool_type[str] = Field(default_factory=set)
    """ 
    Objects with the same pool type will be grouped in a simpy filter store so it can be picked up by a selector 
    see `libsyn_tools.sim.selector.FilterStoreRegistry` for more info
    """

    @staticmethod
    def get_directly_contained_individuals(
            container: LabObject,
            instance_class: Type[T] = None,
            only_present: bool = True
    ) -> list[T]:
        # TODO pls tell me there is a faster way...
        # TODO can we have a function returns SPARQL results as BaseClass instances?
        directly_contains = []
        if instance_class is None:
            raise NotImplementedError(
                "instance_class = None is not working right now: targe_class cannot set to be BaseClass, "
                "this should be fixed in the future"
            )
        else:
            target_class = instance_class
        for instance in target_class.object_lookup.values():
            if only_present and instance.is_present != {True}:
                continue
            if container in instance.is_directly_contained_by:
                directly_contains.append(instance)
        return directly_contains

    def directly_contained_individuals(self, instance_class: Type[T], only_present: bool = True):
        return LabObject.get_directly_contained_individuals(self, instance_class, only_present)

    def all_contained_individuals(self, instance_class: Type[T], only_present: bool = True):
        raise NotImplementedError

    def get_all_contained_individuals(self):
        # TODO get all contained individuals
        raise NotImplementedError

    @property
    def directly_contained_pom_volume(self) -> float:
        vol = 0
        for pom in PortionOfMaterial.object_lookup.values():
            pom: PortionOfMaterial
            if self in pom.is_directly_contained_by and pom.is_present == {True}:
                vol += pom.volume
        return vol

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


class Has_capacity(SimFunctionalDataProperty):
    pass


class MaterialContainer(LabObject):
    has_capacity: Has_capacity[float] = Field(default_factory=set)
    """
    Note this **does not** map to simpy resource capacity: simpy resource capacity is always 1.
    We don't use simpy container even it supports continuous material quantity: it only deals with constant composition.
    """

    @property
    def capacity(self):
        return next(self.has_capacity)


LabObject.model_rebuild()
PortionOfMaterial.model_rebuild()
MaterialContainer.model_rebuild()
