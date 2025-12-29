from __future__ import annotations

import json
from typing import Iterable, Type, TypeVar

from pydantic import Field
from twa.data_model.base_ontology import BaseClass

from libsyn_tools.chem_schema import Chemical
from .base import Individual, SimDataProperty, SimObjectProperty, SimFunctionalDataProperty


# ======================================================================
# Core SPPT classes (upper ontology)
# ======================================================================

class Substance(Individual):
    """
    Material continuant: vials, plates, devices, portions of material, etc.
    """
    pass


class Process(Individual):
    """
    Occurrent/event: transfer, heat, move, interaction, placement, etc.
    The *participants* and *time interval* of a Process live on relations.
    """
    pass


class Place(Individual):
    """
    Site/environment: heater slot, glovebox, ambient air, solvent phase.
    """
    pass


class TimeInterval(Individual):
    """
    Time span. Attach begin/end via functional data properties.
    """
    pass


# ======================================================================
# Core SPPT relations (minimal primitives)
# ======================================================================

class Part_of(SimObjectProperty):
    """
    Transitive mereology: x part_of y.
    • Intended semantics: transitive closure over assemblies (e.g., well → plate).
    • Do not use for dynamic contents; use Process(+Place) or the convenience
      dynamic location relation(s) declared below.
    """
    pass


class Has_participant(SimObjectProperty):
    """
    Event participation: process has_participant substance.
    • Roles (instrument/source/destination/patient/...) travel as data attributes
      on the participation in overlays or via auxiliary properties.
    """
    pass


class Occurs_in(SimObjectProperty):
    """
    Process occurs_in Place (environment/slot/zone).
    """
    pass


class Precedes(SimObjectProperty):
    """
    Process precedes Process (temporal ordering). Causal variants can be layered
    later if desired; the default is ordering only.
    """
    pass


class Has_interval(SimObjectProperty):
    """
    Process has_interval TimeInterval.
    """
    pass


class Has_begin_time(SimFunctionalDataProperty):
    """
    Begin time (float seconds or ISO timestamp string).
    """
    pass


class Has_end_time(SimFunctionalDataProperty):
    """
    End time (float seconds or ISO timestamp string).
    """
    pass


# Attach interval endpoints on TimeInterval individuals
TimeInterval.has_begin_time: Has_begin_time[float | str] = Field(default_factory=set)  # type: ignore[attr-defined]
TimeInterval.has_end_time: Has_end_time[float | str] = Field(default_factory=set)  # type: ignore[attr-defined]


# ======================================================================
# Domain conveniences (subclasses & convenience relations)
# ======================================================================

# --- Data properties ---

class Has_pool_type(SimFunctionalDataProperty):
    """
    Functional data property.
    • Used to group Substances (e.g., devices/containers) into selection pools.
    """
    pass


class Has_capacity(SimFunctionalDataProperty):
    """
    Functional data property (float).
    • Domain: MaterialContainer
    """
    pass


class Has_ingredient(SimDataProperty):
    """
    Data property carrying JSON-serialized Chemical payload(s).
    • Domain: PortionOfMaterial
    """
    pass


class Has_interrupt_events(SimDataProperty):
    """
    Data property of strings, used to record interrupt reasons on participants.
    """
    pass


class Is_made_of(SimDataProperty):
    """
    Textual tags for material-of-construction, coatings, etc. (domain convenience).
    """
    pass


# --- Object properties (domain conveniences) ---

class Is_immediate_part_of(SimObjectProperty):
    """
    Non-transitive structural link (domain convenience).
    • Conceptual alignment: subPropertyOf Part_of (structural mereology).
    • Use for direct assembly (e.g., well is_immediate_part_of plate).
    """
    pass


class Is_directly_contained_by(SimObjectProperty):
    """
    Non-transitive dynamic location link (domain convenience).
    • Conceptual alignment: akin to *located_in* (not parthood) and *not* transitive.
    • Used by operations for dynamic contents (POM in container).
    • Temporal semantics (exposure, movement) are represented via SPPT events + overlays.
    """
    pass


T = TypeVar("T", bound=BaseClass)
T_co = TypeVar("T_co", bound=BaseClass, covariant=True)


# ======================================================================
# Domain classes (built on SPPT: all are Substances)
# ======================================================================

class PortionOfMaterial(Substance):
    """
    A portion of some material that has a chemical composition.

    Notes
    -----
    • Composition travels via Has_ingredient (JSON-serialized Chemical).
    • Dynamic location: Is_directly_contained_by (domain convenience).
    """
    has_ingredient: Has_ingredient[str] = Field(default_factory=set)
    is_directly_contained_by: Is_directly_contained_by["LabObject"] = Field(default_factory=set)

    @property
    def volume(self) -> float:
        v = 0.0
        for chem in self.get_ingredients():
            if chem.volume is None:
                raise ValueError(f"Chemical {chem!r} lacks a `volume` value")
            v += float(chem.volume)
        return v

    def add_chemical(self, chemical: Chemical) -> None:
        self.has_ingredient.add(json.dumps(chemical.model_dump()))

    def get_portion_by_volume(self, volume: float) -> "PortionOfMaterial":
        assert 0 < volume <= self.volume + 1e-9
        portion_size = volume / self.volume
        return self.get_portion(portion_size)

    def get_portion(self, portion_size: float) -> "PortionOfMaterial":
        assert 0 <= portion_size <= 1, f"cannot expand a portion of material: {portion_size}"
        assert len(self.has_ingredient), "the portion of material has no ingredient"
        new_pom = PortionOfMaterial()
        for chemical in self.get_ingredients():
            new_chemical = chemical.split([portion_size, 1 - portion_size])[0]
            new_pom.has_ingredient.add(json.dumps(new_chemical.model_dump()))
        return new_pom

    def get_ingredients(self) -> list[Chemical]:
        return [Chemical(**json.loads(s)) for s in self.has_ingredient]

    def mix_with(self, other: "PortionOfMaterial") -> "PortionOfMaterial":
        """
        Mixing with another POM yields a new POM whose ingredients are
        the union (concatenation) of the two sets of chemicals.
        Quantities are preserved because we carry the serialized Chemicals.
        """
        new_pom = PortionOfMaterial()
        chemicals = self.get_ingredients() + other.get_ingredients()
        for chemical in chemicals:
            new_pom.has_ingredient.add(json.dumps(chemical.model_dump()))
        return new_pom


class LabObject(Substance):
    """
    Pure-data representation of anything that can appear in the lab KG (devices, containers, tools).
    """
    is_made_of: Is_made_of[str] = Field(default_factory=set)
    is_directly_contained_by: Is_directly_contained_by["LabObject"] = Field(default_factory=set)
    is_immediate_part_of: Is_immediate_part_of["LabObject"] = Field(default_factory=set)
    has_pool_type: Has_pool_type[str] = Field(default_factory=set)
    has_interrupt_events: Has_interrupt_events[str] = Field(default_factory=set)

    @staticmethod
    def get_directly_contained_individuals(
            container: "LabObject",
            instance_class: Type[T_co] | None = None,
            only_present: bool = True,
    ) -> list[T_co]:
        """
        Utility to collect all instances of a class *directly* contained in `container`
        via the dynamic location convenience relation.
        """
        if instance_class is None:
            raise NotImplementedError(
                "instance_class=None is not supported here; supply a concrete class."
            )
        target_class = instance_class
        out: list[T_co] = []
        for inst in target_class.object_lookup.values():
            if only_present and inst.is_present != {True}:
                continue
            if container in inst.is_directly_contained_by:
                out.append(inst)
        return out

    def directly_contained_individuals(self, instance_class: Type[T], only_present: bool = True) -> list[T]:
        return LabObject.get_directly_contained_individuals(self, instance_class, only_present)

    @property
    def directly_contained_pom_volume(self) -> float:
        vol = 0.0
        for pom in PortionOfMaterial.object_lookup.values():
            if self in pom.is_directly_contained_by and pom.is_present == {True}:
                vol += pom.volume
        return vol

    # assembly helpers (structural)
    @staticmethod
    def get_immediate_parts(lab_object: "LabObject") -> list["LabObject"]:
        return [inst for inst in LabObject.all_instances() if lab_object in inst.is_immediate_part_of]

    @staticmethod
    def get_all_parts(lab_object: "LabObject", visited=None):
        if visited is None:
            visited = set()
        direct_parts = LabObject.get_immediate_parts(lab_object)
        for part in direct_parts:
            if part not in visited:
                visited.add(part)
                LabObject.get_all_parts(part, visited)
        return visited

    @classmethod
    def all_instances(cls) -> Iterable["LabObject"]:
        """
        Yield every LabObject (of any subclass) that currently exists in memory — once per instance.
        """
        stack = [cls]
        seen_ids: set[str] = set()
        while stack:
            k = stack.pop()
            for obj in k.object_lookup.values():
                if obj.instance_iri not in seen_ids:
                    seen_ids.add(obj.instance_iri)
                    yield obj
            stack.extend(k.__subclasses__())


class MaterialContainer(LabObject):
    """
    Container/holder that can carry capacity (float).
    """
    has_capacity: Has_capacity[float] = Field(default_factory=set)

    @property
    def capacity(self) -> float:
        cap = next(iter(self.has_capacity), None)
        if cap is None:
            raise AttributeError(f"{self.identifier} has no `has_capacity` set")
        return float(cap)


# Ensure models are rebuilt (TWA)
Substance.model_rebuild()
Process.model_rebuild()
Place.model_rebuild()
TimeInterval.model_rebuild()

PortionOfMaterial.model_rebuild()
LabObject.model_rebuild()
MaterialContainer.model_rebuild()
