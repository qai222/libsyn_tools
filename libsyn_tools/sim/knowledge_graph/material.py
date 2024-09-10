from __future__ import annotations

from typing import Optional

from twa.data_model.base_ontology import ObjectProperty

from libsyn_tools import chem_schema
from .base import Individual, SimOntology
from .utils import pydantic_to_ontology_class

Chemical = pydantic_to_ontology_class(pydantic_class=chem_schema.Chemical, base_ontology_class=Individual)


class Has_ingreadient(ObjectProperty):
    rdfs_isDefinedBy = SimOntology


class PortionOfMaterial(Individual):
    has_ingredient: Optional[Has_ingreadient[Chemical]]
