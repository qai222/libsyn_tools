from __future__ import annotations

from typing import Any
from uuid import uuid4

from pydantic import Field
from rdflib import URIRef
from twa.data_model.base_ontology import BaseClass, BaseOntology, DatatypeProperty, ObjectProperty


class SimOntology(BaseOntology):
    base_url = "https://libsyn-sim/kg/"
    namespace = ""
    owl_versionInfo = "2"
    rdfs_comment = 'This is an ontology for the chemistry simulator in library synthesis tools.'


class SimDataProperty(DatatypeProperty):
    rdfs_isDefinedBy = SimOntology


class SimObjectProperty(ObjectProperty):
    rdfs_isDefinedBy = SimOntology


class SimFunctionalDataProperty(SimDataProperty):
    owl_maxQualifiedCardinality = 1


class SimFunctionalObjectProperty(SimObjectProperty):
    owl_maxQualifiedCardinality = 1


class Is_present(SimFunctionalDataProperty):
    """ if a lab object is present or has been annihilated """
    pass


class Individual(BaseClass):
    """ a thing in the knowledge graph """

    def model_post_init(self, __context: Any) -> None:
        # NOTE adding this as it seems to be necessary for other actually overwritten methods to
        # work when multi-inheritance is used
        return super().model_post_init(__context)

    rdfs_isDefinedBy = SimOntology
    """ set default ontology """

    instance_iri: str = Field(default=None, alias='identifier')
    """ instance iri, by default this generated using uuid4 """

    is_present: Is_present[bool] = Field(default={False, })

    model_config = {"arbitrary_types_allowed": True, }

    @property
    def identifier(self) -> str:
        return self.instance_iri

    @classmethod
    def __init_subclass__(cls, **kwargs):
        """
        Every time a concrete subclass is defined, capture *that* class name
        and install a tailored default_factory for `instance_iri`.
        """
        super().__init_subclass__(**kwargs)

        field_info = cls.model_fields['instance_iri']
        field_info.default_factory = lambda c=cls: f"{c.__name__}_{uuid4()}"

        # mark the field non-nullable
        field_info.annotation = str  # removes Optional[...] in schema


Individual.model_rebuild()


def canonical_iri(instance_iri: str) -> URIRef:
    """
    Convert an instance IRI or identifier into the canonical KG URIRef.
    """
    base_url = SimOntology.base_url
    if instance_iri.startswith(base_url):
        return URIRef(instance_iri)
    return URIRef(f"{base_url}{instance_iri}")


def identifier_from_iri(instance_iri: str) -> str:
    """
    Normalize a canonical KG IRI back to its identifier form.
    """
    base_url = SimOntology.base_url
    if instance_iri.startswith(base_url):
        return instance_iri[len(base_url):]
    return instance_iri
