from __future__ import annotations

from pydantic import Field
from twa.data_model.base_ontology import BaseClass
from twa.data_model.base_ontology import BaseOntology


class SimOntology(BaseOntology):
    base_url = "https://libsyn-sim/kg/"
    namespace = "libsyn-sim"
    owl_versionInfo = "0.0.1"
    rdfs_comment = 'This is an ontology for the platform simulator in library synthesis tools.'


class Individual(BaseClass):
    """ a thing in the knowledge graph """

    rdfs_isDefinedBy = SimOntology
    """ set default ontology """

    instance_iri: str = Field(default='', alias='identifier')
    """ instance iri """

    @property
    def identifier(self) -> str:
        return self.instance_iri
