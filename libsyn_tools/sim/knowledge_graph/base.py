from __future__ import annotations

from typing import Any

from pydantic import Field
from twa.data_model.base_ontology import BaseClass, BaseOntology, DatatypeProperty

from libsyn_tools.utils import str_uuid


class SimOntology(BaseOntology):
    base_url = "https://libsyn-sim/kg/"
    namespace = "libsyn-sim"
    owl_versionInfo = "2"
    rdfs_comment = 'This is an ontology for the chemistry simulator in library synthesis tools.'


class Individual(BaseClass):
    """ a thing in the knowledge graph """

    def model_post_init(self, __context: Any) -> None:
        # NOTE adding this as it seems to be necessary for other actually overwritten methods to
        # work when multi-inheritance is used
        return super().model_post_init(__context)

    rdfs_isDefinedBy = SimOntology
    """ set default ontology """

    instance_iri: str = Field(default_factory=str_uuid, alias='identifier')
    """ instance iri, by default this generated using uuid4 """

    is_present: Is_present[bool] = {False, }

    model_config = {"arbitrary_types_allowed": True, }

    @property
    def identifier(self) -> str:
        return self.instance_iri


class Is_present(DatatypeProperty):
    """ if a lab object is present or has been annihilated """
    rdfs_isDefinedBy = SimOntology
    owl_maxQualifiedCardinality = 1


Individual.model_rebuild()
