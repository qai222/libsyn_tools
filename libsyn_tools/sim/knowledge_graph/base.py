from __future__ import annotations

from typing import Any

from pydantic import Field
from twa.data_model.base_ontology import BaseClass, BaseOntology

from libsyn_tools.utils import str_uuid


class SimOntology(BaseOntology):
    base_url = "https://libsyn-sim/kg/"
    namespace = "libsyn-sim"
    owl_versionInfo = "0.0.1"
    rdfs_comment = 'This is an ontology for the chemistry simulator in library synthesis tools.'


class Individual(BaseClass):
    """ a thing in the knowledge graph """

    def model_post_init(self, __context: Any) -> None:
        # TODO: do we want to put this in `Individual`?
        # NOTE adding this as it seems to be necessary for other actually overwritten methods to
        # work when multi-inheritance is used
        # i.e. JuniorLabObject and JuniorInstruction
        return super().model_post_init(__context)

    rdfs_isDefinedBy = SimOntology
    """ set default ontology """

    instance_iri: str = Field(default_factory=str_uuid, alias='identifier')
    """ instance iri, by default this generated using uuid4 """

    model_config = {"arbitrary_types_allowed": True, }

    @property
    def identifier(self) -> str:
        return self.instance_iri
