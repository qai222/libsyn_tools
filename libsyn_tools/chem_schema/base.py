from __future__ import annotations

from pydantic import BaseModel, Field

from libsyn_tools.utils import str_uuid


class Entity(BaseModel):
    """ a thing that deserves a unique identifier and a node in the knowledge graph """

    identifier: str = Field(default_factory=str_uuid)
    """ a unique identifier for this chemical """
