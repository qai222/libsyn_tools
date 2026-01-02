from __future__ import annotations

from typing import Any, List

from .unitary_edit import (
    AddDataProperty,
    AddObjectProperty,
    Annihilate,
    ChangeDataProperty,
    Create,
    RemoveDataProperty,
    RemoveObjectProperty,
    UnitaryEdit,
)


class EffectsBuilder:
    """
    Fluent helper for composing ordered lists of UnitaryEdit instances.
    """

    def __init__(self) -> None:
        self._edits: List[UnitaryEdit] = []

    def create(self, iri: str) -> "EffectsBuilder":
        self._edits.append(Create(instance_1_iri=iri))
        return self

    def annihilate(self, iri: str) -> "EffectsBuilder":
        self._edits.append(Annihilate(instance_1_iri=iri))
        return self

    def set_data(self, subj: str, pred_iri: str, value: Any) -> "EffectsBuilder":
        self._edits.append(
            ChangeDataProperty(
                instance_1_iri=subj,
                property_iri=pred_iri,
                data_value=value,
            )
        )
        return self

    def add_data(self, subj: str, pred_iri: str, value: Any) -> "EffectsBuilder":
        self._edits.append(
            AddDataProperty(
                instance_1_iri=subj,
                property_iri=pred_iri,
                data_value=value,
            )
        )
        return self

    def remove_data(self, subj: str, pred_iri: str, value: Any) -> "EffectsBuilder":
        self._edits.append(
            RemoveDataProperty(
                instance_1_iri=subj,
                property_iri=pred_iri,
                data_value=value,
            )
        )
        return self

    def link(self, subj: str, pred_iri: str, obj: str) -> "EffectsBuilder":
        self._edits.append(
            AddObjectProperty(
                instance_1_iri=subj,
                instance_2_iri=obj,
                property_iri=pred_iri,
            )
        )
        return self

    def unlink(self, subj: str, pred_iri: str, obj: str) -> "EffectsBuilder":
        self._edits.append(
            RemoveObjectProperty(
                instance_1_iri=subj,
                instance_2_iri=obj,
                property_iri=pred_iri,
            )
        )
        return self

    def build(self) -> List[UnitaryEdit]:
        return list(self._edits)
