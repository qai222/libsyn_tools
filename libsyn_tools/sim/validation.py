from __future__ import annotations

from collections.abc import Iterable
from typing import TypeVar

T = TypeVar("T")


def require_singleton(values: Iterable[T] | None, field_name: str, owner_id: str | None = None) -> T:
    if values is None:
        values_list: list[T] = []
    else:
        values_list = list(values)
    owner_label = f" for {owner_id}" if owner_id else ""
    if not values_list:
        raise ValueError(f"{field_name}{owner_label} must have exactly one value; found none")
    if len(values_list) > 1:
        formatted_values = ", ".join(str(value) for value in values_list)
        raise ValueError(
            f"{field_name}{owner_label} must have exactly one value; found {len(values_list)} values: "
            f"[{formatted_values}]"
        )
    return values_list[0]


def require_singleton_or_error(
    values: Iterable[T] | None,
    field_name: str,
    owner_id: str | None = None,
    *,
    context: str | None = None,
) -> T:
    try:
        return require_singleton(values, field_name, owner_id)
    except ValueError as exc:
        owner_label = f"{owner_id}" if owner_id else "unknown"
        context_label = f" ({context})" if context else ""
        raise ValueError(f"{field_name} for {owner_label} is invalid{context_label}: {exc}") from exc
