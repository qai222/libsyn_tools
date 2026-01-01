from __future__ import annotations

"""
Policy configuration for SHACL contract evaluation.

Defaults are audit-only (soft/committed). Supply a PolicyBundle with per-shape
overrides to mark violations as hard/aborted for enforcement.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass(frozen=True)
class PolicyRule:
    """Policy for a single shape (severity + disposition)."""
    severity: str
    disposition: str


@dataclass
class PolicyBundle:
    """
    Policy bundle for SHACL evaluation.

    `defaults` apply to any shape without an explicit override in `per_shape`.
    """
    defaults: PolicyRule = field(default_factory=lambda: PolicyRule(severity="soft", disposition="committed"))
    per_shape: dict[str, PolicyRule] = field(default_factory=dict)

    def rule_for_shape(self, shape_iri: Optional[str]) -> PolicyRule:
        if shape_iri and shape_iri in self.per_shape:
            return self.per_shape[shape_iri]
        return self.defaults
