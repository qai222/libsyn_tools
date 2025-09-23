from __future__ import annotations

"""
SPPT core for the simulator – Substance–Process–Place–Time

Minimal primitives (materialist, Aristotelian-flavored):
- Substance  (material continuants): lab objects (vials, plates, devices, POMs)
- Process    (occurrents/events): transfer, heat, move, interaction
- Place      (sites/environments): heater slot, glovebox, ambient
- Interval   (time spans): [t0, t1]

Primitive relations:
- part_of(Substance, Substance)                [transitive mereology]
- has_participant(Process, Substance)          [bridge things→events; roles travel as data]
- occurs_in(Process, Place)                    [spatial/environmental context]
- precedes(Process, Process)                   [ordering; causal variants can be layered]
- has_interval(Process, Interval)              [begin, end times]

Notes
-----
• We keep domain conveniences like "is_directly_contained_by" elsewhere; for SPPT, we align
  dynamic location to Process+Place via occurs_in, and to convenience in overlay.
• This module defines **OWL-like** classes/props via the TWA/Pydantic base. We don't
  assert OWL imports; instead we echo BFO/RO/Time semantics in docstrings, and keep code simple.
"""

from pydantic import Field
from twa.data_model.base_ontology import BaseClass
from .base import Individual, SimObjectProperty, SimFunctionalDataProperty


# ------------------------
# SPPT Classes
# ------------------------

class Process(Individual):
    """An occurrent/event (e.g., transfer, heat, interaction)."""
    pass


class Place(Individual):
    """A site or environment (e.g., heater slot, glovebox, ambient air)."""
    pass


class TimeInterval(Individual):
    """A time interval; attach begin and end as data properties."""
    pass


# ------------------------
# SPPT Data properties (interval)
# ------------------------

class Has_begin_time(SimFunctionalDataProperty):
    """Begin time (float seconds or ISO timestamp string)."""
    pass


class Has_end_time(SimFunctionalDataProperty):
    """End time (float seconds or ISO timestamp string)."""
    pass


# Attach begin/end to TimeInterval
TimeInterval.has_begin_time: Has_begin_time[float | str] = Field(default_factory=set)  # type: ignore[attr-defined]
TimeInterval.has_end_time: Has_end_time[float | str] = Field(default_factory=set)  # type: ignore[attr-defined]


# ------------------------
# SPPT Object properties (primitives)
# ------------------------

class Part_of(SimObjectProperty):
    """Transitive mereology: x part_of y."""
    # NOTE: we don't mark transitive here programmatically; treated as such in overlay/queries.
    pass


class Has_participant(SimObjectProperty):
    """Event participation: process has_participant substance (roles carried elsewhere)."""
    pass


class Occurs_in(SimObjectProperty):
    """Process occurs_in Place (environment, site)."""
    pass


class Precedes(SimObjectProperty):
    """Process precedes Process (ordering; causal variants can be layered separately)."""
    pass


class Has_interval(SimObjectProperty):
    """Process has_interval TimeInterval."""
    pass


# ------------------------
# Role attachment (lightweight)
# ------------------------

class Has_role(SimFunctionalDataProperty):
    """
    Role label for a (process, substance) pair, stored on the Substance for convenience in this codebase.
    For richer n-ary participation, reify a Participation node; here we keep it minimal and pragmatic.
    """
    pass


# For convenience, let Substance carry a set of (string) roles it played most recently.
# (If you need true n-ary participation with per-event role, reify a Participation class.)
BaseClass.has_role = Field(default_factory=set)  # type: ignore[attr-defined]


# ------------------------
# Compatibility helpers
# ------------------------

# Dynamic content/location convenience – align your existing "is_directly_contained_by" as a domain-level convenience.
# Structural assemblies (well→plate) should use Part_of; dynamic location is better derived from Process+Place.
# We declare it here to make the alignment explicit for future refactors.

class Directly_contained_in(SimObjectProperty):
    """
    Domain convenience for dynamic location (POM in container).
    Conceptually aligns with a "located_in" pattern (not parthood).
    """
    pass

Process.model_rebuild()
Place.model_rebuild()
TimeInterval.model_rebuild()