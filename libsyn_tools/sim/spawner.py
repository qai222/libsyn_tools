"""
Spawners for endogenous event generation.
A spawner attaches itself to one `Simulation`,
runs as a SimPy coroutine and creates new `Operation`s when its own trigger fires.
"""

from __future__ import annotations

import weakref
from abc import ABC, abstractmethod
from typing import Optional, Callable, Any

import simpy
from pydantic import BaseModel, Field, PrivateAttr
from rdflib.namespace import SH, RDF

from .simulation import Simulation, Operation, OperationProcess


class Spawner(BaseModel, ABC):
    _sim_ref: Optional[weakref.ReferenceType] = PrivateAttr(default=None)  # avoid leaks

    def attach(self, sim: Simulation) -> None:
        """
        Register the spawner with a `Simulation` and start the coroutine.
        Each instance can be attached **once**.
        """
        if self._sim_ref is not None:
            raise RuntimeError("Spawner already attached")
        self._sim_ref = weakref.ref(sim)
        sim.env.process(self._run(sim))  # fire-and-forget

    @property
    def sim(self) -> Simulation:
        s = None if self._sim_ref is None else self._sim_ref()
        if s is None:
            raise RuntimeError("Spawner has not been attached yet")
        return s

    @abstractmethod
    def _run(self, sim: Simulation) -> "simpy.events.Event":
        """Implement the trigger logic as a SimPy process."""


class TimerSpawner(Spawner):
    op_factory: Callable[[Any], Operation]
    interval: float = Field(..., description="delta t in sim time")
    start_offset: float = 0.0
    alive: bool = True  # runtime flag

    def cancel(self) -> None:
        """Stop the timer after the current cycle."""
        self.alive = False

    def _sample_dt(self) -> float:
        # TODO add randomness
        rv = self.interval
        return float(rv)

    def _run(self, sim: Simulation):
        env, rng = sim.env, sim.rng
        if self.start_offset > 0:
            yield env.timeout(self.start_offset)

        while self.alive:
            op = self.op_factory(sim)
            sim.spawn_operation(op)
            yield env.timeout(self._sample_dt())


class KGInspectorSpawner(Spawner):
    """
    Periodically (or on every operation) validate the KG against a
    user-supplied SHACL *sub-set* and spawn maintenance / corrective
    Operations when any of those shapes become non-conformant.

    Parameters
    ----------
    shape_dispatch : dict[str, Callable[[str], Operation]]
        Mapping **sourceShape IRI ➜ factory**.  The factory receives the
        *focusNode IRI* string and must return a fully constructed
        Operation object ready for `spawn_operation()`.
    inspect_interval : float
        • > 0   → run every `inspect_interval` seconds of **sim-time**.
        • == 0  → run right after every `Operation_END`.
    """

    shape_dispatch: dict[str, Callable[[str], Operation]]
    inspect_interval: float = Field(0.0, ge=0.0)

    def _spawn_for_violations(self, sim: Simulation, report) -> None:
        """Iterate ValidationResults and spawn Operations as needed."""
        g = report
        for vr in g.subjects(RDF.type, SH.ValidationResult):
            shape_iri = str(g.value(vr, SH.sourceShape))
            focus_iri = str(g.value(vr, SH.focusNode))
            factory = self.shape_dispatch.get(shape_iri)
            if factory is None:
                continue  # shape not in our interest list
            op = factory(focus_iri)
            sim.spawn_operation(op)

    def _run_every_dt(self, sim: Simulation):
        env = sim.env
        dt = self.inspect_interval
        while True:
            yield env.timeout(dt)
            conforms, report, _ = sim.effect_engine.validate_now()
            if not conforms:
                self._spawn_for_violations(sim, report)

    def _run_on_every_operation(self, sim: Simulation):
        """
        Monkey-patch OperationProcess.add_event_log so we get called
        exactly once per OPERATION_END.
        """
        from types import MethodType

        def _wrap_add_event_log(proc_self, event_type, data=None, *, _orig=sim.operation_registry):
            OperationProcess.add_event_log(proc_self, event_type, data)
            if event_type == "OPERATION_END":
                conforms, report, _ = sim.effect_engine.validate_now()
                if not conforms:
                    self._spawn_for_violations(sim, report)

        # patch *new* processes as well
        def _patch(proc):
            proc.add_event_log = MethodType(_wrap_add_event_log, proc)

        for proc in sim.operation_registry.values():
            _patch(proc)

        # remember for future spawned ops
        OperationProcess.add_event_log = _wrap_add_event_log  # type: ignore
        OperationProcess.add_event_log_patched = True  # flag

        yield sim.env.timeout(float("inf"))  # keep coroutine alive

    def _run(self, sim: Simulation):
        if self.inspect_interval > 0:
            return self._run_every_dt(sim)
        # inspect_interval == 0
        return self._run_on_every_operation(sim)
