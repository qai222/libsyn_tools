"""
Spawners for endogenous event generation.
A spawner attaches itself to one `Simulation`,
runs as a SimPy coroutine and creates new `Operation`s when its own trigger fires.
"""

from __future__ import annotations

import weakref
from abc import ABC, abstractmethod
from types import MethodType
from typing import Optional, Callable, Any
from loguru import logger

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
            logger.critical(f"DBG sourceShape = {shape_iri}")
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
        original_add = OperationProcess.add_event_log

        def _wrap(proc_self, event_type, data=None, *, _orig=original_add):
            _orig(proc_self, event_type, data)  # call previous impl.
            if event_type == "OPERATION_END":
                conforms, report, _ = sim.effect_engine.validate_now()
                if not conforms:
                    self._spawn_for_violations(sim, report)

        # class-level patch
        OperationProcess.add_event_log = _wrap
        OperationProcess._kginsp_patched = True

        # re-bind every existing instance
        for p in sim.operation_registry.values():
            p.add_event_log = MethodType(_wrap, p)

        yield sim.env.timeout(float("inf"))

    def _run(self, sim: Simulation):
        if self.inspect_interval > 0:
            yield from self._run_every_dt(sim)
        # inspect_interval == 0
        yield from self._run_on_every_operation(sim)


class ProcessInterruptSpawner(Spawner):
    """
    Listen for `simpy.Interrupt` events that abort running operations and
    launch corrective or resumption Operations.

    Parameters
    ----------
    interrupt_dispatch : dict[str, Callable[[OperationProcess, str], Operation]]
        Mapping **reason string ➜ factory**.
        * The *reason* is what the interrupted process stored under
          `data["reason"]` when it called
          `add_event_log("OPERATION_INTERRUPT", {"reason": ...})`.
        * The factory receives `(proc, reason)` and must return a fully
          constructed `Operation` ready for `sim.spawn_operation()`.
          Return `None` to ignore the interrupt.
    """

    interrupt_dispatch: dict[str, Callable[[OperationProcess, str], Optional[Operation]]]

    def _install_wrapper(self, sim: Simulation):
        """
        Monkey-patch `OperationProcess.add_event_log` so that every time an
        interrupt is logged we can react.  The wrapper is installed **once**
        per Python process, no matter how many spawners are attached.
        """
        if getattr(OperationProcess, "_interrupt_wrapper_installed", False):
            return  # somebody else already patched

        original_add = OperationProcess.add_event_log
        spawner_ref = weakref.ref(self)  # avoid cycles

        def _wrapped(proc_self, event_type, data=None, *, _orig=original_add):
            # keep original behaviour first
            _orig(proc_self, event_type, data)

            if event_type != "OPERATION_INTERRUPT":
                return

            spawner = spawner_ref()
            if spawner is None:  # spawner GC'ed
                return

            reason = str((data or {}).get("reason", ""))
            factory = spawner.interrupt_dispatch.get(reason)
            if factory is None:
                return  # reason not mapped

            op = factory(proc_self, reason)
            if op is not None:
                spawner.sim.spawn_operation(op)

        # class-level patch (affects future OperationProcess instances)
        OperationProcess.add_event_log = _wrapped
        OperationProcess._interrupt_wrapper_installed = True  # flag

        # also patch every *existing* instance so they go through wrapper
        for proc in sim.operation_registry.values():
            proc.add_event_log = MethodType(_wrapped, proc)

    def _run(self, sim: Simulation):
        # one-time installation, then just keep the coroutine alive
        self._install_wrapper(sim)
        yield sim.env.timeout(float("inf"))  # dormant forever
