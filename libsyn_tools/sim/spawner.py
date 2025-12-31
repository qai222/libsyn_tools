"""
Spawners for endogenous event generation — without monkey-patching.

Refactor summary
----------------
- Spawners now **subscribe** to `Simulation.callbacks` lifecycle hooks
  rather than patching `OperationProcess.add_event_log`.
- TimerSpawner still runs as a SimPy coroutine.
- KGInspectorSpawner can validate either (a) after every OPERATION_END
  (inspect_interval==0) or (b) periodically (inspect_interval>0).
- ProcessInterruptSpawner subscribes to `on_operation_interrupt`.

Importantly, no method monkey-patching remains here.
"""

from __future__ import annotations

import weakref
from abc import ABC, abstractmethod
from typing import Optional, Callable, Any, Dict

from loguru import logger
from pydantic import BaseModel, Field, PrivateAttr
from rdflib.namespace import SH, RDF

from .simulation import Simulation, Operation, OperationProcess
from .knowledge_graph import identifier_from_iri


class Spawner(BaseModel, ABC):
    _sim_ref: Optional[weakref.ReferenceType] = PrivateAttr(default=None)

    def attach(self, sim: Simulation) -> None:
        if self._sim_ref is not None:
            raise RuntimeError("Spawner already attached")
        self._sim_ref = weakref.ref(sim)
        self._on_attach(sim)

    @property
    def sim(self) -> Simulation:
        s = None if self._sim_ref is None else self._sim_ref()
        if s is None:
            raise RuntimeError("Spawner has not been attached yet")
        return s

    @staticmethod
    def _sim_time(sim: Simulation, dt: float) -> float:
        """Apply simulation_speed_factor to a base dt for spawner timeouts."""
        return dt * sim.speed_factor

    @abstractmethod
    def _on_attach(self, sim: Simulation) -> None:
        """Install callbacks or processes into the Simulation."""

    @abstractmethod
    def _detach_safe(self, sim: Simulation) -> None:
        """Optional: remove callbacks (not strictly required for most runs)."""


class TimerSpawner(Spawner):
    op_factory: Callable[[Any], Operation]
    interval: float = Field(..., description="delta t in base time (scaled by simulation_speed_factor)")
    start_offset: float = 0.0
    alive: bool = True

    def cancel(self) -> None:
        self.alive = False

    def _sample_dt(self) -> float:
        return float(self.interval)

    def _on_attach(self, sim: Simulation) -> None:
        sim.env.process(self._run(sim))

    def _detach_safe(self, sim: Simulation) -> None:
        self.alive = False

    def _run(self, sim: Simulation):
        env = sim.env
        if self.start_offset > 0:
            yield env.timeout(self._sim_time(sim, self.start_offset))
        while self.alive:
            op = self.op_factory(sim)
            sim.spawn_operation(op)
            yield env.timeout(self._sim_time(sim, self._sample_dt()))


class KGInspectorSpawner(Spawner):
    """
    Validate the KG against a SHACL sub-set and spawn corrective Operations.

    shape_dispatch : dict[str, Callable[[str], Operation]]
        Mapping **sourceShape IRI ➜ factory**. The factory receives the
        *focusNode IRI* and must return a fully constructed Operation
        ready for spawn_operation(). Return None to ignore.
    inspect_interval : float
        • > 0   → run every `inspect_interval` seconds of base time
                 (scaled by simulation_speed_factor)
        • == 0  → run right after every OPERATION_END (via callbacks)
    """

    shape_dispatch: Dict[str, Callable[[str], Optional[Operation]]]
    inspect_interval: float = Field(0.0, ge=0.0)
    _subscribed: bool = PrivateAttr(default=False)

    def _spawn_for_violations(self, sim: Simulation, report) -> None:
        g = report
        for vr in g.subjects(RDF.type, SH.ValidationResult):
            shape_iri = str(g.value(vr, SH.sourceShape))
            focus_iri = identifier_from_iri(str(g.value(vr, SH.focusNode)))
            logger.critical(f"DBG sourceShape = {shape_iri}")
            factory = self.shape_dispatch.get(shape_iri)
            if factory is None:
                continue
            op = factory(focus_iri)
            if op is not None:
                sim.spawn_operation(op)

    def _on_attach(self, sim: Simulation) -> None:
        if self.inspect_interval > 0:
            sim.env.process(self._run_every_dt(sim))
        else:
            # subscribe to operation_end lifecycle
            if not self._subscribed:
                def _on_end(proc: OperationProcess):
                    conforms, report, _ = sim.effect_engine.validate_now()
                    if not conforms:
                        self._spawn_for_violations(sim, report)

                sim.callbacks.on_operation_end.append(_on_end)
                self._on_end = _on_end  # keep reference for detach
                self._subscribed = True

    def _detach_safe(self, sim: Simulation) -> None:
        if self._subscribed and hasattr(self, "_on_end"):
            try:
                sim.callbacks.on_operation_end.remove(self._on_end)
            except ValueError:
                pass
        self._subscribed = False

    def _run_every_dt(self, sim: Simulation):
        env = sim.env
        dt = self.inspect_interval
        while True:
            yield env.timeout(self._sim_time(sim, dt))
            conforms, report, _ = sim.effect_engine.validate_now()
            if not conforms:
                self._spawn_for_violations(sim, report)


class ProcessInterruptSpawner(Spawner):
    """
    React to `simpy.Interrupt` events that abort running operations.

    interrupt_dispatch : dict[str, Callable[[OperationProcess, str], Optional[Operation]]]
        Mapping **reason string ➜ factory**.  The factory receives
        `(proc, reason)` and must return a fully constructed Operation
        or None to ignore the interrupt.
    """

    interrupt_dispatch: Dict[str, Callable[[OperationProcess, str], Optional[Operation]]]

    def _on_attach(self, sim: Simulation) -> None:
        def _on_interrupt(proc: OperationProcess, reason: str):
            factory = self.interrupt_dispatch.get(str(reason))
            if factory is None:
                return
            op = factory(proc, reason)
            if op is not None:
                self.sim.spawn_operation(op)

        sim.callbacks.on_operation_interrupt.append(_on_interrupt)
        self._on_interrupt = _on_interrupt

    def _detach_safe(self, sim: Simulation) -> None:
        if hasattr(self, "_on_interrupt"):
            try:
                sim.callbacks.on_operation_interrupt.remove(self._on_interrupt)
            except ValueError:
                pass
