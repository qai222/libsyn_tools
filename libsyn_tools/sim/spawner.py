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

from .simulation import Simulation, Operation


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
