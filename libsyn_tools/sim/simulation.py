from __future__ import annotations

import random
from collections import defaultdict
from typing import Dict, List, Optional

import pandas as pd
import simpy
from loguru import logger
from pandas._typing import FilePath
from pydantic import BaseModel

from .effect_engine import EffectEngine
from .knowledge_graph import LabObject
from .operation.operation import Operation


class OperationEventRecord(BaseModel):
    operation_id: str
    timestamp: float
    event_type: str
    operation_data: dict

    def __repr__(self) -> str:  # pragma: no cover – for interactive debug only
        return (
            f"OperationEventRecord("
            f"operation_id={self.operation_id}, "
            f"timestamp={self.timestamp:.3f}, "
            f"event_type={self.event_type})"
        )


class OperationProcess:
    """
    Wraps an Operation so that it can be executed as a SimPy process.
    """

    simpy_process: Optional[simpy.events.Process] = None

    def __init__(
            self,
            *,
            env: simpy.Environment,
            operation: Operation,
            operation_registry: Dict[str, "OperationProcess"],
            dependents: Dict[str, List[str]],
            history_log: List[OperationEventRecord],
            effect_engine: EffectEngine,
            speed_factor: float,
    ):
        self.env = env
        self.operation = operation
        self.operation_registry = operation_registry
        self.dependents = dependents

        self.history_log = history_log
        self.effect_engine = effect_engine
        self.speed_factor = speed_factor

        # Public event that predecessors / dependents can `yield`
        self.done_event = env.event()

    def sim_time(self, dt: float) -> float:
        """Scale a wall-clock delta by the *simulation* speed factor."""
        return dt * self.speed_factor

    def add_event_log(self, event_type: str, data: dict | None = None) -> None:
        self.history_log.append(
            OperationEventRecord(
                operation_id=self.operation.identifier,
                timestamp=self.env.now,
                event_type=event_type,
                operation_data=data or self.operation.model_dump(),
            )
        )

    def run(self):
        """
        High-level lifecycle
        --------------------
        1.  Wait for schedule window + precedent operations
        2.  `operation.pre_act()` → participant resolution & locking
        3.  Simulated execution delay (`temporal_cost`)
        4.  Apply edits atomically (rollback on failure)
        5.  Release locks / cascade interrupts as needed
        """
        try:
            # 1) if scheduled, minimum delay to the scheduled time
            if self.operation.scheduled_start_time is not None:
                delay = self.operation.scheduled_start_time - self.env.now
                if delay > 0:
                    yield self.env.timeout(self.sim_time(delay))

            # wait for all precedents at once
            precedent_events = [self.operation_registry[pid].done_event for pid in self.operation.required_precedents]
            if precedent_events:
                yield simpy.events.AllOf(self.env, precedent_events)

            # 2) Resolve dynamic participants + acquire locks
            yield self.operation.pre_act(self.env)
            # check presumptions
            for pres in self.operation.presumptions:
                if not pres.validate(...):  # supply KG or context
                    raise ValueError(f"Presumption failed for {self.operation.identifier}: {pres}")

            # 3) Log start and simulate intrinsic duration
            self.add_event_log("OPERATION_START")
            logger.debug(f"[t={self.env.now:.2f}] Start {self.operation.__class__.__name__}={self.operation.identifier}")

            if self.operation.temporal_cost:
                yield self.env.timeout(self.sim_time(self.operation.temporal_cost))

            # 4) apply edits atomically
            staged = self.effect_engine.prepare(self.operation)
            self.effect_engine.apply(staged, self.env)

            # Normal completion – mark process done
            self.done_event.succeed()
            self.add_event_log("OPERATION_END")
            logger.debug(f"[t={self.env.now:.2f}] Finished {self.operation.identifier}")

        # -------------------------------------------------------------- #
        # Error / interrupt handling
        # -------------------------------------------------------------- #
        except simpy.Interrupt as intr:  # downstream abort
            logger.warning(f"{self.operation.identifier} interrupted: {intr.cause!r}")
            self.done_event.fail(intr)
            self.add_event_log("OPERATION_ABORTED")

        except Exception as err:
            self.done_event.fail(err)
            self.add_event_log("OPERATION_ERROR")
            logger.error(f"[t={self.env.now:.2f}] {self.operation.identifier} failed: {err!r}")
            self._cascade_interrupt(err)
            raise

        finally:
            # Always release locks obtained in *pre_act*
            self.operation.post_act()

    def _cascade_interrupt(self, cause: Exception):
        """Interrupt all dependent SimPy processes."""
        for dep_id in self.dependents[self.operation.identifier]:
            dep_proc = self.operation_registry[dep_id]
            if dep_proc.simpy_process and dep_proc.simpy_process.is_alive:
                try:
                    dep_proc.simpy_process.interrupt(cause)
                except RuntimeError:  # already terminated
                    pass


class Simulation:
    """
    Orchestrates a collection of `Operation` instances inside a SimPy
    `Environment`.  One `Simulation` ≈ one *experimental run*.
    """

    def __init__(
            self,
            operations: List[Operation],
            *,
            simulation_speed_factor: float = 1.0,
            random_seed: int | None = None,
    ):
        # 0) SimPy env + RNG
        self.env = simpy.Environment()
        self.rng = random.Random(random_seed)

        # 1) Core data
        self.operations = operations
        self.speed_factor = simulation_speed_factor
        self.effect_engine = EffectEngine()

        # 2) Runtime bookkeeping
        self.operation_registry: Dict[str, OperationProcess] = {}
        self.dependents: Dict[str, List[str]] = defaultdict(list)
        self.history_log: List[OperationEventRecord] = []

        # 3) Build runtime artefacts
        self._build_dependency_map()
        self._build_resources()
        self._build_processes()

    # .................................................................. #
    # Construction helpers
    # .................................................................. #
    def _build_dependency_map(self) -> None:
        for op in self.operations:
            for pred in op.required_precedents:
                self.dependents[pred].append(op.identifier)

    def _build_resources(self) -> None:
        """
        For **every** `LabObject` currently known in memory:

        • Create a capacity-1 `simpy.Resource`
        • Store it in *both* the per-simulation map *and*
          the global runtime map used by selectors
        • If the object carries a `pool_type`, drop it into the
          corresponding `FilterStore` so Attribute/History selectors work
          without manual intervention.
        """
        for obj in LabObject.object_lookup.values():
            self.effect_engine._register_if_new(obj, self.env)

    def _build_processes(self) -> None:
        for op in self.operations:
            if op.identifier in self.operation_registry:
                raise ValueError(f"Duplicate Operation identifier {op.identifier}")
            self.operation_registry[op.identifier] = OperationProcess(
                env=self.env,
                operation=op,
                operation_registry=self.operation_registry,
                dependents=self.dependents,
                history_log=self.history_log,
                effect_engine=self.effect_engine,
                speed_factor=self.speed_factor,
            )

    def run(self, until: float | None = None) -> None:
        """Kick off all `OperationProcess` coroutines and block until done."""
        for proc in self.operation_registry.values():
            proc.simpy_process = self.env.process(proc.run())

        logger.info("Simulation start")
        self.env.run(until=until)
        logger.info(f"Simulation end @ t = {self.env.now}")

    def export_event_log(self, filename: FilePath) -> None:
        """Persist the in-memory history to CSV/JSON downstream."""
        df_log = pd.DataFrame.from_records([r.model_dump() for r in self.history_log])
        df_log.to_csv(filename, index=False)
        logger.info(f"Event log exported → {filename}")

    # Convenience factory ------------------------------------------------ #
    @classmethod
    def compile_actions(cls, *actions: Operation, **kwargs) -> "Simulation":
        """Sugar for `Simulation(list(actions), **kwargs)`."""
        return cls(list(actions), **kwargs)
