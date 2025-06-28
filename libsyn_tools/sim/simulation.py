from __future__ import annotations

import random
from collections import defaultdict
from typing import Dict, List, Optional

import pandas as pd
import simpy
from loguru import logger
from pandas._typing import FilePath
from pydantic import BaseModel

from .action import Action
from .effect_engine import EffectEngine


class ActionEventRecord(BaseModel):
    action_id: str
    timestamp: float
    event_type: str
    action_data: dict

    def __repr__(self):
        return (
            f"ActionEventRecord(action_id={self.action_id}, timestamp={self.timestamp}, event_type={self.event_type})"
        )


class ActionProcess:
    """
    Wraps an Action so that it can be executed as a SimPy process.
    """

    simpy_process: Optional[simpy.events.Process] = None

    def __init__(
            self,
            *,
            env: simpy.Environment,
            action: Action,
            action_registry: Dict[str, ActionProcess],
            dependents: Dict[str, List[str]],  # pre‑computed dependency map
            resource_map: Dict[str, simpy.Resource],
            history_log: List[ActionEventRecord],
            effect_engine: EffectEngine,
            speed_factor: float,
    ):
        """
        :param env: A SimPy Environment.
        :param action: The Action we want to simulate.
        :param action_registry: Maps action.identifier -> ActionProcess.
        :param resource_map: Maps resource IRI -> simpy.Resource for concurrency control.
        :param history_log: A list that will store ActionEventRecord objects.
        """
        self.dependents = dependents
        self.env = env
        self.action = action
        self.action_registry = action_registry
        self.resource_map = resource_map

        # We store the reference to a shared or global event log
        self.history_log = history_log

        # Signaling event to indicate this process has finished.
        self.done_event = env.event()

        self.effect_engine = effect_engine

        self.speed_factor = speed_factor

    def sim_time(self, dt: float) -> float:
        """Scale `dt` by the global *simulation_speed_factor* (P0‑5)."""
        return dt * self.speed_factor

    def add_event_log(self, event_type: str, data: dict | None = None) -> None:
        self.history_log.append(
            ActionEventRecord(
                action_id=self.action.identifier,
                timestamp=self.env.now,
                event_type=event_type,
                action_data=data or self.action.model_dump(),
            )
        )

    def run(self):
        """
        The main generator function describing the logic for running an Action in the simulation.
        """

        resource_reqs: dict[simpy.Resource, simpy.events.Request] = dict()
        try:

            # if scheduled, minimum delay to the scheduled time
            if self.action.scheduled_start_time is not None:
                delay = self.action.scheduled_start_time - self.env.now
                if delay > 0:
                    yield self.env.timeout(self.sim_time(delay))

            # wait for all precedents at once
            precedent_events = [self.action_registry[pid].done_event for pid in self.action.required_precedents]
            if precedent_events:
                yield self.env.all_of(precedent_events)

            # check presumptions
            for pres in self.action.presumptions:
                if not pres.validate(...):  # supply KG or context
                    raise ValueError(f"Presumption failed for {self.action.identifier}: {pres}")

            # acquire resource
            for iri in sorted(set(self.action.get_resources())):
                try:
                    res = self.resource_map[iri]
                except KeyError as exc:
                    raise RuntimeError(f"Resource {iri} is invalid; build_resources incomplete") from exc

                if res not in resource_reqs:
                    # choose request type based on resource class
                    resource_reqs[res] = res.request()
            yield self.env.all_of(resource_reqs.values())

            # action start
            self.add_event_log(event_type="ACTION_START")
            logger.debug(f"[t={self.env.now:.2f}] Start {self.action.identifier}")

            # simulated execution delay
            if self.action.temporal_cost:
                yield self.env.timeout(self.sim_time(self.action.temporal_cost))

            # Stage → apply → finalize
            staged = self.effect_engine.prepare(self.action)
            self.effect_engine.apply(staged)
            self.effect_engine.finalize(self.action)

            # Normal completion
            self.done_event.succeed()
            self.add_event_log("ACTION_END")
            logger.debug(f"[t={self.env.now:.2f}] Finished {self.action.identifier}")

        except simpy.Interrupt as intr:  # downstream abort (P0‑2)
            logger.warning(f"{self.action.identifier} interrupted: {intr.cause!r}")
            self.done_event.fail(intr)
            self.add_event_log("ACTION_ABORTED")

        except Exception as err:
            self.done_event.fail(err)
            self.add_event_log("ACTION_ERROR")
            logger.error(f"[t={self.env.now:.2f}] {self.action.identifier} failed: {err!r}")
            self._cascade_interrupt(err)
            raise
        finally:
            for req in resource_reqs.values():
                if req.triggered:
                    req.resource.release(req)

    def _cascade_interrupt(self, cause: Exception):
        """Interrupt all dependent SimPy processes (P0‑2)."""
        for dep_id in self.dependents[self.action.identifier]:
            dep_proc = self.action_registry[dep_id]
            if dep_proc.simpy_process and dep_proc.simpy_process.is_alive:
                try:
                    dep_proc.simpy_process.interrupt(cause)
                except RuntimeError:  # process already terminated
                    pass


class ActionSimulation:
    def __init__(
            self,
            actions: List[Action],
            *,
            simulation_speed_factor: float = 1.0,
            random_seed: int | None = None,
    ):
        self.env = simpy.Environment()
        self.actions = actions

        self.rng = random.Random(random_seed)
        self.speed_factor = simulation_speed_factor

        self.action_registry: Dict[str, ActionProcess] = {}
        self.resource_map: Dict[str, simpy.Resource] = {}
        # A single list to store all event records from the entire simulation
        self.history_log: List[ActionEventRecord] = []

        self.effect_engine = EffectEngine()
        self.dependents: Dict[str, List[str]] = defaultdict(list)

        self.build_dependency_map()
        self.build_resources()
        self.build_processes()

    def build_dependency_map(self):
        for act in self.actions:
            for pred in act.required_precedents:
                self.dependents[pred].append(act.identifier)

    def build_resources(self):
        for act in self.actions:
            for iri in act.get_resources():
                if iri not in self.resource_map:
                    self.resource_map[iri] = simpy.Resource(self.env, capacity=1)

    def build_processes(self):
        for act in self.actions:
            if act.identifier in self.action_registry:
                raise ValueError(f"Duplicate Action identifier {act.identifier}.")
            self.action_registry[act.identifier] = ActionProcess(
                env=self.env,
                action=act,
                action_registry=self.action_registry,
                dependents=self.dependents,
                resource_map=self.resource_map,
                history_log=self.history_log,  # pass the shared log
                effect_engine=self.effect_engine,
                speed_factor=self.speed_factor,
            )

    def run(self, until: float = None):
        for proc in self.action_registry.values():
            proc.simpy_process = self.env.process(proc.run())

        logger.info("Starting the simulation environment.")
        self.env.run(until=until)
        logger.info(f"Simulation ended at time: {self.env.now}")

    def export_event_log(self, filename: FilePath):
        """
        Example function to export the event log to a CSV or JSON file.
        """
        df_log = pd.DataFrame.from_records([r.model_dump() for r in self.history_log])
        df_log.to_csv(filename, index=False)
        logger.info(f"Event log exported to: {filename}")

    @classmethod
    def compile_actions(cls, *actions: "Action", **kwargs) -> "ActionSimulation":
        """
        Factory wrapper that instantiates :class:`ActionSimulation` directly.
        """
        return cls(list(actions), **kwargs)

    def get_resource_pool(self, iris: list[str]) -> list[simpy.Resource]:
        """
        Return the *Resource* objects corresponding to *iris* (ordered).
        """
        return [self.resource_map[i] for i in iris]
