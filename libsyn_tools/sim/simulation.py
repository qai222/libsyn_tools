from __future__ import annotations

from typing import Dict, List

import pandas as pd
import simpy
from loguru import logger
from pandas._typing import FilePath
from pydantic import BaseModel

from .action import Action, UnitaryEdit
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

    def __init__(
            self,
            env: simpy.Environment,
            action: Action,
            action_registry: Dict[str, ActionProcess],
            resource_map: Dict[str, simpy.Resource],
            history_log: List[ActionEventRecord],
            effect_engine: EffectEngine,
    ):
        """
        :param env: A SimPy Environment.
        :param action: The Action we want to simulate.
        :param action_registry: Maps action.identifier -> ActionProcess.
        :param resource_map: Maps resource IRI -> simpy.Resource for concurrency control.
        :param history_log: A list that will store ActionEventRecord objects.
        """
        self.env = env
        self.action = action
        self.action_registry = action_registry
        self.resource_map = resource_map

        # We store the reference to a shared or global event log
        self.history_log = history_log

        # Signaling event to indicate this process has finished.
        self.done_event = env.event()

        self.effect_engine = effect_engine

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
        # if scheduled, minimum delay to the scheduled time
        if self.action.scheduled_start_time is not None:
            delay = self.action.scheduled_start_time - self.env.now
            if delay > 0:
                yield self.env.timeout(delay)

        # wait for all precedents at once
        precedent_events = [self.action_registry[pid].done_event for pid in self.action.required_precedents]
        if precedent_events:
            try:
                yield self.env.all_of(precedent_events)
            except Exception:  # keep original traceback
                logger.error(f"{self.action.identifier}: precedent failure", exc_info=True)
                raise  # bubbles to outer handler

        # check presumptions
        for pres in self.action.presumptions:
            if not pres.validate(...):  # supply KG or context
                raise ValueError(
                    f"Presumption failed for {self.action.identifier}: {pres}"
                )

        # acquire resource
        resource_reqs: dict[simpy.Resource, simpy.events.Request] = dict()
        for iri in sorted(set(self.action.get_resources())):
            res = self.resource_map.setdefault(iri, simpy.Resource(self.env, 1))
            if res not in resource_reqs:
                resource_reqs[res] = res.request()
        try:
            yield self.env.all_of(resource_reqs.values())
        except Exception as err:
            # acquire failure: mark action failed for dependents
            self.done_event.fail(err)
            raise

        applied_edits: List[UnitaryEdit] = []
        try:
            # Log ACTION_START
            self.add_event_log(event_type="ACTION_START")
            logger.debug(f"[t={self.env.now:.2f}] Start {self.action.identifier}")

            # simulated execution delay
            if self.action.temporal_cost:
                yield self.env.timeout(self.action.temporal_cost)

            # Stage → apply → finalize
            staged = self.effect_engine.prepare(self.action)
            self.effect_engine.apply(staged)
            applied_edits.extend(staged)
            self.effect_engine.finalize(self.action)

            # Normal completion
            self.done_event.succeed()
            self.add_event_log("ACTION_END")
            logger.debug(f"[t={self.env.now:.2f}] Finished {self.action.identifier}")

        except Exception as err:
            # Failure path: rollback & propagate
            self.effect_engine.rollback(applied_edits)
            self.done_event.fail(err)
            self.add_event_log("ACTION_ERROR")
            logger.error(f"[t={self.env.now:.2f}] {self.action.identifier} aborted: {err!r}")
            raise

        finally:
            for req in resource_reqs.values():
                if req.triggered:
                    req.resource.release(req)


class ActionSimulation:
    def __init__(self, actions: List[Action]):
        self.env = simpy.Environment()
        self.actions = actions
        self.action_registry: Dict[str, ActionProcess] = {}
        self.resource_map: Dict[str, simpy.Resource] = {}

        # A single list to store all event records from the entire simulation
        self.history_log: List[ActionEventRecord] = []

        self.effect_engine = EffectEngine()

    def build_processes(self):
        for act in self.actions:
            if act.identifier in self.action_registry:
                raise ValueError(f"Duplicate Action identifier {act.identifier}.")
            self.action_registry[act.identifier] = ActionProcess(
                env=self.env,
                action=act,
                action_registry=self.action_registry,
                resource_map=self.resource_map,
                history_log=self.history_log,  # pass the shared log
                effect_engine=self.effect_engine
            )

    def schedule_all(self):
        for proc in self.action_registry.values():
            self.env.process(proc.run())

    def run(self, until: float = None):
        logger.info("Starting the simulation environment.")
        self.build_processes()
        self.schedule_all()
        self.env.run(until=until)
        logger.info(f"Simulation ended at time {self.env.now}")

    def export_event_log(self, filename: FilePath):
        """
        Example function to export the event log to a CSV or JSON file.
        """

        df_log = pd.DataFrame.from_records([r.model_dump() for r in self.history_log])
        df_log.to_csv(filename, index=False)
        logger.info(f"Event log exported to {filename}")
