from __future__ import annotations
from typing import Dict, List

import pandas as pd
import simpy
from loguru import logger
from pandas._typing import FilePath
from pydantic import BaseModel

from .action import Action


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

    def run(self):
        """
        The main generator function describing the logic for running an Action in the simulation.
        """
        # 1) Wait until scheduled_start_time if specified
        if self.action.scheduled_start_time is not None:
            delay = self.action.scheduled_start_time - self.env.now
            if delay > 0:
                yield self.env.timeout(delay)

        # 2) Wait for all required precedents to finish
        for precedent_id in self.action.required_precedents:
            precedent_proc = self.action_registry.get(precedent_id)
            if not precedent_proc:
                raise ValueError(
                    f"Action {self.action.identifier} declares a required precedent {precedent_id}, "
                    "but no such action exists in the simulation registry."
                )
            yield precedent_proc.done_event

        # 3) Request needed resources
        # dynamically get resources
        self.action.resources = self.action.get_resources()
        requests = []
        for resource_iri in self.action.resources:
            if resource_iri not in self.resource_map:
                self.resource_map[resource_iri] = simpy.Resource(self.env, capacity=1)
            req = self.resource_map[resource_iri].request()
            requests.append(req)
        yield self.env.all_of(requests)

        # 4) The action is about to start:
        self.history_log.append(
            ActionEventRecord(
                action_id=self.action.identifier,
                timestamp=self.env.now,
                event_type="ACTION_START",
                action_data=self.action.model_dump(),
            )
        )
        logger.debug(f"[t={self.env.now:.2f}] Starting action {self.action.identifier}.")
        self.action.pre_act()
        # TODO double check if we should run preact here
        # TODO a better way may be run preact twice (before and after env.timeout and compare if they are identical)
        #  since we are assuming the (inferred) effects stay unchanged before and after timeout

        # 5) "Execute" the action, which might have a temporal_cost
        duration = self.action.temporal_cost or 0.0
        if duration > 0:
            yield self.env.timeout(duration)

        # 6) Actually apply the UnitaryEdits, logging each one:
        for edit in self.action.action_effects:
            # Log the application of each edit
            edit.apply()
            self.history_log.append(
                ActionEventRecord(
                    action_id=self.action.identifier,
                    timestamp=self.env.now,
                    event_type="EFFECT_APPLIED",
                    action_data=self.action.model_dump(),
                    # details={"unitary_edit_type": edit.type, "unitary_edit_data": edit.model_dump()},
                )
            )

        # 7) Action is finished; release the resources
        for i, rsrc in enumerate(self.action.resources):
            self.resource_map[rsrc].release(requests[i])

        # 8) Mark action as complete
        self.done_event.succeed()
        self.history_log.append(
            ActionEventRecord(
                action_id=self.action.identifier,
                timestamp=self.env.now,
                event_type="ACTION_END",
                action_data=self.action.model_dump(),
            )
        )
        logger.debug(f"[t={self.env.now:.2f}] Finished action {self.action.identifier}.")
        self.action.post_act()


class ActionSimulation:
    def __init__(self, actions: List[Action]):
        self.env = simpy.Environment()
        self.actions = actions
        self.action_registry: Dict[str, ActionProcess] = {}
        self.resource_map: Dict[str, simpy.Resource] = {}

        # A single list to store all event records from the entire simulation
        self.history_log: List[ActionEventRecord] = []

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
