from __future__ import annotations

import simpy
from pydantic import Field

from libsyn_tools.sim.operation.operation import Operation, StrOrSelector, _OpState
from libsyn_tools.sim.operation.unitary_edit import UnitaryEdit


class Wait(Operation):
    """
    An operation that waits for a duration, optionally locking a resource.
    """

    participant_resource: StrOrSelector | None = Field(default=None)

    def pre_act(self, env: simpy.Environment):
        if self.participant_resource is None:
            if self.sim_state is not _OpState.NEW:
                raise RuntimeError(f"{self.identifier}: pre_act called in state {self.sim_state}")
            self.sim_state = _OpState.PREPARED
            return env.process(self._pre_act_noop())
        return super().pre_act(env)

    def _pre_act_noop(self):
        self.resolved_resources = {}
        self.resources = []
        self.operation_effects = self.get_operation_effects()
        if False:  # pragma: no cover
            yield None

    def get_operation_effects(self) -> list[UnitaryEdit]:
        return []
