from __future__ import annotations

import math
from enum import Enum
from typing import Any

from .base import Entity


class OperationType(str, Enum):
    """ default operation types """

    TransferContainer = 'TransferContainer'

    TransferLiquid = 'TransferLiquid'

    TransferSolid = 'TransferSolid'

    Heating = 'Heating'

    CleanContainer = 'CleanContainer'

    Purification = 'Purification'

    Concentration = 'Concentration'

    ConcentrationAndPurification = 'ConcentrationAndPurification'

    MakeSolution = 'MakeSolution'


class Operation(Entity):
    """ a process whose temporal consumption can be estimated """

    type: OperationType
    """ the type of the operation """

    process_times: dict[str, float] = dict()
    """ map the functional modules' identifiers to the times required for processing this operation """

    from_reaction: str | None = None
    """ the identifier of the reaction that includes this operation """

    annotations: dict[str, Any] = dict()
    """ annotations for this operation """

    precedents: list[str] = []
    """ identifiers of the preceding operations, these are unrestricted precedence relationships in contrast to 
    those explicitly defined in lmax/lmin """

    @property
    def can_be_processed_by(self) -> list[str]:
        """ the identifiers of the functional modules that can process this operation """
        return sorted([k for k, v in self.process_times.items() if v < 1e8])

    @property
    def finite_process_times(self) -> dict[str, float]:
        return {k: v for k, v in self.process_times.items() if v < math.inf}

class FunctionalModule(Entity):
    """ a functional module is a set of hardware units that can perform one or more operations """

    name: str
    """ name of the functional module, only for annotation purposes """

    capacity: int = 1
    """ can process at most this number of operations at the same time """

    hardware_units: list[str] = []
    """ the constituents of this functional module, only for annotation purposes for now """

    fifo: bool = False
    """ if first-in-first-out policy applies to this machine """

    can_process: list[OperationType]
    """ what types of operations this module can process """
