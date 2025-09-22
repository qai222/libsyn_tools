from __future__ import annotations

from typing import Callable, List, Any


class LifecycleCallbacks:
    """A minimal pub/sub surface for simulation lifecycle events."""

    def __init__(self) -> None:
        self.on_operation_start: List[Callable[[Any], None]] = []
        self.on_operation_end: List[Callable[[Any], None]] = []
        self.on_operation_interrupt: List[Callable[[Any, str], None]] = []
        self.on_violation: List[Callable[[Any], None]] = []

    # Emitters ---------------------------------------------------------

    def emit_operation_start(self, proc: Any) -> None:
        for cb in list(self.on_operation_start):
            cb(proc)

    def emit_operation_end(self, proc: Any) -> None:
        for cb in list(self.on_operation_end):
            cb(proc)

    def emit_operation_interrupt(self, proc: Any, reason: str) -> None:
        for cb in list(self.on_operation_interrupt):
            cb(proc, reason)

    def emit_violation(self, record: Any) -> None:
        for cb in list(self.on_violation):
            cb(record)
