from __future__ import annotations

from typing import Callable, List, Any, Optional


class LifecycleCallbacks:
    """A minimal pub/sub surface for simulation lifecycle events."""

    def __init__(
        self,
        *,
        history_logger: Optional[Callable[[str, float, str, dict], None]] = None,
    ) -> None:
        self.on_operation_start: List[Callable[[Any], None]] = []
        self.on_operation_end: List[Callable[[Any], None]] = []
        self.on_operation_interrupt: List[Callable[[Any, str], None]] = []
        self.on_violation: List[Callable[[Any], None]] = []
        self._history_logger = history_logger

    def _log_callback_exception(
        self,
        *,
        operation_id: str,
        timestamp: float,
        event_name: str,
        callback: Callable[..., Any],
        exc: Exception,
    ) -> None:
        if self._history_logger is None:
            return
        callback_name = getattr(callback, "__name__", repr(callback))
        self._history_logger(
            operation_id,
            timestamp,
            "CALLBACK_ERROR",
            {
                "event": event_name,
                "callback": callback_name,
                "error_type": type(exc).__name__,
                "error": str(exc),
            },
        )

    # Emitters ---------------------------------------------------------

    def emit_operation_start(self, proc: Any) -> None:
        for cb in list(self.on_operation_start):
            try:
                cb(proc)
            except Exception as exc:
                self._log_callback_exception(
                    operation_id=getattr(proc.operation, "identifier", "UNKNOWN"),
                    timestamp=getattr(proc.env, "now", 0.0),
                    event_name="operation_start",
                    callback=cb,
                    exc=exc,
                )

    def emit_operation_end(self, proc: Any) -> None:
        for cb in list(self.on_operation_end):
            try:
                cb(proc)
            except Exception as exc:
                self._log_callback_exception(
                    operation_id=getattr(proc.operation, "identifier", "UNKNOWN"),
                    timestamp=getattr(proc.env, "now", 0.0),
                    event_name="operation_end",
                    callback=cb,
                    exc=exc,
                )

    def emit_operation_interrupt(self, proc: Any, reason: str) -> None:
        for cb in list(self.on_operation_interrupt):
            try:
                cb(proc, reason)
            except Exception as exc:
                self._log_callback_exception(
                    operation_id=getattr(proc.operation, "identifier", "UNKNOWN"),
                    timestamp=getattr(proc.env, "now", 0.0),
                    event_name="operation_interrupt",
                    callback=cb,
                    exc=exc,
                )

    def emit_violation(self, record: Any) -> None:
        for cb in list(self.on_violation):
            try:
                cb(record)
            except Exception as exc:
                self._log_callback_exception(
                    operation_id=getattr(record, "operation_id", "UNKNOWN"),
                    timestamp=getattr(record, "sim_time", 0.0),
                    event_name="violation",
                    callback=cb,
                    exc=exc,
                )
