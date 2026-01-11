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
import math
from typing import Optional, Callable, Any, Dict, Tuple, Set

import simpy

from loguru import logger
from pydantic import BaseModel, Field, PrivateAttr, field_validator
from rdflib.namespace import SH, RDF
from rdflib.term import URIRef

from .simulation import Simulation, Operation, OperationProcess
from .knowledge_graph import identifier_from_iri
from .effect_engine import KnowledgeGraph
from .effect_shacl import SHACLViolationRecord, _iter_validation_results, _first


class Spawner(BaseModel, ABC):
    _sim_ref: Optional[weakref.ReferenceType] = PrivateAttr(default=None)
    _alive: bool = PrivateAttr(default=True)
    _process: Optional[simpy.events.Process] = PrivateAttr(default=None)

    def attach(self, sim: Simulation) -> None:
        if self._sim_ref is not None:
            raise RuntimeError("Spawner already attached")
        self._sim_ref = weakref.ref(sim)
        if hasattr(sim, "_register_spawner"):
            sim._register_spawner(self)
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

    def detach(self) -> None:
        sim = None if self._sim_ref is None else self._sim_ref()
        if sim is None:
            self._sim_ref = None
            return
        self._detach_safe(sim)
        self._sim_ref = None

    @property
    def requires_until(self) -> bool:
        return False


class TimerSpawner(Spawner):
    op_factory: Callable[[Any], Operation]
    interval: float = Field(..., description="delta t in base time (scaled by simulation_speed_factor)")
    start_offset: float = 0.0
    alive: bool = True

    @field_validator("interval")
    @classmethod
    def _validate_interval(cls, value: float) -> float:
        try:
            is_finite = math.isfinite(value)
        except TypeError as exc:
            raise ValueError("interval must be finite and > 0") from exc
        if not is_finite or value <= 0:
            raise ValueError("interval must be finite and > 0")
        return value

    @field_validator("start_offset")
    @classmethod
    def _validate_start_offset(cls, value: float) -> float:
        try:
            is_finite = math.isfinite(value)
        except TypeError as exc:
            raise ValueError("start_offset must be finite and >= 0") from exc
        if not is_finite or value < 0:
            raise ValueError("start_offset must be finite and >= 0")
        return value

    def cancel(self) -> None:
        self.alive = False

    def _sample_dt(self) -> float:
        return float(self.interval)

    def _on_attach(self, sim: Simulation) -> None:
        self._alive = True
        self._process = sim.env.process(self._run(sim))

    def _detach_safe(self, sim: Simulation) -> None:
        self.alive = False
        self._alive = False
        if self._process is not None:
            try:
                if not self._process.triggered:
                    self._process.interrupt("TimerSpawner detached")
            except Exception:
                pass
            self._process = None

    @property
    def requires_until(self) -> bool:
        return True

    def _run(self, sim: Simulation):
        env = sim.env
        try:
            if self.start_offset > 0:
                yield env.timeout(self._sim_time(sim, self.start_offset))
            while self.alive and self._alive:
                try:
                    op = self.op_factory(sim)
                    sim.spawn_operation(op)
                except Exception as err:
                    logger.warning(f"TimerSpawner failed to spawn operation: {err}")
                yield env.timeout(self._sim_time(sim, self._sample_dt()))
        except simpy.Interrupt:
            return


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

    Note: PolicyEnforcerSpawner + ValidationAuditSpawner is the preferred
    closed-loop composition for violation-driven remediation.
    """

    shape_dispatch: Dict[str, Callable[[str], Optional[Operation]]]
    inspect_interval: float = Field(0.0, ge=0.0)
    _subscribed: bool = PrivateAttr(default=False)

    @field_validator("inspect_interval")
    @classmethod
    def _validate_inspect_interval(cls, value: float) -> float:
        try:
            is_finite = math.isfinite(value)
        except TypeError as exc:
            raise ValueError("inspect_interval must be finite and >= 0") from exc
        if not is_finite or value < 0:
            raise ValueError("inspect_interval must be finite and >= 0")
        return value

    def _spawn_for_violations(self, sim: Simulation, report) -> None:
        g = report
        for vr in g.subjects(RDF.type, SH.ValidationResult):
            shape_node = g.value(vr, SH.sourceShape)
            focus_node = g.value(vr, SH.focusNode)
            if shape_node is None or focus_node is None:
                logger.debug(
                    "Skipping SHACL violation with missing sourceShape/focusNode",
                )
                continue
            if not isinstance(shape_node, URIRef):
                logger.debug(
                    "Skipping SHACL violation with non-IRI sourceShape",
                )
                continue
            if not isinstance(focus_node, URIRef):
                logger.debug(
                    "Skipping SHACL violation with non-IRI focusNode",
                )
                continue
            shape_iri = str(shape_node)
            focus_iri = identifier_from_iri(str(focus_node))
            focus_obj = KnowledgeGraph.get_object_from_lookup(focus_iri)
            if focus_obj is None:
                focus_obj = KnowledgeGraph.get_object_from_lookup(str(focus_node))
            if focus_obj is None:
                logger.debug(
                    "Skipping SHACL violation with unknown focusNode",
                )
                continue
            focus_iri = focus_obj.identifier
            # Debug logging for shape dispatch; keep at debug level to avoid
            # polluting normal test/output runs.
            logger.debug(f"sourceShape = {shape_iri}")
            factory = self.shape_dispatch.get(shape_iri)
            if factory is None:
                continue
            try:
                op = factory(focus_iri)
                if op is not None:
                    op.remediation = True
                    sim.spawn_operation(op)
            except Exception as err:
                logger.warning(f"KGInspectorSpawner factory/spawn failed: {err}")

    def _on_attach(self, sim: Simulation) -> None:
        if self.inspect_interval > 0:
            self._alive = True
            self._process = sim.env.process(self._run_every_dt(sim))
        else:
            # subscribe to operation_end lifecycle
            if not self._subscribed:
                def _on_end(proc: OperationProcess):
                    try:
                        conforms, report, _ = sim.effect_engine.validate_now()
                        if not conforms:
                            self._spawn_for_violations(sim, report)
                    except Exception as err:
                        logger.warning(f"KGInspectorSpawner validate_now failed: {err}")

                sim.callbacks.on_operation_end.append(_on_end)
                self._on_end = _on_end  # keep reference for detach
                self._subscribed = True

    def _detach_safe(self, sim: Simulation) -> None:
        self._alive = False
        if self._subscribed and hasattr(self, "_on_end"):
            try:
                sim.callbacks.on_operation_end.remove(self._on_end)
            except ValueError:
                pass
        self._subscribed = False
        if self._process is not None:
            try:
                if not self._process.triggered:
                    self._process.interrupt("KGInspectorSpawner detached")
            except Exception:
                pass
            self._process = None

    @property
    def requires_until(self) -> bool:
        return self.inspect_interval > 0

    def _run_every_dt(self, sim: Simulation):
        env = sim.env
        dt = self.inspect_interval
        try:
            while self._alive:
                yield env.timeout(self._sim_time(sim, dt))
                try:
                    conforms, report, _ = sim.effect_engine.validate_now()
                    if not conforms:
                        self._spawn_for_violations(sim, report)
                except Exception as err:
                    logger.warning(f"KGInspectorSpawner validate_now failed: {err}")
        except simpy.Interrupt:
            return


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
            try:
                factory = self.interrupt_dispatch.get(str(reason))
                if factory is None:
                    return
                op = factory(proc, reason)
                if op is not None:
                    op.remediation = True
                    self.sim.spawn_operation(op)
            except Exception as err:
                logger.warning(f"ProcessInterruptSpawner factory/spawn failed: {err}")

        sim.callbacks.on_operation_interrupt.append(_on_interrupt)
        self._on_interrupt = _on_interrupt

    def _detach_safe(self, sim: Simulation) -> None:
        if hasattr(self, "_on_interrupt"):
            try:
                sim.callbacks.on_operation_interrupt.remove(self._on_interrupt)
            except ValueError:
                pass


class PolicyEnforcerSpawner(Spawner):
    """
    React to SHACLViolationRecord events and spawn remediation operations.

    shape_dispatch : dict[str, Callable[[SHACLViolationRecord], Optional[Operation]]]
        Mapping **shape IRI ➜ factory**. The factory receives the violation record.
    """

    shape_dispatch: Dict[str, Callable[[SHACLViolationRecord], Optional[Operation]]]
    origins: Tuple[str, ...] = ("SHACL",)
    dispositions: Tuple[str, ...] = ("committed", "aborted")
    dedupe: bool = True
    max_remediations_per_focus: int = 1
    _seen_violation_ids: Set[str] = PrivateAttr(default_factory=set)
    _seen_op_shapes: Set[Tuple[str, Optional[str]]] = PrivateAttr(default_factory=set)
    _remediation_counts: Dict[Tuple[str, Optional[str]], int] = PrivateAttr(default_factory=dict)

    def _should_skip(self, record: SHACLViolationRecord) -> bool:
        if record.origin not in self.origins:
            return True
        if record.disposition not in self.dispositions:
            return True
        if record.shape_iri is None:
            return True
        op_shape = (record.shape_iri, record.focus_iri)
        if self.max_remediations_per_focus >= 0:
            count = self._remediation_counts.get(op_shape, 0)
            if count >= self.max_remediations_per_focus:
                return True
        if not self.dedupe:
            return False
        if record.violation_id in self._seen_violation_ids:
            return True
        if op_shape in self._seen_op_shapes:
            return True
        return False

    def _mark_seen(self, record: SHACLViolationRecord) -> None:
        if record.shape_iri is not None:
            op_shape = (record.shape_iri, record.focus_iri)
            self._remediation_counts[op_shape] = self._remediation_counts.get(op_shape, 0) + 1
            if self.dedupe:
                self._seen_op_shapes.add(op_shape)
        if self.dedupe:
            self._seen_violation_ids.add(record.violation_id)

    def _on_attach(self, sim: Simulation) -> None:
        def _on_violation(record: SHACLViolationRecord):
            try:
                if self._should_skip(record):
                    return
                factory = self.shape_dispatch.get(record.shape_iri)
                if factory is None:
                    return
                op = factory(record)
                if op is None:
                    return
                op.remediation = True
                self._mark_seen(record)
                precedents = None
                if record.operation_id in sim.operation_registry:
                    precedents = [record.operation_id]
                sim.spawn_operation(op, precedents=precedents)
            except Exception as err:
                logger.warning(f"PolicyEnforcerSpawner factory/spawn failed: {err}")

        sim.callbacks.on_violation.append(_on_violation)
        self._on_violation = _on_violation

    def _detach_safe(self, sim: Simulation) -> None:
        if hasattr(self, "_on_violation"):
            try:
                sim.callbacks.on_violation.remove(self._on_violation)
            except ValueError:
                pass


class ValidationAuditSpawner(Spawner):
    """
    Periodically audit SHACL and emit SHACLViolationRecord events.

    inspect_interval : float
        • > 0   → run every `inspect_interval` seconds of base time
        • == 0  → run right after every OPERATION_END (via callbacks)
    """

    inspect_interval: float = Field(0.0, ge=0.0)
    _subscribed: bool = PrivateAttr(default=False)

    @field_validator("inspect_interval")
    @classmethod
    def _validate_inspect_interval(cls, value: float) -> float:
        try:
            is_finite = math.isfinite(value)
        except TypeError as exc:
            raise ValueError("inspect_interval must be finite and >= 0") from exc
        if not is_finite or value < 0:
            raise ValueError("inspect_interval must be finite and >= 0")
        return value

    def _emit_violation_records(self, sim: Simulation, report) -> None:
        for vr in _iter_validation_results(report):
            shape_node = report.value(vr, SH.sourceShape)
            focus_node = report.value(vr, SH.focusNode)
            if shape_node is None or focus_node is None:
                logger.debug(
                    "Skipping SHACL validation result with missing sourceShape/focusNode",
                )
                continue
            if not isinstance(shape_node, URIRef):
                logger.debug(
                    "Skipping SHACL validation result with non-IRI sourceShape",
                )
                continue
            if not isinstance(focus_node, URIRef):
                logger.debug(
                    "Skipping SHACL validation result with non-IRI focusNode",
                )
                continue
            focus_iri = identifier_from_iri(str(focus_node))
            focus_obj = KnowledgeGraph.get_object_from_lookup(focus_iri)
            if focus_obj is None:
                focus_obj = KnowledgeGraph.get_object_from_lookup(str(focus_node))
            if focus_obj is None:
                logger.debug(
                    "Skipping SHACL validation result with unknown focusNode",
                )
                continue
            focus_iri = focus_obj.identifier
            shape_iri = str(shape_node)
            policy = sim.effect_engine.policy.rule_for_shape(shape_iri) if sim.effect_engine.policy else None
            rec = SHACLViolationRecord(
                sim_time=sim.env.now,
                operation_id="VALIDATION_AUDIT",
                origin="SHACL",
                severity=policy.severity if policy else "soft",
                disposition=policy.disposition if policy else "committed",
                shape_iri=shape_iri,
                focus_iri=focus_iri,
                message=_first(report, vr, SH.resultMessage),
                report_graph_ttl=report.serialize(format="turtle"),
            )
            sim.callbacks.emit_violation(rec)

    def _on_attach(self, sim: Simulation) -> None:
        if self.inspect_interval > 0:
            self._alive = True
            self._process = sim.env.process(self._run_every_dt(sim))
        else:
            if not self._subscribed:
                def _on_end(proc: OperationProcess):
                    try:
                        conforms, report, _ = sim.effect_engine.validate_now()
                        if not conforms:
                            self._emit_violation_records(sim, report)
                    except Exception as err:
                        logger.warning(f"ValidationAuditSpawner validate_now failed: {err}")

                sim.callbacks.on_operation_end.append(_on_end)
                self._on_end = _on_end
                self._subscribed = True

    def _detach_safe(self, sim: Simulation) -> None:
        self._alive = False
        if self._subscribed and hasattr(self, "_on_end"):
            try:
                sim.callbacks.on_operation_end.remove(self._on_end)
            except ValueError:
                pass
        self._subscribed = False
        if self._process is not None:
            try:
                if not self._process.triggered:
                    self._process.interrupt("ValidationAuditSpawner detached")
            except Exception:
                pass
            self._process = None

    @property
    def requires_until(self) -> bool:
        return self.inspect_interval > 0

    def _run_every_dt(self, sim: Simulation):
        env = sim.env
        dt = self.inspect_interval
        try:
            while self._alive:
                yield env.timeout(self._sim_time(sim, dt))
                try:
                    conforms, report, _ = sim.effect_engine.validate_now()
                    if not conforms:
                        self._emit_violation_records(sim, report)
                except Exception as err:
                    logger.warning(f"ValidationAuditSpawner validate_now failed: {err}")
        except simpy.Interrupt:
            return
