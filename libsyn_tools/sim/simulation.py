from __future__ import annotations

import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import simpy
from loguru import logger
from pandas._typing import FilePath
from pydantic import BaseModel
from rdflib import Graph, ConjunctiveGraph
from tqdm import tqdm

from .effect_engine import (
    EffectEngine,
    KnowledgeGraph,
    EngineMechanicalError,
    ContractViolationError,
    SHACLValidationError,
)
from .overlay import SPPTOverlayProvider, CurrentVolumeOverlayProvider
from .knowledge_graph import LabObject, Has_interrupt_events, identifier_from_iri
from .lifecycle import LifecycleCallbacks
from .operation.operation import Operation, _OpState
from .operation.selector import SelectorCancelled
from .operation.runtime import (
    _needs_runtime_tracking,
    get_object_for_resource,
    get_runtime_context,
    get_runtime_state,
)
from .operation.unitary_edit import AddDataProperty
from .report import RunReport
from .validation import require_singleton_or_error


class OperationEventRecord(BaseModel):
    operation_id: str
    timestamp: float
    event_type: str
    operation_data: dict

    def __repr__(self) -> str:  # pragma: no cover
        return f"OperationEventRecord(operation_id={self.operation_id}, timestamp={self.timestamp:.3f}, event_type={self.event_type})"


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
            callbacks: LifecycleCallbacks,  # NEW
    ):
        self.env = env
        self.operation = operation
        self.operation_registry = operation_registry
        self.dependents = dependents

        self.history_log = history_log
        self.effect_engine = effect_engine
        self.speed_factor = speed_factor
        self.callbacks = callbacks  # NEW

        self.done_event = env.event()

        # Track the in-flight pre_act SimPy process so interrupts can cancel it.
        self._pre_act_process: Optional[simpy.events.Process] = None
        self._pending_interrupt: str | None = None
        self.deferred_start: bool = False

    def sim_time(self, dt: float) -> float:
        return dt * self.speed_factor

    def add_event_log(self, event_type: str, data: dict | None = None) -> None:
        # Timestamp is in simulated time (scaled by simulation_speed_factor).
        # Base time can be derived by dividing by speed_factor when needed.
        self.history_log.append(
            OperationEventRecord(
                operation_id=identifier_from_iri(self.operation.identifier),
                timestamp=self.env.now,
                event_type=event_type,
                operation_data=data or self.operation.model_dump(),
            )
        )

    def request_interrupt(self, reason: str) -> None:
        if self.done_event.triggered:
            return
        if self.simpy_process is None:
            self._pending_interrupt = str(reason)
            self._handle_interrupt(str(reason))
            return
        if not getattr(self.simpy_process, "triggered", True):
            self.simpy_process.interrupt(reason)

    def _safe_cleanup(self) -> None:
        try:
            self.operation.cleanup(self.env)
        except Exception as err:
            logger.warning(
                f"{self.operation.identifier}: cleanup failed after abort/interrupt: {err}"
            )

    def _finish_abort(self, *, reason: str, error: Exception | None = None) -> None:
        self._safe_cleanup()
        data = {"reason": reason}
        if error is not None:
            data["error_type"] = type(error).__name__
            data["error"] = str(error)
        self.add_event_log("OPERATION_ABORT", data)
        if not self.done_event.triggered:
            self.done_event.succeed()
        self.callbacks.emit_operation_end(self)

    def _handle_interrupt(self, reason: str) -> None:
        # If we were interrupted while waiting on pre_act(), cancel the
        # child process so it cannot continue acquiring locks in the
        # background after we have cleaned up.
        if self._pre_act_process is not None and not self._pre_act_process.triggered:
            try:
                self._pre_act_process.interrupt(reason)
                self._pre_act_process.defused = True
            except Exception:
                pass
            self._pre_act_process = None

        edits = []
        reason_txt = f"{self.operation.identifier}:{reason}"
        # record interrupt on participants (only if literal IRIs)
        for participant_name in self.operation.model_fields:
            if not participant_name.startswith("participant_"):
                continue
            iri = getattr(self.operation, participant_name)
            if not isinstance(iri, str):
                continue
            edits.append(
                AddDataProperty(
                    instance_1_iri=iri,
                    property_iri=Has_interrupt_events.predicate_iri,
                    data_value=reason_txt,
                )
            )
        if edits:
            # Interrupt bookkeeping should never crash the simulation.
            # If a policy marks interrupt events as forbidden (aborted),
            # swallow the resulting ContractViolationError so locks are
            # still released and the op can terminate cleanly.
            try:
                locked_iris: list[str] = []
                for req in self.operation.locks:
                    resource = getattr(req, "resource", None)
                    if resource is None:
                        continue
                    try:
                        obj = get_object_for_resource(resource)
                    except Exception:
                        continue
                    if obj is not None:
                        locked_iris.append(obj.identifier)
                if not locked_iris:
                    locked_iris = list(self.operation.resources)
                self.effect_engine.apply(
                    edits,
                    self.env,
                    operation_id=self.operation.identifier,
                    locked_iris=locked_iris,
                )
            except Exception as err:
                try:
                    logger.warning(f"Interrupt bookkeeping failed: {err}")
                except Exception:
                    pass
        # Always release locks / finish the op cleanly.
        self._safe_cleanup()

        self.add_event_log("OPERATION_INTERRUPT", {"reason": str(reason)})
        self.callbacks.emit_operation_interrupt(self, str(reason))
        if not self.done_event.triggered:
            self.done_event.succeed()
        self.callbacks.emit_operation_end(self)

    def run(self):
        if self._pending_interrupt is not None and self.done_event.triggered:
            return
        try:
            yield from self._run_core()
        except simpy.Interrupt as interrupt:
            self._handle_interrupt(str(interrupt.cause))
        except SelectorCancelled as cancel:
            self._handle_interrupt(str(cancel))
        except SHACLValidationError:
            self._safe_cleanup()
            raise
        except (EngineMechanicalError, ContractViolationError) as err:
            self._finish_abort(reason=str(err), error=err)
        except Exception as err:
            self._finish_abort(reason="unexpected exception", error=err)

    def _run_core(self):
        try:
            if self._pending_interrupt is not None:
                raise simpy.Interrupt(self._pending_interrupt)
            # scheduled start gate
            if self.operation.scheduled_start_time is not None:
                try:
                    is_finite = math.isfinite(self.operation.scheduled_start_time)
                except TypeError as exc:
                    raise ValueError(
                        f"{self.operation.identifier}: scheduled_start_time must be finite and >= 0"
                    ) from exc
                if not is_finite or self.operation.scheduled_start_time < 0:
                    raise ValueError(
                        f"{self.operation.identifier}: scheduled_start_time must be finite and >= 0"
                    )
                base_now = self.env.now / self.speed_factor
                delay_base = self.operation.scheduled_start_time - base_now
                if delay_base > 0:
                    yield self.env.timeout(self.sim_time(delay_base))

            # wait for precedents
            precedent_events = [self.operation_registry[pid].done_event for pid in self.operation.required_precedents]
            if precedent_events:
                yield simpy.events.AllOf(self.env, precedent_events)

            # pre-act
            self._pre_act_process = self.operation.pre_act(self.env)
            yield self._pre_act_process
            self._pre_act_process = None
            self.operation._mark_running()

            # START
            self.callbacks.emit_operation_start(self)
            self.add_event_log("OPERATION_START")
            logger.debug(f"[t={self.env.now:.2f}] Start {self.operation.__class__.__name__}={self.operation.identifier}")

            # intrinsic duration
            if self.operation.temporal_cost:
                yield self.env.timeout(self.sim_time(self.operation.temporal_cost))

            # apply edits (mechanical checks + SHACL audit)
            staged = self.effect_engine.prepare(self.operation)
            self.effect_engine.apply(
                staged,
                self.env,
                operation_id=self.operation.identifier,
                locked_iris=self.operation.resources,
            )

            # provenance: remember which ops touched each runtime-tracked object
            # Guard against ANNIHILATE ⇒ resource removed; skip non-present objects
            for iri in self.operation.resources:
                obj = KnowledgeGraph.get_object_from_lookup(iri)
                if _needs_runtime_tracking(obj) and getattr(obj, "is_present", {False}) == {True}:
                    get_runtime_state(obj, self.env).recent_operations.append(self.operation)

            # FINISH
            self.operation.post_act(self.env)
            self.add_event_log("OPERATION_END")
            if not self.done_event.triggered:
                self.done_event.succeed()
            self.callbacks.emit_operation_end(self)
            logger.debug(f"[t={self.env.now:.2f}] Finished {self.operation.identifier}")
        except Exception:
            self._safe_cleanup()
            raise


class Simulation:
    """
    Orchestrates a collection of `Operation` instances inside a SimPy environment.

    simulation_speed_factor scales all simulated durations (operations and spawner timeouts)
    by multiplying base time values to obtain sim-time delays.

    Example (enforced SHACL policy):
        policy = PolicyBundle(per_shape={shape_iri: PolicyRule(severity="hard", disposition="aborted")})
        sim = Simulation(ops, shacl_shapes=shape_graph)
        sim.effect_engine.policy = policy
    """

    def __init__(
            self,
            operations: List[Operation],
            *,
            simulation_speed_factor: float = 1.0,
            random_seed: int | None = None,
            shacl_shapes: str | Path | Graph | None = None,
            shacl_inference: str = "owlrl"
    ):
        self.env = simpy.Environment()
        get_runtime_context(self.env)
        self.rng = random.Random(random_seed)

        self.operations = operations
        try:
            is_finite = math.isfinite(simulation_speed_factor)
        except TypeError as exc:
            raise ValueError("simulation_speed_factor must be finite and > 0") from exc
        if not is_finite or simulation_speed_factor <= 0:
            raise ValueError("simulation_speed_factor must be finite and > 0")
        self.speed_factor = simulation_speed_factor

        if shacl_shapes is None:
            shapes_graph = None
        elif isinstance(shacl_shapes, (Graph, ConjunctiveGraph)):
            shapes_graph = shacl_shapes
        else:
            shapes_graph = Graph().parse(str(shacl_shapes), format="turtle")

        self.history_log: List[OperationEventRecord] = []

        # NEW: lifecycle callbacks (shared among engine + processes)
        self.callbacks = LifecycleCallbacks(history_logger=self._log_history_event)

        self.effect_engine = EffectEngine(
            shapes_graph=shapes_graph,
            inference=shacl_inference,
            callbacks=self.callbacks,  # NEW
        )
        self.env._libsyn_effect_engine = self.effect_engine
        self._sppt_overlay = SPPTOverlayProvider(self.callbacks)
        self._current_volume_overlay = CurrentVolumeOverlayProvider()
        self.effect_engine.register_overlay_provider(self._sppt_overlay.snapshot)
        self.effect_engine.register_overlay_provider(self._current_volume_overlay.snapshot)

        self.operation_registry: Dict[str, OperationProcess] = {}
        self.dependents: Dict[str, List[str]] = defaultdict(list)
        self._spawners: list[object] = []

        self._validate_precedents(self.operations)
        self._build_dependency_map()
        self._build_resources()
        self._build_processes()

    @staticmethod
    def _ensure_identifier_form(operation_id: str) -> None:
        normalized = identifier_from_iri(operation_id)
        if normalized != operation_id:
            raise ValueError(
                f"Operation identifiers must use identifier form, got {operation_id!r}"
            )

    @staticmethod
    def _normalize_precedent_ids(precedents: list[str]) -> list[str]:
        normalized: list[str] = []
        seen: set[str] = set()
        for pred in precedents:
            pred_id = identifier_from_iri(pred)
            if pred_id in seen:
                continue
            seen.add(pred_id)
            normalized.append(pred_id)
        return normalized

    @staticmethod
    def _validate_precedents(
            operations: List[Operation],
            available_ids: set[str] | None = None,
    ) -> None:
        if available_ids is None:
            available_ids = {op.identifier for op in operations}
        for op in operations:
            Simulation._ensure_identifier_form(op.identifier)
            op.required_precedents = Simulation._normalize_precedent_ids(
                op.required_precedents
            )
        for op in operations:
            for pred in op.required_precedents:
                if pred not in available_ids:
                    raise ValueError(
                        f"Operation {op.identifier!r} requires precedent {pred!r} "
                        "which is not registered."
                    )
        graph: dict[str, list[str]] = {
            op.identifier: [pred for pred in op.required_precedents if pred in available_ids]
            for op in operations
            if op.identifier in available_ids
        }
        visited: set[str] = set()
        visiting: set[str] = set()

        for node in graph:
            if node in visited:
                continue
            stack: list[tuple[str, iter[str]]] = [(node, iter(graph.get(node, [])))]
            path: list[str] = [node]
            visiting.add(node)
            while stack:
                current, neighbors = stack[-1]
                try:
                    neighbor = next(neighbors)
                except StopIteration:
                    stack.pop()
                    visiting.remove(current)
                    visited.add(current)
                    path.pop()
                    continue
                if neighbor in visited:
                    continue
                if neighbor in visiting:
                    if neighbor in path:
                        cycle_start = path.index(neighbor)
                        cycle = path[cycle_start:] + [neighbor]
                    else:
                        cycle = [neighbor, current, neighbor]
                    raise ValueError(
                        f"Precedent cycle detected: {' -> '.join(cycle)}"
                    )
                visiting.add(neighbor)
                stack.append((neighbor, iter(graph.get(neighbor, []))))
                path.append(neighbor)

    def _build_dependency_map(self) -> None:
        for op in self.operations:
            for pred in op.required_precedents:
                self.dependents[pred].append(op.identifier)

    def _build_resources(self) -> None:
        for obj in LabObject.all_instances():
            self.effect_engine._register_if_new(obj, self.env)

    def _build_processes(self) -> None:
        for op in self.operations:
            if op.sim_state is not _OpState.NEW:
                raise ValueError(f"Operation {op.identifier} already used in another Simulation")
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
                callbacks=self.callbacks,  # NEW
            )

    def _log_history_event(self, operation_id: str, timestamp: float, event_type: str, data: dict) -> None:
        self.history_log.append(
            OperationEventRecord(
                operation_id=identifier_from_iri(operation_id),
                timestamp=timestamp,
                event_type=event_type,
                operation_data=data,
            )
        )

    def _register_spawner(self, spawner: object) -> None:
        if spawner not in self._spawners:
            self._spawners.append(spawner)

    def _attached_spawners(self) -> list[object]:
        attached: list[object] = []
        for spawner in self._spawners:
            sim_ref = getattr(spawner, "_sim_ref", None)
            if sim_ref is None:
                continue
            sim = sim_ref()
            if sim is self:
                attached.append(spawner)
        return attached

    def run(self, until: float | None = None) -> None:
        """Start all processes and block until done or until time limit."""
        for proc in self.operation_registry.values():
            if proc.simpy_process is None:
                proc.simpy_process = self.env.process(proc.run())

        bar = None
        def _bar_on_end(proc: OperationProcess):
            if bar is not None:
                bar.update()

        try:
            # progress bar via lifecycle callback (no monkey-patch)
            bar = tqdm(total=len(self.operation_registry), desc="Sim", unit="op")
            self.callbacks.on_operation_end.append(_bar_on_end)

            logger.info("Simulation start")
            if until is None:
                for spawner in self._attached_spawners():
                    requires_until = getattr(spawner, "requires_until", False)
                    if callable(requires_until):
                        requires_until = requires_until()
                    if requires_until:
                        raise ValueError(
                            "Simulation.run(until=None) requires an explicit until when a periodic spawner is attached"
                        )
            self.env.run(until=until)
            logger.info(f"Simulation end @ t = {self.env.now}")
        finally:
            if _bar_on_end in self.callbacks.on_operation_end:
                self.callbacks.on_operation_end.remove(_bar_on_end)
            if bar is not None:
                bar.close()
            for spawner in list(self._attached_spawners()):
                try:
                    if hasattr(spawner, "detach"):
                        spawner.detach()
                    else:
                        spawner._detach_safe(self)
                except Exception:
                    pass

    def export_event_log(self, filename: FilePath) -> None:
        df_log = pd.DataFrame.from_records([r.model_dump() for r in self.history_log])
        df_log.to_csv(filename, index=False)
        logger.info(f"Event log exported → {filename}")

    def _build_instance_history_dataframe(self) -> pd.DataFrame:
        terminal_events = {"OPERATION_END", "OPERATION_ABORT", "OPERATION_INTERRUPT"}
        end_time_index = {
            r.operation_id: r.timestamp
            for r in self.history_log
            if r.event_type in terminal_events
        }
        rows: list[dict] = []
        ctx = get_runtime_context(self.env, create=False)
        for rs in ctx.runtime_cache.values():
            if not _needs_runtime_tracking(rs.obj):
                continue
            pool_type = None
            if getattr(rs.obj, "has_pool_type", set()):
                pool_type = require_singleton_or_error(
                    getattr(rs.obj, "has_pool_type", set()),
                    "has_pool_type",
                    rs.obj.identifier,
                    context="instance history",
                )
            for operation in rs.recent_operations:
                op_id = identifier_from_iri(operation.identifier)
                rows.append(
                    {
                        "instance_iri": rs.obj.identifier,
                        "instance_type": rs.obj.__class__.__name__,
                        "pool_type": pool_type,
                        "operation_id": op_id,
                        "operation_type": operation.__class__.__name__,
                        "sim_timestamp": end_time_index.get(op_id, None),
                    }
                )
        return pd.DataFrame(rows)

    def spawn_operation(
            self,
            op: Operation,
            *,
            precedents: list[str] | None = None,
            start_immediately: bool = True,
    ) -> "OperationProcess":
        if op.sim_state is not _OpState.NEW:
            raise ValueError(f"Operation {op.identifier!r} is already in use")
        self._ensure_identifier_form(op.identifier)
        if op.identifier in self.operation_registry:
            raise ValueError(f"Operation id {op.identifier!r} already exists")

        merged_precedents = list(op.required_precedents)
        if precedents:
            merged_precedents.extend(precedents)
        op.required_precedents = self._normalize_precedent_ids(merged_precedents)

        existing_ops = [proc.operation for proc in self.operation_registry.values()]
        self._validate_precedents(existing_ops + [op])
        for pred in op.required_precedents:
            proc = self.operation_registry.get(pred)
            if proc is None:
                continue
            if proc.deferred_start and not proc.done_event.triggered:
                raise ValueError(
                    f"Operation {op.identifier!r} depends on deferred operation {pred!r}"
                )
        proc = OperationProcess(
            env=self.env,
            operation=op,
            operation_registry=self.operation_registry,
            dependents=self.dependents,
            history_log=self.history_log,
            effect_engine=self.effect_engine,
            speed_factor=self.speed_factor,
            callbacks=self.callbacks,  # NEW
        )
        proc.deferred_start = not start_immediately
        self.operation_registry[op.identifier] = proc
        for pred in op.required_precedents:
            self.dependents[pred].append(op.identifier)

        if start_immediately:
            proc.simpy_process = self.env.process(proc.run())
        return proc

    @classmethod
    def compile_actions(cls, *actions: Operation, **kwargs) -> "Simulation":
        return cls(list(actions), **kwargs)

    def export_instance_history(self, filename: FilePath) -> None:
        df = self._build_instance_history_dataframe()
        df.to_csv(filename, index=False)
        logger.info(f"Instance history exported → {filename}")

    def build_report(self, include_ttl: bool = False) -> RunReport:
        terminal_event_types = {"OPERATION_END", "OPERATION_ABORT", "OPERATION_INTERRUPT"}
        event_log_df = pd.DataFrame.from_records([r.model_dump() for r in self.history_log])
        if event_log_df.empty:
            event_log_df = pd.DataFrame(
                columns=["operation_id", "timestamp", "event_type", "operation_data"]
            )
        effects_by_op = {
            identifier_from_iri(op_id): proc.operation.describe_effects()
            for op_id, proc in self.operation_registry.items()
        }
        if not event_log_df.empty and "operation_id" in event_log_df:
            event_log_df["effect_descriptions"] = event_log_df["operation_id"].map(
                lambda op_id: effects_by_op.get(op_id, [])
            )
        violation_records = self.effect_engine.get_violation_records()
        if include_ttl:
            shacl_df = pd.DataFrame.from_records([rec.model_dump() for rec in violation_records])
        else:
            shacl_df = pd.DataFrame.from_records(
                [rec.model_dump(exclude={"report_graph_ttl"}) for rec in violation_records]
            )
        instance_history_df = self._build_instance_history_dataframe()

        start_times = [r.timestamp for r in self.history_log if r.event_type == "OPERATION_START"]
        terminal_times = [r.timestamp for r in self.history_log if r.event_type in terminal_event_types]
        if terminal_times:
            earliest_start = min(start_times) if start_times else min(terminal_times)
            makespan = max(terminal_times) - earliest_start
        else:
            makespan = None

        op_counts = {
            "start": sum(1 for r in self.history_log if r.event_type == "OPERATION_START"),
            "end": sum(1 for r in self.history_log if r.event_type == "OPERATION_END"),
            "abort": sum(1 for r in self.history_log if r.event_type == "OPERATION_ABORT"),
            "interrupt": sum(1 for r in self.history_log if r.event_type == "OPERATION_INTERRUPT"),
        }
        started_ops = {r.operation_id for r in self.history_log if r.event_type == "OPERATION_START"}
        terminal_ops = {r.operation_id for r in self.history_log if r.event_type in terminal_event_types}
        in_progress_ops = started_ops - terminal_ops
        op_counts["in_progress"] = len(in_progress_ops)

        violation_counts = {
            "by_origin": {},
            "by_severity": {},
            "by_shape_iri": {},
        }
        violation_by_shape_disposition: dict[tuple[str | None, str], int] = {}
        for rec in violation_records:
            violation_counts["by_origin"][rec.origin] = violation_counts["by_origin"].get(rec.origin, 0) + 1
            violation_counts["by_severity"][rec.severity] = violation_counts["by_severity"].get(rec.severity, 0) + 1
            shape_key = rec.shape_iri if rec.shape_iri is not None else "None"
            violation_counts["by_shape_iri"][shape_key] = violation_counts["by_shape_iri"].get(shape_key, 0) + 1
            disp_key = (rec.shape_iri, rec.disposition)
            violation_by_shape_disposition[disp_key] = violation_by_shape_disposition.get(disp_key, 0) + 1

        if started_ops and "operation_id" in instance_history_df.columns:
            utilization_df = instance_history_df[instance_history_df["operation_id"].isin(started_ops)]
        else:
            utilization_df = instance_history_df

        utilization_by_pool_type: dict[str, int] = {}
        for pool_type in utilization_df.get("pool_type", pd.Series([], dtype=object)).dropna():
            utilization_by_pool_type[pool_type] = utilization_by_pool_type.get(pool_type, 0) + 1

        utilization_by_module: dict[str, int] = {}
        if "pool_type" in utilization_df.columns:
            module_rows = utilization_df[utilization_df["pool_type"] == "MODULE"]
            for module_iri in module_rows["instance_iri"]:
                utilization_by_module[module_iri] = utilization_by_module.get(module_iri, 0) + 1

        remediation_ops_spawned = 0
        if not event_log_df.empty and "operation_data" in event_log_df:
            start_rows = event_log_df[event_log_df["event_type"] == "OPERATION_START"]
            for data in start_rows["operation_data"]:
                if isinstance(data, dict) and data.get("remediation"):
                    remediation_ops_spawned += 1

        summary = {
            "makespan": makespan,
            "operation_counts": op_counts,
            "violation_counts": violation_counts,
            "violations_by_shape_disposition": [
                {"shape_iri": shape_iri, "disposition": disposition, "count": count}
                for (shape_iri, disposition), count in violation_by_shape_disposition.items()
            ],
            "resource_utilization": {
                "by_pool_type": utilization_by_pool_type,
                "by_module": utilization_by_module,
            },
            "remediation_ops_spawned": remediation_ops_spawned,
        }

        return RunReport(
            summary=summary,
            event_log=event_log_df,
            shacl_violations=shacl_df,
            instance_history=instance_history_df,
        )

    def run_and_report(
            self,
            out_dir: str | Path,
            until: float | None = None,
            include_ttl: bool = False,
    ) -> RunReport:
        self.run(until=until)
        report = self.build_report(include_ttl=include_ttl)
        report.write_dir(out_dir, include_ttl=include_ttl)
        return report
