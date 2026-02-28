# ### THIS IS THE START OF CONTENT OF examples/sim/parallel_flavor_blending.py ###
#!/usr/bin/env python3
"""
Parallel flavor blending — endogenous spawners, SPPT + environment overlays,
coverage-safe effects, enforced "rinse-on-ingredient-switch" policy, SHAKE participant
augmentation, ambient sampler, and JSON exports for visualization.

New in this version
-------------------
• **Hard policy:** Before any device (needle) dispenses a *different* ingredient than its
  previous one, an explicit **RINSE** is inserted and required to complete first.
  This is implemented by building a **device-level queued plan**: every op that touches a device
  (RINSE or ING) is chained with device-level precedents so the actual runtime order on the device
  matches the plan and the rinse barrier always happens before the switch.
• **Rinse operation:** `RinseNeedleGeneric` **annihilates all residue POMs** on the device.
  (No POM locking required; present flag becomes False and filters respect it.)
• SHAKE participant augmentation overlay and ambient sampler are kept so the exported
  contacts/exposures plots are informative.

Design notes
------------
• We no longer rely on analyte-switch deep-clean; with enforced rinses the switch is clean.
  (You can still add a periodic deep-clean spawner if desired.)
• We removed the threshold-based Idle Rinse spawner to avoid duplicate rinses; rinses are scheduled explicitly.
• We keep device residue 'film' stamping in transfers so rinses act on tangible residues.
• Build-time device assignment is used (round-robin over needles) to make the policy precise and race-free.
"""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Dict, Tuple, List, Optional, Any
from uuid import uuid4

from pydantic import PrivateAttr, Field
from rdflib import Graph, Namespace
from rdflib.namespace import RDF
from rdflib.term import URIRef
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim import Simulation
from libsyn_tools.sim.lifecycle import LifecycleCallbacks
from libsyn_tools.sim.operation.operation import Operation, StrOrSelector
from libsyn_tools.sim.operation.selector import AttributeSelector
from libsyn_tools.sim.operation.unitary_edit import (
    Create, AddObjectProperty, RemoveObjectProperty, Annihilate, UnitaryEdit
)
from libsyn_tools.sim.spawner import Spawner, ProcessInterruptSpawner
from libsyn_tools.sim.overlay.sppt_overlay import SPPTOverlayProvider
from libsyn_tools.sim.knowledge_graph.ontology import (
    LabObject, MaterialContainer, PortionOfMaterial,
    # SPPT predicates
    Has_participant, Has_interval, Has_begin_time, Has_end_time, Occurs_in,
    # domain convenience relations
    Is_immediate_part_of, Is_directly_contained_by,
)
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByVolume
from libsyn_tools.chem_schema import Chemical

LIB = Namespace("https://libsyn-sim/kg/")


# =============================================================================
# Registration & world seeding (Create-free; consistent with engine)
# =============================================================================

def _ensure_object_lookups() -> None:
    """Ensure every LabObject subclass has an `object_lookup` registry."""
    def walk(c):
        yield c
        for sc in c.__subclasses__():
            yield from walk(sc)
    for cls in set(walk(LabObject)):
        if getattr(cls, "object_lookup", None) is None:
            setattr(cls, "object_lookup", {})


def _register_obj(o: LabObject) -> None:
    cls = o.__class__
    d = getattr(cls, "object_lookup", None)
    if d is None:
        setattr(cls, "object_lookup", {})
        d = getattr(cls, "object_lookup")
    d[o.identifier] = o


def _present(o: LabObject) -> None:
    # direct world-present flag (reserve Create/Annihilate for effects)
    o.is_present = {True}  # type: ignore[attr-defined]


def new_object(cls, *, iri: Optional[str] = None, pool: Optional[str] = None) -> LabObject:
    o = cls(identifier=iri or f"{cls.__name__}_{uuid4().hex[:6]}")
    if pool:
        o.has_pool_type.add(pool)
    _register_obj(o)
    _present(o)
    return o


def new_container(*, iri: Optional[str] = None, pool: Optional[str] = None) -> MaterialContainer:
    return new_object(MaterialContainer, iri=iri, pool=pool)  # type: ignore


def seed_ingredient(pool: str, total_ml: float, n: int) -> List[MaterialContainer]:
    """Create n reservoirs, each holding a POM with volume==mass (density=1.0)."""
    rs: List[MaterialContainer] = []
    per = total_ml / n
    for i in range(n):
        r = new_container(iri=f"{pool}_{i+1}", pool=pool)
        pom = PortionOfMaterial(identifier=f"POM_{pool}_{i+1}")
        pom.add_chemical(Chemical(mass=per, density=1.0))
        _register_obj(pom)
        _present(pom)
        AddObjectProperty(
            instance_1_iri=pom.identifier,
            instance_2_iri=r.identifier,
            property_iri=Is_directly_contained_by.predicate_iri,
        ).apply()
        rs.append(r)
    return rs


def build_racks_and_vials() -> Tuple[LabObject, LabObject, List[MaterialContainer]]:
    rack1 = new_object(LabObject, iri="RACK1", pool="RACK")
    rack2 = new_object(LabObject, iri="RACK2", pool="RACK")
    vials: List[MaterialContainer] = []
    for i in range(6):
        v = new_container(iri=f"F{i+1}", pool="FLAVOR_VIAL")
        AddObjectProperty(
            instance_1_iri=v.identifier,
            instance_2_iri=rack1.identifier,
            property_iri=Is_immediate_part_of.predicate_iri,
        ).apply()
        vials.append(v)
    for i in range(6, 12):
        v = new_container(iri=f"F{i+1}", pool="FLAVOR_VIAL")
        AddObjectProperty(
            instance_1_iri=v.identifier,
            instance_2_iri=rack2.identifier,
            property_iri=Is_immediate_part_of.predicate_iri,
        ).apply()
        vials.append(v)
    return rack1, rack2, vials


# =============================================================================
# Operations (coverage-safe)
# =============================================================================

class ShakeRack(Operation):
    """Resource-holding only; interval captured by SPPT."""
    participant_rack: StrOrSelector
    def get_operation_effects(self) -> List[UnitaryEdit]: return []


class ResumeShake(Operation):
    """Resource-holding only; resumes a SHAKE after interruption."""
    participant_rack: StrOrSelector
    def get_operation_effects(self) -> List[UnitaryEdit]: return []


class RinseNeedleGeneric(Operation):
    """
    Enforced rinse: annihilate all residue POMs currently on the device.
    (No POM locking required; they become non-present and are ignored by queries.)
    """
    participant_device: StrOrSelector
    participant_waste: StrOrSelector | None = None  # optional, marker only

    def get_operation_effects(self) -> List[UnitaryEdit]:
        dev = self.participant_device
        edits: List[UnitaryEdit] = []
        poms = LabObject.get_directly_contained_individuals(dev, PortionOfMaterial, only_present=True)
        for pom in poms:
            edits.append(Annihilate(instance_1_iri=pom.identifier))
        # Optional marker to make rinses visible in instance history
        if self.participant_waste:
            mid = f"POM_RINSEMARK_{uuid4().hex[:6]}"
            _register_obj(PortionOfMaterial(identifier=mid))
            edits += [
                Create(instance_1_iri=mid),
                AddObjectProperty(
                    instance_1_iri=mid,
                    instance_2_iri=self.participant_waste,
                    property_iri=Is_directly_contained_by.predicate_iri,
                ),
            ]
        return edits


class SealVial(Operation):
    """Create CAP_<vial> and attach via immediate-part-of (interpreted as “sealed”)."""
    participant_vial: StrOrSelector

    def get_operation_effects(self) -> List[UnitaryEdit]:
        vial = self.participant_vial
        cap_id = f"CAP_{vial}"
        _register_obj(LabObject(identifier=cap_id))
        return [
            Create(instance_1_iri=cap_id),
            AddObjectProperty(
                instance_1_iri=cap_id,
                instance_2_iri=vial,
                property_iri=Is_immediate_part_of.predicate_iri,
            ),
        ]


class TransferMaterialByVolumeWithResidue(TransferMaterialByVolume):
    """
    Stock transfer + small semantic “film” POM created on the device so rinses have
    something concrete to remove. Subject is created here; device is locked → OK.
    """
    residue_fraction: float = 0.002

    def get_operation_effects(self) -> List[UnitaryEdit]:
        edits = super().get_operation_effects()
        dev = getattr(self, "participant_device", None)
        if isinstance(dev, str) and self.residue_fraction > 0:
            film_id = f"POM_FILM_{dev}_{uuid4().hex[:6]}"
            _register_obj(PortionOfMaterial(identifier=film_id))
            edits += [
                Create(instance_1_iri=film_id),
                AddObjectProperty(
                    instance_1_iri=film_id,
                    instance_2_iri=dev,
                    property_iri=Is_directly_contained_by.predicate_iri,
                ),
            ]
        return edits


class AmbientExposure(Operation):
    """Expose a single vial to Ambient for `temporal_cost` seconds (no edits)."""
    participant_vial: StrOrSelector
    def get_operation_effects(self) -> List[UnitaryEdit]: return []


# =============================================================================
# Environment overlay + SPPT helpers + SHAKE participant augmentation
# =============================================================================

class EnvironmentTagProvider:
    """
    Process→occurs_in→Place/* via op-id prefix:
      SHAKE_* → Shaker;  RINSE_* → RinseStation;  else → Ambient.
    """
    def __init__(self, callbacks: LifecycleCallbacks) -> None:
        self._spans: Dict[str, Tuple[float | None, float | None]] = {}
        callbacks.on_operation_start.append(self._on_start)
        callbacks.on_operation_end.append(self._on_end)

    def _on_start(self, proc) -> None:
        self._spans.setdefault(proc.operation.identifier, (float(proc.env.now), None))

    def _on_end(self, proc) -> None:
        op_id = proc.operation.identifier
        t0, _ = self._spans.get(op_id, (None, None))
        self._spans[op_id] = (t0, float(proc.env.now))

    @staticmethod
    def _env_for_op(op_id: str) -> str:
        if op_id.startswith("SHAKE_"): return "Shaker"
        if op_id.startswith("RINSE_"): return "RinseStation"
        return "Ambient"

    def snapshot(self) -> Graph:
        g = Graph()
        for nm in ("Ambient", "Shaker", "RinseStation"):
            g.add((LIB[f"Place/{nm}"], RDF.type, LIB.Place))
        for op_id, (t0, t1) in self._spans.items():
            if t0 is None or t1 is None: continue
            env_iri = LIB[f"Place/{self._env_for_op(op_id)}"]
            g.add((LIB[f"Process/{op_id}"], URIRef(Occurs_in.predicate_iri), env_iri))
        return g


def _interval_of(g: Graph, proc_iri: URIRef) -> Optional[Tuple[float, float]]:
    for _, _, int_iri in g.triples((proc_iri, URIRef(Has_interval.predicate_iri), None)):
        t0 = t1 = None
        for _, _, val in g.triples((int_iri, URIRef(Has_begin_time.predicate_iri), None)):
            t0 = float(val)
        for _, _, val in g.triples((int_iri, URIRef(Has_end_time.predicate_iri), None)):
            t1 = float(val)
        if t0 is not None and t1 is not None:
            return (t0, t1)
    return None


def _occurs_in_of(g: Graph, proc_iri: URIRef) -> Optional[str]:
    for _, _, env in g.triples((proc_iri, URIRef(Occurs_in.predicate_iri), None)):
        return str(env).split("/")[-1]
    return None


def contacts_of(entity_iri: str, overlay: Graph) -> List[Tuple[str, float, float, str]]:
    out: List[Tuple[str, float, float, str]] = []
    x = LIB[entity_iri]
    procs = set(p for p, _, _ in overlay.triples((None, URIRef(Has_participant.predicate_iri), x)))
    for p in procs:
        span = _interval_of(overlay, p)
        if span is None: continue
        t0, t1 = span
        for _, _, y in overlay.triples((p, URIRef(Has_participant.predicate_iri), None)):
            if y == x: continue
            out.append((str(y).split(str(LIB))[-1], t0, t1, str(p).split(str(LIB))[-1]))
    return out


def exposures_of(entity_iri: str, overlay: Graph) -> Dict[str, float]:
    acc: Dict[str, float] = {}
    x = LIB[entity_iri]
    procs = set(p for p, _, _ in overlay.triples((None, URIRef(Has_participant.predicate_iri), x)))
    for p in procs:
        span = _interval_of(overlay, p)
        if span is None: continue
        t0, t1 = span
        env = _occurs_in_of(overlay, p) or "Ambient"
        acc[env] = acc.get(env, 0.0) + (t1 - t0)
    return acc


class AssemblyParticipantAugmenter:
    """
    Adds Has_participant triples for immediate parts of any participant LabObject.
    Effect: SHAKE_<rack> processes now include all vials on that rack.
    """
    def __init__(self, callbacks):
        self._spans: dict[str, dict] = {}
        callbacks.on_operation_start.append(self._on_start)
        callbacks.on_operation_end.append(self._on_end)

    def _on_start(self, proc) -> None:
        op = proc.operation
        self._spans.setdefault(op.identifier, {"t0": float(proc.env.now), "t1": None, "parts": list(op.resources)})

    def _on_end(self, proc) -> None:
        op = proc.operation
        s = self._spans.setdefault(op.identifier, {"t0": None, "t1": None, "parts": list(op.resources)})
        s["t1"] = float(proc.env.now)
        if not s.get("parts") and op.resources:
            s["parts"] = list(op.resources)

    def snapshot(self) -> Graph:
        g = Graph()
        for op_id, span in self._spans.items():
            parts = span.get("parts") or []
            for iri in parts:
                obj = KnowledgeGraph.get_object_from_lookup(iri)
                if not isinstance(obj, LabObject): continue
                for v in LabObject.get_immediate_parts(obj):
                    g.add((LIB[f"Process/{op_id}"], URIRef(Has_participant.predicate_iri), LIB[v.identifier]))
        return g


# =============================================================================
# Spawners (ambient sampler only; rinse is now enforced in the queued plan)
# =============================================================================

class AutoShakeSpawner(Spawner):
    """Spawn SHAKE_<rack> when all vials on that rack have received all recipe keys."""
    vial_to_rack: Dict[str, str]
    rack_to_vials: Dict[str, set[str]]
    recipe_keys: List[str]
    _seen: Dict[str, set[str]] = PrivateAttr(default_factory=dict)
    _shook: set[str] = PrivateAttr(default_factory=set)

    def _on_attach(self, sim: Simulation) -> None:
        self._sim = sim
        self._seen = {v: set() for v in self.vial_to_rack}
        sim.callbacks.on_operation_end.append(self._on_end)

    def _detach_safe(self, sim: Simulation) -> None:
        try: sim.callbacks.on_operation_end.remove(self._on_end)
        except ValueError: pass

    def _on_end(self, proc) -> None:
        opid = proc.operation.identifier
        if not opid.startswith("ING_"): return
        parts = opid.split("_")
        if len(parts) < 2: return
        ing  = "_".join(parts[:-1])
        vial = parts[-1]
        if ing not in self.recipe_keys: return
        self._seen.setdefault(vial, set()).add(ing)
        rack = self.vial_to_rack.get(vial)
        if not rack or rack in self._shook: return
        vials = self.rack_to_vials.get(rack, set())
        done_all = all(self._seen.get(v, set()) >= set(self.recipe_keys) for v in vials)
        if done_all:
            rid = f"SHAKE_{rack}"
            if rid not in self._sim.operation_registry:
                self._sim.spawn_operation(ShakeRack(identifier=rid, participant_rack=rack, temporal_cost=3.0))
                self._shook.add(rack)


class AmbientSamplerSpawner(Spawner):
    """
    After a vial has received all recipe keys, optionally spawn up to N short AmbientExposure
    ops for that vial with randomized durations — diversifies exposure profiles.
    """
    recipe_keys: List[str]
    max_events_per_vial: int = 2
    dur_range: Tuple[float, float] = (0.2, 1.2)
    _seen: Dict[str, set[str]] = PrivateAttr(default_factory=dict)
    _spawned: Dict[str, int] = PrivateAttr(default_factory=dict)

    def _on_attach(self, sim: Simulation) -> None:
        self._sim = sim
        sim.callbacks.on_operation_end.append(self._on_end)

    def _detach_safe(self, sim: Simulation) -> None:
        try: sim.callbacks.on_operation_end.remove(self._on_end)
        except ValueError: pass

    def _on_end(self, proc) -> None:
        opid = proc.operation.identifier
        if not opid.startswith("ING_"): return
        ing, vial = "_".join(opid.split("_")[:-1]), opid.split("_")[-1]
        s = self._seen.setdefault(vial, set()); s.add(ing)
        if s >= set(self.recipe_keys):
            n_already = self._spawned.get(vial, 0)
            if n_already >= self.max_events_per_vial: return
            remain = self.max_events_per_vial - n_already
            k = random.randint(0, remain)
            for _ in range(k):
                dur = random.uniform(*self.dur_range)
                eid = f"AMB_{vial}_{uuid4().hex[:6]}"
                self._sim.spawn_operation(AmbientExposure(
                    identifier=eid, participant_vial=vial, temporal_cost=dur
                ))
                self._spawned[vial] = self._spawned.get(vial, 0) + 1


# =============================================================================
# Interrupt mapping + randomized shake interruption
# =============================================================================

def attach_interrupt_spawner(sim: Simulation) -> None:
    def _factory(proc, reason):
        rack = getattr(proc.operation, "participant_rack", None)
        if isinstance(rack, str):
            return ResumeShake(identifier=f"RESUME_{rack}_{uuid4().hex[:6]}",
                               participant_rack=rack, temporal_cost=0.6)
        return None
    ProcessInterruptSpawner(interrupt_dispatch={"USR": _factory}).attach(sim)


def attach_random_interrupt_for_shake(sim: Simulation, rack_ids: List[str], prob: float = 0.6, seed: int = 9) -> None:
    rng = random.Random(seed)
    targets = [f"SHAKE_{r}" for r in rack_ids]
    chosen = rng.choice(targets) if rng.random() < prob and targets else None
    if not chosen: return

    def _kick():
        for _ in range(6):
            yield sim.env.timeout(0.5)
            proc = sim.operation_registry.get(chosen)
            if proc and proc.simpy_process:
                proc.simpy_process.interrupt("USR")
                break
    sim.env.process(_kick())


# =============================================================================
# "Why am I waiting?" recorder (also exported as waits.json)
# =============================================================================

def attach_wait_recorder(sim: Simulation):
    waits: Dict[str, Dict[str, Any]] = {}
    end_times: Dict[str, float] = {}
    last_device_usage: Dict[str, Tuple[str, float]] = {}

    def _on_end(proc) -> None:
        op = proc.operation
        end_times[op.identifier] = proc.env.now
        dev = getattr(op, "participant_device", None)
        if isinstance(dev, str):
            last_device_usage[dev] = (op.identifier, proc.env.now)

    def _on_start(proc) -> None:
        op = proc.operation
        now = proc.env.now
        preds = op.required_precedents
        earliest = max((end_times.get(p, 0.0) for p in preds), default=0.0)
        waited = max(0.0, now - earliest)

        dev = getattr(op, "participant_device", None)
        prev_holder, prev_holder_end = None, None
        if waited > 0 and isinstance(dev, str) and dev in last_device_usage:
            prev_holder, prev_holder_end = last_device_usage[dev]

        waits[op.identifier] = {
            "waited_s": float(waited),
            "waiting_on_device": dev if isinstance(dev, str) else None,
            "prev_holder": prev_holder,
            "prev_holder_end": float(prev_holder_end) if prev_holder_end is not None else None,
        }

        print(f"[wait] {op.identifier}: waited {waited:.2f}"
              + (f" for {dev} (prev by {prev_holder} @ {prev_holder_end:.2f})"
                 if waited > 0 and isinstance(dev, str) and prev_holder else ""))

    sim.callbacks.on_operation_end.append(_on_end)
    sim.callbacks.on_operation_start.append(_on_start)
    return waits


# =============================================================================
# Build formulation ops (device-queued plan with enforced rinse on ingredient switch)
# =============================================================================

def build_formulation_ops_enforced_rinse(
    *,
    vials: List[MaterialContainer],
    ing_sel: Dict[str, AttributeSelector],
    recipe_ml: Dict[str, float],
    device_ids: List[str],
    leave_residue: bool = True,
) -> List[Operation]:
    """
    Build a plan that **queues operations per device** and enforces a **rinse barrier**
    whenever the device switches to a different ingredient.

    Strategy:
      • Deterministic device assignment (round-robin over needles) per vial.
      • Two kinds of precedents:
          - per-vial chain (dose order inside each vial)
          - per-device chain (FIFO queue on each device)
      • Insert an explicit RINSE op when device's last ingredient != current ingredient.
    """
    OpCls = TransferMaterialByVolumeWithResidue if leave_residue else TransferMaterialByVolume

    assign: Dict[str, str] = {}
    for i, v in enumerate(vials):
        assign[v.identifier] = device_ids[i % max(1, len(device_ids))]

    last_ing_for_dev: Dict[str, Optional[str]] = {dev: None for dev in device_ids}
    last_op_for_dev: Dict[str, Optional[str]] = {dev: None for dev in device_ids}
    last_op_for_vial: Dict[str, Optional[str]] = {v.identifier: None for v in vials}

    ops: List[Operation] = []
    for vial in vials:
        viri = vial.identifier
        dev = assign[viri]
        for ing, vol in recipe_ml.items():
            # If device is switching analyte, insert a RINSE gate first
            if last_ing_for_dev[dev] is not None and last_ing_for_dev[dev] != ing:
                rid = f"RINSE_{dev}_{uuid4().hex[:6]}"
                preds = [x for x in (last_op_for_dev[dev],) if x is not None]
                ops.append(RinseNeedleGeneric(
                    identifier=rid,
                    participant_device=dev,
                    participant_waste="WASTE_RINSE",
                    temporal_cost=0.5,
                    required_precedents=preds.copy(),
                ))
                last_op_for_dev[dev] = rid  # device queue advances

            # Now the actual transfer op
            opid = f"{ing}_{viri}"
            vial_preds = [x for x in (last_op_for_vial[viri],) if x is not None]
            dev_preds  = [x for x in (last_op_for_dev[dev],)  if x is not None]
            all_preds = list(dict.fromkeys(vial_preds + dev_preds))  # unique, keep order

            ops.append(OpCls(
                identifier=opid,
                participant_source=ing_sel[ing],
                participant_destination=viri,
                participant_device=dev,           # literal device → deterministic locking and queueing
                transfer_volume=vol,
                temporal_cost=random.uniform(0.12, 0.30),  # jitter for varied total times
                required_precedents=all_preds,
            ))

            last_op_for_vial[viri] = opid
            last_op_for_dev[dev]   = opid
            last_ing_for_dev[dev]  = ing

    return ops


# =============================================================================
# JSON export helpers
# =============================================================================

def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def export_visualization_jsons(
    out_dir: Path,
    *,
    sim: Simulation,
    overlay: Graph,
    waits: Dict[str, Dict[str, Any]],
    rack_to_vials: Dict[str, set[str]],
    vials: List[MaterialContainer],
    random_seed: Optional[int] = None,
) -> None:
    def _san(v):
        if isinstance(v, (str, int, float, bool)) or v is None: return v
        if isinstance(v, list): return [_san(x) for x in v]
        return getattr(v, "identifier", None) or str(v)

    # 1) Processes (intervals + participants + occurs_in + op_type)
    starts: Dict[str, float] = {}
    ends: Dict[str, float] = {}
    participants: Dict[str, Dict[str, Any]] = {}

    for rec in sim.history_log:
        if rec.event_type == "OPERATION_START":
            starts.setdefault(rec.operation_id, float(rec.timestamp))
            participants.setdefault(rec.operation_id, {k: _san(v) for k, v in rec.operation_data.items()})
        elif rec.event_type == "OPERATION_END":
            ends[rec.operation_id] = float(rec.timestamp)
            if rec.operation_id not in participants:
                participants[rec.operation_id] = {k: _san(v) for k, v in rec.operation_data.items()}

    processes: List[Dict[str, Any]] = []
    for op_id, t0 in starts.items():
        t1 = ends.get(op_id, None)
        proc_iri = LIB[f"Process/{op_id}"]
        env = _occurs_in_of(overlay, proc_iri)
        op_type = op_id.split("_", 1)[0] if "_" in op_id else op_id
        processes.append({
            "op_id": op_id,
            "t0": t0,
            "t1": t1,
            "duration": (t1 - t0) if (t1 is not None) else None,
            "op_type": op_type,                 # for plots
            "occurs_in": env,                   # Shaker/RinseStation/Ambient
            "participants": participants.get(op_id, {}),
        })
    _write_json(out_dir / "processes.json", processes)

    # 2) Waits
    waits_list = [{"op_id": k, **v} for k, v in waits.items()]
    _write_json(out_dir / "waits.json", waits_list)

    # 3) Device timeline (segments per device from ended processes)
    device_segments: Dict[str, List[Dict[str, Any]]] = {}
    for p in processes:
        op_id = p["op_id"]
        t0, t1 = p["t0"], p["t1"]
        dev = p["participants"].get("participant_device")
        if dev and (t0 is not None) and (t1 is not None):
            device_segments.setdefault(dev, []).append({
                "device": dev, "op_id": op_id, "t0": t0, "t1": t1
            })
    _write_json(out_dir / "device_timeline.json", device_segments)

    # 4) Exposures & contacts per vial (via overlay)
    exposures_by_vial: Dict[str, Dict[str, float]] = {}
    contacts_by_vial: Dict[str, List[Dict[str, Any]]] = {}
    for v in vials:
        viri = v.identifier
        exposures_by_vial[viri] = exposures_of(viri, overlay)
        contacts = []
        for other, t0, t1, pid in contacts_of(viri, overlay):
            contacts.append({"other": other, "t0": t0, "t1": t1, "process": pid})
        contacts_by_vial[viri] = contacts

    _write_json(out_dir / "exposures_by_vial.json", exposures_by_vial)
    _write_json(out_dir / "contacts_by_vial.json", contacts_by_vial)

    # 5) Deck layout
    deck = {rack: sorted(list(vs)) for rack, vs in rack_to_vials.items()}
    _write_json(out_dir / "deck_layout.json", deck)

    # 6) Run meta
    meta = {
        "random_seed": random_seed,
        "n_ops": len(sim.operation_registry),
        "n_history_events": len(sim.history_log),
    }
    _write_json(out_dir / "run_meta.json", meta)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    _ensure_object_lookups()

    out_dir = Path(__file__).resolve().parent / "_out_blend_realistic"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Stable Places referenced by environment overlay
    for nm in ("Ambient", "Shaker", "RinseStation"):
        _ = new_object(LabObject, iri=f"Place/{nm}")

    # Waste sinks for materialized effects
    new_container(iri="WASTE_RINSE", pool="WASTE")
    new_container(iri="WASTE_DEEP",  pool="WASTE")

    # Devices (needles) — literal device ids to support policy queueing
    needles = [new_object(LabObject, iri=f"NEEDLE{i+1}", pool="LH") for i in range(2)]
    device_ids = [n.identifier for n in needles]

    # Ingredients
    pools = ["ING_VAN", "ING_CIT", "ING_LIM", "ING_SWE", "ING_TOP"]
    for p in pools:
        seed_ingredient(p, total_ml=200.0, n=3)

    # Racks + vials
    rack1, rack2, vials = build_racks_and_vials()

    # AutoShake maps
    vial_to_rack: Dict[str, str] = {}
    rack_to_vials: Dict[str, set[str]] = {rack1.identifier: set(), rack2.identifier: set()}
    for v in vials[:6]:
        vial_to_rack[v.identifier] = rack1.identifier; rack_to_vials[rack1.identifier].add(v.identifier)
    for v in vials[6:]:
        vial_to_rack[v.identifier] = rack2.identifier; rack_to_vials[rack2.identifier].add(v.identifier)

    # Recipe & selectors
    recipe_ml = {"ING_VAN": 4.0, "ING_CIT": 6.0, "ING_LIM": 5.0, "ING_SWE": 3.0, "ING_TOP": 2.0}
    ing_sel = {p: AttributeSelector(pool_type=p, predicate=lambda o: True) for p in pools}

    # Build dosing ops with enforced rinse-on-ingredient-switch (device-queued plan)
    ops = build_formulation_ops_enforced_rinse(
        vials=vials,
        ing_sel=ing_sel,
        recipe_ml=recipe_ml,
        device_ids=device_ids,
        leave_residue=True,
    )

    # --- Simulation seed (owned by this script) ---
    SEED = 17

    # Simulation + overlays
    sim = Simulation(ops, random_seed=SEED)
    sppt_provider = SPPTOverlayProvider(sim.callbacks)
    env_provider  = EnvironmentTagProvider(sim.callbacks)
    sim.effect_engine.register_overlay_provider(sppt_provider.snapshot)
    sim.effect_engine.register_overlay_provider(env_provider.snapshot)

    # Augment SHAKE participants so vials appear in that process
    assembly_aug = AssemblyParticipantAugmenter(sim.callbacks)
    sim.effect_engine.register_overlay_provider(assembly_aug.snapshot)

    # Spawners
    AutoShakeSpawner(vial_to_rack=vial_to_rack,
                     rack_to_vials=rack_to_vials,
                     recipe_keys=list(recipe_ml.keys())).attach(sim)
    AmbientSamplerSpawner(
        recipe_keys=list(recipe_ml.keys()),
        max_events_per_vial=2,
        dur_range=(0.2, 1.0),
    ).attach(sim)

    # Interrupt mapping + randomized shake interruption
    attach_interrupt_spawner(sim)
    attach_random_interrupt_for_shake(sim, rack_ids=[rack1.identifier, rack2.identifier], prob=0.5, seed=9)

    # Wait recorder
    waits = attach_wait_recorder(sim)

    print("Running realistic flavor blending with enforced rinse-on-switch…")
    sim.run()

    # Overlay snapshot
    overlay = Graph()
    overlay += sppt_provider.snapshot()
    overlay += env_provider.snapshot()
    overlay += assembly_aug.snapshot()

    # Exports
    sim.export_event_log(out_dir / "event_log.csv")
    sim.export_instance_history(out_dir / "instance_history.csv")
    export_visualization_jsons(out_dir, sim=sim, overlay=overlay, waits=waits,
                               rack_to_vials=rack_to_vials, vials=vials,
                               random_seed=SEED)

    print(f"\nExports written to: {out_dir}")
    print("Done.")


if __name__ == "__main__":
    main()
# ### THIS IS THE END OF CONTENT OF examples/sim/parallel_flavor_blending.py ###
