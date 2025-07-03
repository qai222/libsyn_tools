from __future__ import annotations

from collections.abc import Iterable
from typing import List

import simpy
from loguru import logger
from pyshacl import validate
from rdflib import ConjunctiveGraph
from rdflib import Graph, Literal, URIRef, Namespace
from rdflib.namespace import XSD
from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.sim.knowledge_graph.physical_entities import BaseClass, MaterialContainer
from libsyn_tools.sim.operation import Operation, UnitaryEdit, UnitaryEditType, FilterStoreRegistry, get_runtime_state
from libsyn_tools.sim.operation.runtime import _RESOURCE_MAP, _RUNTIME_CACHE, _needs_runtime_tracking

LIB_SYN = Namespace("https://libsyn-sim/kg/")


class SHACLValidationError(RuntimeError):
    """
    Raised when the knowledge graph violates one or more user-supplied
    SHACL shapes *after* an Operation has been applied.
    """
    pass


class EffectEngine:
    """
    This module provides **transaction-like semantics** around
    `UnitaryEdit` application:

    1.  *prepare()* – return a **copy** of the staged edits that an
        `Operation` built in its *pre-act* phase.
    2.  *apply()*  – iterate through the edits **atomically**.
        If any single edit raises, all previously-applied edits are rolled
        back in *reverse* order and the original exception is re-raised.
    """

    def __init__(
            self,
            *,
            shapes_graph: ConjunctiveGraph | None = None,
            raise_shacl: bool = False,
    ):
        """
        Parameters
        ----------
        shapes_graph
            RDF graph containing user-defined SHACL shapes.
        raise_shacl
            if raise when shacl validation fails.
        """
        self.shapes_graph = shapes_graph
        self.inference = "rdfs"
        self.raise_shacl = raise_shacl

    def _build_overlay_graph(self) -> Graph:
        """
        Construct an *ephemeral* rdflib.Graph that contains **derived**
        triples needed for SHACL validation.

        This graph is **not** written back to the KnowledgeGraph.
        """
        g = Graph()

        # add current volume
        for c in MaterialContainer.object_lookup.values():
            if c.is_present != {True}:  # skip annihilated objects
                continue
            vol = c.directly_contained_pom_volume  # existing helper
            g.add(
                (
                    URIRef(c.instance_iri),
                    LIB_SYN.currentVolume,
                    Literal(vol, datatype=XSD.double),
                )
            )
        return g

    def _run_shacl_validation(self) -> None:
        if self.shapes_graph is None:
            return

        # 1) base data graph (asserted triples)
        data_graph: Graph = KnowledgeGraph.graph()

        # 2) overlay graph with derived facts
        overlay_graph: Graph = self._build_overlay_graph()

        # 3) merged view for validation  (ConjunctiveGraph |= overlay)
        union_graph = ConjunctiveGraph()
        for triple in data_graph.triples((None, None, None)):
            union_graph.add(triple)
        for triple in overlay_graph.triples((None, None, None)):
            union_graph.add(triple)

        conforms, report_graph, _ = validate(
            union_graph,
            shacl_graph=self.shapes_graph,
            ont_graph=None,
            inference=self.inference,
            advanced=True,
            debug=False,
        )
        logger.debug(f"shapes conform: {conforms}")
        if not conforms:
            logger.error(report_graph.serialize(format="turtle"))
            if self.raise_shacl:
                raise SHACLValidationError(report_graph)

    def prepare(self, action: Operation) -> List[UnitaryEdit]:
        """Return a **defensive copy** of the staged edits for *action*.

        The `Operation` must already have executed its *pre-act* phase
        so the list is fully populated.
        """
        return action.operation_effects.copy()

    def apply(self, edits: Iterable[UnitaryEdit], env: simpy.Environment) -> None:
        """
        Apply each `UnitaryEdit` **in order**.

        If an edit fails, everything applied so far is *rolled back* and
        the original exception bubbles up.
        """
        applied: list[UnitaryEdit] = []

        for edit in edits:
            subj = KnowledgeGraph.get_object_from_lookup(edit.instance_1_iri)
            obj2 = (
                KnowledgeGraph.get_object_from_lookup(edit.instance_2_iri)
                if edit.type in (UnitaryEditType.ADD_OBJECT_PROPERTY,
                                 UnitaryEditType.REMOVE_OBJECT_PROPERTY)
                   and edit.instance_2_iri
                else None
            )

            self._register_if_new(subj, env)
            # ensure every referenced LabObject owns a Resource
            if edit.type is UnitaryEditType.CREATE:
                pass  # postpone until after .apply()
            elif obj2:  # only for ADD, not REMOVE
                self._register_if_new(obj2, env)

            logger.debug(f"Applying edit: {edit.type} – {subj.__class__.__name__}={edit.instance_1_iri}")
            edit.apply()
            applied.append(edit)

            if edit.type is UnitaryEditType.CREATE:
                self._register_if_new(subj, env)

            self._log_and_sync(subj, edit, env)
            self._log_and_sync(obj2, edit, env)

            if edit.type is UnitaryEditType.ANNIHILATE:
                # remove resource + filter entry only after logging
                self._unregister_object(subj)

        self._run_shacl_validation()

    def _log_and_sync(self, obj: BaseClass, edit: UnitaryEdit,
                      env: simpy.Environment) -> None:
        """Append *edit* to obj.runtime history and normalise FilterStore."""
        if obj is None or not _needs_runtime_tracking(obj):
            return
        get_runtime_state(obj, env).recent_edits.append(edit)
        self._sync_filter_stores(obj, env)

    @staticmethod
    def _register_if_new(obj: BaseClass, env: simpy.Environment):
        if not _needs_runtime_tracking(obj):
            return
        if obj.instance_iri not in _RESOURCE_MAP:
            _RESOURCE_MAP[obj.instance_iri] = simpy.Resource(env, capacity=1)
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)
            logger.debug(f"auto register new object: {obj.__class__.__name__}={obj.instance_iri}")

    def _unregister_object(self, obj: BaseClass):
        """Remove *all* runtime artefacts for a vanished object."""
        if not _needs_runtime_tracking(obj):
            return
        _RESOURCE_MAP.pop(obj.instance_iri, None)
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        _RUNTIME_CACHE.pop(obj.instance_iri, None)
        setattr(obj, "_runtime", None)  # clear pointer

    def _sync_filter_stores(self, obj: BaseClass, env: simpy.Environment):
        """
        After *any* edit we normalise FilterStore membership so rollbacks
        never leave an object in the wrong pool.
        """
        if not _needs_runtime_tracking(obj):
            return
        FilterStoreRegistry.remove_obj_from_filter_store(obj)
        if obj.is_present == {True}:  # ← guard
            FilterStoreRegistry.put_obj_into_filter_store(obj, env)
