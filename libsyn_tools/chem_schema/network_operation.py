from __future__ import annotations

import math
import random
from collections import defaultdict, Counter
import numpy as np
import networkx as nx
from loguru import logger
from networkx.readwrite import json_graph
from plotly.colors import DEFAULT_PLOTLY_COLORS, unlabel_rgb

from .base import Entity
from .chemical import StateOfMatter, Chemical
from .network_reaction import ReactionNetwork
from .operation import Operation, OperationType, FunctionalModule
from .reaction import ChemicalReaction, ChemicalReactionSpecificationLevel
from ..utils import CytoNode, CytoEdge, CytoEdgeData, CytoNodeData, drawing_url


def _default_process_time(operation: Operation, rng: random.Random, multiplier: float = 1.0, ):
    """
    gives the default estimates by the operation's type and some of its properties

    :param operation:
    :param multiplier: a multiplier to be applied to the estimate
    :return:
    """

    setup_time = 0.6  # min
    solid_dispensing_rate = 15  # g/min
    liquid_dispensing_rate = 43  # mL/min
    heating_period_min = 60  # min
    heating_period_max = 360  # min
    purification_cost_per_batch = 49  # min
    purification_batch_volume = 100  # mL
    clean_container_cost = 2.3  # min
    concentration_cost_min = 16  # min
    concentration_cost_max = 60  # min
    transfer_container_cost = 0.5  # min
    dissolution_min = 0  # min
    dissolution_max = 15  # min

    if operation.type == OperationType.TransferSolid:
        # for solid transfer, time = setup time + amount / dispensing rate
        chemical = Chemical(**operation.annotations['chemical'])
        estimate = setup_time + chemical.mass / solid_dispensing_rate

    elif operation.type == OperationType.TransferLiquid:
        # no setup time needed for liquid
        chemical = Chemical(**operation.annotations['chemical'])
        estimate = chemical.volume / liquid_dispensing_rate

    elif operation.type == OperationType.TransferContainer:
        # a flat time cost
        estimate = transfer_container_cost

    elif operation.type == OperationType.Heating:
        # a random heating time
        estimate = rng.uniform(heating_period_min, heating_period_max)

    elif operation.type == OperationType.Purification:
        # purification is limited by the batch size
        chemical = Chemical(**operation.annotations['chemical'])
        n_batch = chemical.volume // purification_batch_volume + 1
        estimate = purification_cost_per_batch * n_batch

    elif operation.type == OperationType.CleanContainer:
        # a flat cost for container cleaning
        estimate = clean_container_cost

    elif operation.type == OperationType.Concentration:
        # a random cost for concentration
        estimate = rng.uniform(concentration_cost_min, concentration_cost_max)

    elif operation.type == OperationType.ConcentrationAndPurification:
        chemical = Chemical(**operation.annotations['chemical'])
        n_batch = chemical.volume // purification_batch_volume + 1
        estimate = purification_cost_per_batch * n_batch
        estimate += rng.uniform(concentration_cost_min, concentration_cost_max)

    elif operation.type == OperationType.MakeSolution:
        # chemical = Chemical(**operation.annotations['solid_chemical'])
        estimate = setup_time  # + chemical.mass / solid_dispensing_rate
        chemical = Chemical(**operation.annotations['liquid_chemical'])
        estimate += chemical.volume / liquid_dispensing_rate
        estimate = rng.uniform(dissolution_min, dissolution_max)
    else:
        raise ValueError(f"unsupported operation type: {operation.type}")
    logger.info(f"{operation.type}: {estimate * multiplier} min")
    return estimate * multiplier


def default_process_time(operations: list[Operation], functional_modules: list[FunctionalModule], rng: random.Random):
    """ calculate process times for a list of operations """
    for op in operations:
        op.process_times = dict()
    for fm in functional_modules:
        for op in operations:
            if op.type in fm.can_process:
                pt = _default_process_time(op, rng=rng, multiplier=rng.uniform(0.7, 1.3))
            else:
                pt = math.inf
            op.process_times[fm.identifier] = pt
    for op in operations:
        logger.info(
            f"{op.type}: processing time mean {np.mean(list(op.finite_process_times.values()))} std {np.std(list(op.finite_process_times.values()))}")


class OperationGraph(Entity):
    reaction: ChemicalReaction = None

    functional_modules: list[FunctionalModule] = []

    make_solutions: dict[str, Operation] = dict()
    """ solid solute xor solvent identifier -> (make solution) """

    add_solutions: dict[str, Operation] = dict()
    """ solid solute xor solvent identifier -> transfer liquid """

    add_liquids: dict[str, Operation] = dict()
    """ liquid identifier -> transfer liquid, this excludes the solvent """

    heat: Operation | None = None

    concentrate_and_purify: Operation | None = None

    lmin: dict[str, dict[str, float]] = dict()

    lmax: dict[str, dict[str, float]] = dict()

    @property
    def operation_dictionary(self) -> dict[str, Operation]:
        operations = []
        operations += [*self.make_solutions.values()]
        operations += [*self.add_solutions.values()]
        operations += [*self.add_liquids.values()]
        operations += [self.heat, self.concentrate_and_purify]
        return {o.identifier: o for o in operations if o}

    @classmethod
    def from_reaction(cls, reaction: ChemicalReaction, rng: random.Random, lmin_range=(10, 60), lmax_range=(30, 90)):
        # TODO this is still a rather ad-hoc template,
        #  users should be able to define their own template via a better interface
        """
        construct the operation network for one reaction

        1. liquid additions are preceded by solid additions so the balance reading does not fluctuate
        2. reaction vessel addition sequence:
            a. add all solutions in arbitrary order, if no solutions then add solvent
            b. add all remaining liquids in arbitrary order
        3. solids have to be added after dissolving in the solvent, these are the solutions mentioned in 2.a.
            a. the solvent of a reaction is the reactant/reagent with the largest mole amount
            b. the solvent is split into n folds, the mass of each fold is proportional to the solute mass
            c. max lag time applies between the making of these solutions and the succeeding heating
        4. lmin is applied between heating and concentration and is proportional to the heating temperature
            (simulating the cooling process)

        :param reaction: the source reaction
        :param functional_modules: all functional modules in the lab
        :return:
        """
        assert reaction.specification_level == ChemicalReactionSpecificationLevel.QUANTIFIED

        og = cls()
        og.lmin = defaultdict(dict)
        og.lmax = defaultdict(dict)

        # get solids and liquids
        solids, liquids = [], []
        for chemical in reaction.reactants + reaction.reagents:
            if chemical.state_of_matter == StateOfMatter.SOLID:
                solids.append(chemical)
            elif chemical.state_of_matter == StateOfMatter.LIQUID:
                liquids.append(chemical)
            else:
                raise ValueError(f"Unsupported state of matter: {chemical}")
        solids = sorted(solids, key=lambda x: x.mass)
        liquids = sorted(liquids, key=lambda x: x.volume)

        assert len(liquids), f"no liquid found for this reaction: {reaction.reaction_smiles}"
        solvent = liquids[-1]

        if len(solids):

            masses = [s.mass for s in solids]
            mass_ratios = [m / sum(masses) for m in masses]
            for solid, ratio in zip(solids, mass_ratios):
                solvent_this_fold = solvent * ratio
                # make solution, we then apply lmax between this and heating
                make_solution = Operation(
                    type=OperationType.MakeSolution,
                    annotations={
                        "solid_chemical": solid.model_dump(),
                        "liquid_chemical": solvent_this_fold.model_dump(),
                        "notes": f"making {solvent.smiles} solution for {solid.smiles}"
                    },

                )
                og.make_solutions[solid.identifier] = make_solution

                # solution addition to reaction vessel
                mix = solvent_this_fold + solid
                solution_addition = Operation(
                    type=OperationType.TransferLiquid,
                    from_reaction=reaction.identifier,
                    annotations={
                        "chemical": mix.model_dump(),
                        "notes": f"add solution {solid.smiles} in {solvent.smiles}"
                    },
                )
                og.add_solutions[solid.identifier] = solution_addition

        else:
            logger.warning(f"no solid in this reaction: {reaction}")
            solvent_addition = Operation(
                type=OperationType.TransferLiquid,
                annotations={
                    "chemical": solvent.model_dump(),
                    "notes": f"add all solvent {solvent.smiles}"
                },
            )
            og.make_solutions[solvent.identifier] = solvent_addition

        # additions for liquids in reactants/reagents excluding the solvent
        remaining_liquids = liquids[:-1]
        for liquid in remaining_liquids:
            transfer_liquid = Operation(
                type=OperationType.TransferLiquid,
                annotations={
                    "chemical": liquid.model_dump(),
                    "notes": f"add liquid {liquid.smiles}"
                },
            )
            og.add_liquids[liquid.identifier] = transfer_liquid

        og.heat = Operation(
            type=OperationType.Heating,
            annotations={
                "temperature": reaction.temperature,
                "notes": f"heating the reaction vessel to afford product: {reaction.product_smiles}"
            }
        )

        og.concentrate_and_purify = Operation(
            type=OperationType.ConcentrationAndPurification,
            annotations={"chemical": reaction.products[0].model_dump(),
                         "notes": f"purify and concentrate the product: {reaction.product_smiles}"}
        )

        # add precedence relationships
        og.concentrate_and_purify.precedents.append(og.heat.identifier)  # heating before concentration
        og.lmin[og.heat.identifier][og.concentrate_and_purify.identifier] = rng.uniform(
            *lmin_range)  # minimum lag for cooling
        # additions before heating
        for addition in [*og.add_liquids.values()] + [*og.add_solutions.values()] + [*og.make_solutions.values()]:
            if addition:
                og.heat.precedents.append(addition.identifier)

        # # solvent/solutions addition before liquids
        # for add_liquid in [*og.add_liquids.values()]:
        #     for add_solution in [*og.add_solutions.values()]:
        #         add_liquid.precedents.append(add_solution.identifier)

        for solid_id, make_solution in og.make_solutions.items():
            # need to first make the solution then add to reaction vessel
            try:
                add_solution = og.add_solutions[solid_id]
                add_solution.precedents.append(make_solution.identifier)
                # once the solution is made there is a lmax to the actual heating
                og.lmax[make_solution.identifier][og.heat.identifier] = rng.uniform(*lmax_range)
            except KeyError:  # when solid_id is actually solvent_id
                assert len(solids) == 0

        og.reaction = reaction
        for operation in og.operation_dictionary.values():
            operation.from_reaction = reaction.identifier
        return og

    @property
    def nx_digraph(self) -> nx.DiGraph:
        """
        convert to a networkx directed graph, nodes are **identifiers**, edge (u, v) means u precedes v
        """
        g = nx.DiGraph()
        for o in self.operation_dictionary.values():
            for p in o.precedents:
                g.add_edge(p, o.identifier, lmax=math.inf, lmin=0)

        # Adding an edge that already exists updates the edge data. (for multi-graph this is different?)
        # see https://networkx.org/documentation/stable/reference/classes/generated/networkx.Graph.add_edge.html
        for u in self.lmax:
            for v in self.lmax[u]:
                g.add_edge(u, v, lmax=self.lmax[u][v])

        attrs_lmin = dict()
        for u in self.lmin:
            for v in self.lmin[u]:
                g.add_edge(u, v, lmax=self.lmin[u][v])
        return g

    @property
    def starting_operations(self) -> list[Operation]:
        g = self.nx_digraph
        sops = []
        for n in g.nodes:
            if g.out_degree(n) > 0 and g.in_degree(n) == 0:
                sops.append(n)
        return [self.operation_dictionary[oid] for oid in sops]

    @property
    def ending_operation(self) -> Operation:
        g = self.nx_digraph
        sops = []
        for n in g.nodes:
            if g.out_degree(n) == 0 and g.in_degree(n) > 0:
                sops.append(n)
        assert len(sops) == 1
        return self.operation_dictionary[sops[0]]


class OperationNetwork(Entity):
    operation_graphs: list[OperationGraph]

    adjacency_data: dict

    def randomly_assign_temperature(self, rng: random.Random, temperature_range=(50, 300)):
        # when not using ASKCOS
        for op in self.operations:
            if "temperature" in op.annotations and op.type == OperationType.Heating:
                op.annotations['temperature'] = rng.uniform(*temperature_range)

    def get_default_compatability_ad_hoc(self) -> dict[str, dict[str, bool]]:
        """
        compatability by ad hoc rules

        :return: d[oid][oid] -> bool
        """

        def temperature_rule(temp: float) -> str:
            temp_c = temp - 273.15
            if 35 < temp_c < 85:
                return "mild"
            elif 85 <= temp_c < 150:
                return "hot"
            elif 150 <= temp_c:
                return "hot+"
            elif 20 < temp_c <= 35:
                return "rt"
            elif -5 < temp_c <= 20:
                return "cold"
            elif temp_c <= -5:
                return "cold+"

        compatability = defaultdict(dict)
        for o_i in self.operations:
            oid = o_i.identifier
            for o_j in self.operations:
                ojd = o_j.identifier
                if oid == ojd:
                    c = 1
                elif o_i.type != o_j.type:  # this seems unnecessary as it is disallowed by `can_process`
                    c = 0
                elif o_i.type == o_j.type == OperationType.Heating:
                    o_i_temperature = o_i.annotations['temperature']
                    o_j_temperature = o_j.annotations['temperature']
                    if o_i_temperature is None or o_j_temperature is None:
                        c = 1  # by default, we **allow** them to be processed on the same module
                    elif temperature_rule(o_j_temperature) == temperature_rule(o_i_temperature):
                        c = 1
                    else:
                        c = 0
                else:
                    c = 1
                compatability[oid][ojd] = c
        return compatability

    def get_default_compatability(self, temperature_threshold: float) -> dict[str, dict[str, bool]]:
        """
        default compatability -- if two operations can be processed on the same module

        :param temperature_threshold: tolerance for heating to be processed on the same module
        :return: d[oid][oid] -> bool
        """
        compatability = defaultdict(dict)
        for o_i in self.operations:
            oid = o_i.identifier
            for o_j in self.operations:
                ojd = o_j.identifier
                if oid == ojd:
                    c = 1
                elif o_i.type != o_j.type:  # this seems unnecessary as it is disallowed by `can_process`
                    c = 0
                elif o_i.type == o_j.type == OperationType.Heating:
                    o_i_temperature = o_i.annotations['temperature']
                    o_j_temperature = o_j.annotations['temperature']
                    if o_i_temperature is None or o_j_temperature is None:
                        c = 1  # by default, we **allow** them to be processed on the same module
                    elif abs(o_i_temperature - o_j_temperature) < temperature_threshold:
                        c = 1
                    else:
                        c = 0
                else:
                    c = 1
                compatability[oid][ojd] = c
        return compatability

    @property
    def nx_digraph(self) -> nx.DiGraph:
        return json_graph.adjacency_graph(self.adjacency_data)

    @property
    def operations(self) -> list[Operation]:
        operations = []
        for og in self.operation_graphs:
            operations += [*og.operation_dictionary.values()]
        return operations

    @property
    def max_processing_times(self):
        mpt = dict()
        for o in self.operations:
            try:
                mpt[o.identifier + f'--{o.type}'] = max([pt for pt in o.process_times.values() if pt < math.inf])
            except ValueError:
                mpt[o.identifier + f'--o.type'] = math.inf
        return mpt

    @property
    def operations_by_reaction(self) -> dict[str, list[Operation]]:
        d = defaultdict(list)
        for o in self.operations:
            d[o.from_reaction].append(o)
        return d

    @property
    def operation_dictionary(self) -> dict[str, Operation]:
        d = dict()
        for og in self.operation_graphs:
            d.update(og.operation_dictionary)
        return d

    @classmethod
    def from_reaction_network(cls, reaction_network: ReactionNetwork, rng: random.Random):

        operation_graphs = [
            OperationGraph.from_reaction(r, rng) for r in reaction_network.chemical_reactions
        ]

        # connecting operation graphs: if a starting operation X of a reaction uses SMILES1,
        # then the ending operation Y of the reaction that produces SMILES1 precedes X.
        smiles_to_ending_operation = {
            Chemical(**og.ending_operation.annotations['chemical']).smiles: og.ending_operation
            for og in operation_graphs
        }
        for og in operation_graphs:
            for starting_operation_x in og.starting_operations:
                for chemical_key in ['chemical', 'solid_chemical', 'liquid_chemical']:
                    try:
                        need_smiles = Chemical(**starting_operation_x.annotations[chemical_key]).smiles
                    except KeyError:
                        need_smiles = None
                    if need_smiles is not None and need_smiles in smiles_to_ending_operation:
                        ending_operation_y = smiles_to_ending_operation[need_smiles]
                        starting_operation_x.precedents.append(ending_operation_y.identifier)
        g = nx.DiGraph()
        for og in operation_graphs:
            # first add binary precedence
            for op in og.operation_dictionary.values():
                for p in op.precedents:
                    g.add_edge(p, op.identifier, lmin=0, lmax=math.inf)
            # then lmin and lmax
            for u in og.lmin:
                for v in og.lmin[u]:
                    g.add_edge(u, v, lmin=og.lmin[u][v])
            for u in og.lmax:
                for v in og.lmax[u]:
                    g.add_edge(u, v, lmax=og.lmax[u][v])
        adj_data = json_graph.adjacency_data(g)
        return cls(operation_graphs=operation_graphs, adjacency_data=adj_data)

    def to_cytoscape_elements(self, fig_size=80) -> list[CytoEdge | CytoNode]:
        """
        convert to cytoscape elements where nodes are operations and edges are precedence relations
        """
        nodes = []
        edges = []
        g = self.nx_digraph

        for oid in g.nodes:
            operation = self.operation_dictionary[oid]
            node_data = CytoNodeData(id=oid, label=operation.type, )
            if 'chemical' in operation.annotations:
                c_data = operation.annotations['chemical']
                c = Chemical(**c_data)
                node_data['url'] = drawing_url(c.smiles, size=fig_size)
            node_data['operation'] = operation.model_dump()
            node = CytoNode(data=node_data, classes=f"ON__{operation.type.value}")
            nodes.append(node)

        for u, v, d in g.edges(data=True):
            try:
                lmin = d['lmin']
            except KeyError:
                lmin = 0
            try:
                lmax = d['lmax']
            except KeyError:
                lmax = math.inf

            # TODO not sure why this can be None, it may come from reading the nx graph dump
            if lmin is None:
                lmin = 0
            if lmax is None:
                lmax = math.inf

            edge_data = CytoEdgeData(id=str((u, v)), source=u, target=v)
            edge_data['lmin'] = lmin
            edge_data['lmax'] = lmax

            if lmin < 1e-6 and lmax == math.inf:
                cyto_class = "ON__bp_edge"  # binary precedence
            elif lmin >= 1e-6 and lmax == math.inf:
                cyto_class = "ON__lmin_edge"
            elif lmin < 1e-6 and lmax < math.inf:
                cyto_class = "ON__lmax_edge"
            elif lmin >= 1e-6 and lmax < math.inf:
                cyto_class = "ON__lmax_and_lmin_edge"
            else:
                raise ValueError(f"unexpected lmin and lmax: {lmin}, {lmax}")
            edge = CytoEdge(data=edge_data, classes=cyto_class)
            edges.append(edge)
        return nodes + edges

    @property
    def summary(self):
        d = {
            "# reactions": len(self.operations_by_reaction),
            "# operations": len(self.operations),
        }
        types = [o.type for o in self.operations]
        counter = Counter(types)
        for k in sorted(set(types)):
            d[f'# type {k}'] = counter[k]
        return d

    @staticmethod
    def get_cyto_stylesheet():
        sheet = [
            {
                # global style for edges
                'selector': 'edge',
                'style': {
                    'curve-style': 'unbundled-bezier',
                    'taxi-direction': 'vertical',
                    'target-arrow-shape': 'triangle',
                    'target-arrow-color': 'black',
                    "opacity": "0.9",
                    "line-color": "black",
                    # "width": "mapData(weight, 0, 1, 1, 8)",
                    "overlay-padding": "3px"
                }
            },
            {'selector': '.ON__lmin_edge', 'style': {"line-color": "red", }},
            {'selector': '.ON__lmax_edge', 'style': {"line-color": "blue", }},
            {'selector': '.ON__lmax_and_lmin_edge', 'style': {"line-color": "green", }},
            {
                'selector': ':selected',
                'style': {
                    'z-index': 1000,
                    # 'background-color': 'SteelBlue',
                    'border-opacity': "1.0",
                    "border-color": "SteelBlue",
                    'line-color': 'SteelBlue',
                    "border-width": "8px",
                }
            },
        ]

        assert len(DEFAULT_PLOTLY_COLORS) >= len([*OperationType])
        otypes = [o.value for o in OperationType]
        for otype, color in zip(otypes, DEFAULT_PLOTLY_COLORS[: len(otypes)]):
            rgb_tuple = (int(d) for d in unlabel_rgb(color))
            hex_color = "#{0:02x}{1:02x}{2:02x}".format(*rgb_tuple)

            if otype == OperationType.Purification.value:
                style = {
                    'selector': f'.ON__{otype}',
                    'style': {
                        "background-color": "white",
                        'background-image': 'data(url)',
                        'width': 180,
                        'height': 100,
                        'shape': 'circle',
                        'background-fit': 'contain',
                        "border-width": "6px",
                        "border-color": hex_color,
                        "border-opacity": "1.0",
                    }
                }
            else:
                style = {
                    'selector': f'.ON__{otype}',
                    'style': {
                        "background-color": hex_color,
                    }
                }
            sheet.append(style)
        return sheet
