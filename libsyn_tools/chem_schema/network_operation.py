import math
import random
from collections import defaultdict, Counter

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


def _default_process_time(operation: Operation, multiplier: float = 1.0):
    """
    gives the default estimates by the operation's type and some of its properties

    :param operation:
    :param multiplier: a multiplier to be applied to the estimate
    :return:
    """

    setup_time = 0.6  # min
    solid_dispensing_rate = 4.3  # g/min
    liquid_dispensing_rate = 79  # mL/min
    heating_period = 73  # min
    purification_cost_per_batch = 49  # min
    purification_batch_volume = 100  # mL
    clean_container_cost = 2.3  # min
    concentration_cost = 36  # min
    transfer_container_cost = 0.5  # min

    if operation.type == OperationType.TransferSolid:
        # for solid transfer, time = setup time + amount / dispensing rate
        chemical = Chemical(**operation.annotations['chemical'])
        estimate = setup_time + chemical.mass / solid_dispensing_rate

    elif operation.type == OperationType.TransferLiquid:
        # no setup time needed for liquid
        chemical = Chemical(**operation.annotations['chemical'])
        estimate = chemical.mass / liquid_dispensing_rate

    elif operation.type == OperationType.TransferContainer:
        # a flat time cost
        estimate = transfer_container_cost

    elif operation.type == OperationType.Heating:
        # use default heating time
        estimate = heating_period

    elif operation.type == OperationType.Purification:
        # purification is limited by the batch size
        chemical = Chemical(**operation.annotations['chemical'])
        n_batch = chemical.volume // purification_batch_volume + 1
        estimate = purification_cost_per_batch * n_batch

    elif operation.type == OperationType.CleanContainer:
        # a flat cost for container cleaning
        estimate = clean_container_cost

    elif operation.type == OperationType.Concentration:
        # a flat cost for concentration
        estimate = concentration_cost

    else:
        raise ValueError(f"unsupported operation type: {operation.type}")

    return estimate * multiplier


def default_process_time(operations: list[Operation], functional_modules: list[FunctionalModule], seed: int = 42):
    """ calculate process times for a list of operations """
    random.seed(seed)
    for op in operations:
        op.process_times = dict()
    for fm in functional_modules:
        for op in operations:
            if op.type in fm.can_process:
                pt = _default_process_time(op, multiplier=random.uniform(0.8, 1.2))
            else:
                pt = math.inf
            op.process_times[fm.identifier] = pt


class OperationGraph(Entity):
    reaction: ChemicalReaction = None

    functional_modules: list[FunctionalModule] = []

    make_solutions: dict[str, tuple[Operation, Operation]] = dict()
    """ solid solute identifier -> (transfer solid, transfer liquid) """

    add_solutions: dict[str, Operation] = dict()
    """ solid solute identifier -> transfer liquid """

    add_liquids: dict[str, Operation] = dict()
    """ liquid identifier -> transfer liquid, this excludes the solvent """

    add_solvent: Operation | None = None
    """ if no solid then add solvent first """

    heat: Operation | None = None

    purify: Operation | None = None

    concentrate: Operation | None = None

    lmin: dict[str, dict[str, float]] = dict()

    lmax: dict[str, dict[str, float]] = dict()

    @property
    def operation_dictionary(self) -> dict[str, Operation]:
        operations = []
        for o1, o2 in self.make_solutions.values():
            operations += [o1, o2]
        operations += [*self.add_solutions.values()]
        operations += [*self.add_liquids.values()]
        operations += [self.add_solvent, self.heat, self.purify, self.concentrate]
        return {o.identifier: o for o in operations if o}

    @classmethod
    def from_reaction(cls, reaction: ChemicalReaction, lmin_range=(10, 15), lmax_range=(30, 60)):
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
                # solid transfer to make solution
                transfer_solid = Operation(
                    type=OperationType.TransferSolid,
                    annotations={
                        "chemical": solid.model_dump(),
                        "notes": f"transfer solid {solid.smiles} for making solution for {solid.smiles}"
                    },
                )
                # liquid transfer to make solution, we then apply lmax between this and heating
                transfer_liquid = Operation(
                    type=OperationType.TransferLiquid,
                    annotations={
                        "chemical": solvent_this_fold.model_dump(),
                        "notes": f"transfer solvent {solvent.smiles} for making solution for {solid.smiles}"
                    },
                )
                og.make_solutions[solid.identifier] = (transfer_solid, transfer_liquid)

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
            og.add_solvent = solvent_addition

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
            annotations={"notes": f"heating the reaction vessel to afford product: {reaction.product_smiles}"}
        )

        og.purify = Operation(
            type=OperationType.Purification,
            annotations={"chemical": reaction.products[0].model_dump(),
                         "notes": f"purify the product: {reaction.product_smiles}"}
        )

        og.concentrate = Operation(
            type=OperationType.Concentration,
            annotations={"notes": f"concentrate the raw product: {reaction.product_smiles}"}
        )

        # add precedence relationships
        og.purify.precedents.append(og.concentrate.identifier)  # concentrate before purification
        og.concentrate.precedents.append(og.heat.identifier)  # heating before concentration
        og.lmin[og.heat.identifier][og.concentrate.identifier] = random.uniform(*lmin_range)  # minimum lag for cooling
        # additions before heating
        for addition in [*og.add_liquids.values()] + [*og.add_solutions.values()] + [og.add_solvent, ]:
            if addition:
                og.heat.precedents.append(addition.identifier)
        # solvent/solutions addition before liquids
        for add_liquid in [*og.add_liquids.values()]:
            if og.add_solvent:
                add_liquid.precedents.append(og.add_solvent.identifier)
            else:
                for add_solution in [*og.add_solutions.values()]:
                    add_liquid.precedents.append(add_solution.identifier)
        # first solid then liquid in making solutions
        for solid_id, (transfer_solid, transfer_liquid) in og.make_solutions.items():
            transfer_liquid.precedents.append(transfer_solid.identifier)
            # need to first make the solution then add to reaction vessel
            add_solution = og.add_solutions[solid_id]
            add_solution.precedents.append(transfer_liquid.identifier)
            # once the solution is made there is a lmax to the actual heating
            og.lmax[transfer_liquid.identifier][og.heat.identifier] = random.uniform(*lmax_range)
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
    def from_reaction_network(cls, reaction_network: ReactionNetwork):

        random.seed(42)  # a random.uniform maybe used in `OperationGraph.from_reaction`
        operation_graphs = [
            OperationGraph.from_reaction(r) for r in reaction_network.chemical_reactions
        ]

        # connecting operation graphs: if a starting operation X of a reaction uses SMILES1,
        # then the ending operation Y of the reaction that produces SMILES1 precedes X.
        smiles_to_ending_operation = {
            Chemical(**og.ending_operation.annotations['chemical']).smiles: og.ending_operation
            for og in operation_graphs
        }
        for og in operation_graphs:
            for starting_operation_x in og.starting_operations:
                if 'chemical' in starting_operation_x.annotations:
                    need_smiles = Chemical(**starting_operation_x.annotations['chemical']).smiles
                    if need_smiles in smiles_to_ending_operation:
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
