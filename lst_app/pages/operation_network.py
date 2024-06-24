import glob
import math
import os

import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash import html, get_app, Input, Output, no_update, dcc, State
from dash import register_page

from libsyn_tools.chem_schema import ReactionNetwork, Chemical, OperationNetwork, Operation
from libsyn_tools.utils import json_load, FilePath, CytoEdgeData, CytoNodeData, StateOfMatter

"""
operation network page, very similar to the reaction network page
"""
# TODO DRY against network page
# page setup and constants
register_page(__name__, path='/operation_network', description="Operation")
app = get_app()
cyto.load_extra_layouts()
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PAGE_ID_HEADER = "OperationNetwork__"  # ID header for components in this page

# selector for available routes
component_selector_available_routes_id = PAGE_ID_HEADER + "component_selector_available_routes"
routes_folders = sorted(glob.glob(f"{THIS_DIR}/../routes/*"))
routes_folders = [p for p in routes_folders if os.path.isfile(os.path.join(p, "reaction_network.json"))]
component_selector_available_routes = dbc.Select(
    id=component_selector_available_routes_id,
    options=[{"label": os.path.basename(p), "value": p} for p in routes_folders],
    placeholder="Select a reaction network for operations...",
)

# local store for the network
component_store_network_id = PAGE_ID_HEADER + "component_store_network"
component_store_network = dcc.Store(id=component_store_network_id, storage_type="session")


# TODO decide if we really need session storage


# callback for network store
@app.callback(
    Output(component_store_network_id, 'data'),
    Input(component_selector_available_routes_id, 'value')
)
def load_reaction_network_to_operation_network(network_folder: FilePath) -> dict:
    if not network_folder:
        return no_update
    network_path = os.path.join(network_folder, "reaction_network.json")
    reaction_network = ReactionNetwork(**json_load(network_path))
    try:
        operation_network = OperationNetwork.from_reaction_network(reaction_network)
    except:  # TODO better capture errors from constructor
        # liquid fix
        for r in reaction_network.chemical_reactions:
            if all(c.state_of_matter == StateOfMatter.SOLID for c in r.reactants + r.reagents):
                r.reactants[0].state_of_matter = StateOfMatter.LIQUID
        operation_network = OperationNetwork.from_reaction_network(reaction_network)
    return operation_network.model_dump()


# cytoscape canvas for the network
cyto_stylesheet = OperationNetwork.get_cyto_stylesheet()
component_cytoscape_id = PAGE_ID_HEADER + "component_cytoscape"
component_cytoscape = cyto.Cytoscape(
    id=component_cytoscape_id,
    layout={
        'name': 'dagre',
        'nodeDimensionsIncludeLabels': True,
        'animate': True,
        'animationDuration': 1000,
        # 'rankDir': 'LR',
        'align': 'UL',
    },
    # style={'width': '100%', 'height': '100%'},  # browser exploded when setting height
    className="border-primary border w-100",
    responsive=True,
    stylesheet=cyto_stylesheet,
)


# callback for visualize the network in cytoscape
@app.callback(
    Output(component_cytoscape_id, 'elements'),
    Input(component_store_network_id, 'data'),
)
def visualize_operation_network(network_data: dict):
    if not network_data:
        return no_update
    operation_network = OperationNetwork(**network_data)
    elements = operation_network.to_cytoscape_elements(fig_size=80)
    return elements


# information panel
component_information_panel_id = PAGE_ID_HEADER + "component_information_panel"
component_information_panel = html.Div(id=component_information_panel_id)


def _display_network_summary(network: OperationNetwork) -> dbc.Table:
    rows = []
    for k, v in network.summary.items():
        row = html.Tr([html.Td(k), html.Td(str(v))])
        rows.append(row)
    table_body = [html.Tbody(rows)]
    table = dbc.Table(table_body, bordered=True, hover=True, responsive=True, striped=True, color="light")
    return table


def _display_chemical(chemical: Chemical):
    mass = chemical.mass
    if mass is None:
        mass = "None"
    else:
        mass = "{:.3f} g".format(mass)
    moles = chemical.moles
    if moles is None:
        moles = "None"
    else:
        moles = "{:.6f}".format(moles)
    row0 = html.Tr([html.Td("Identifier"), html.Td(chemical.identifier)])
    row1 = html.Tr([html.Td("Mass"), html.Td(mass)])
    row2 = html.Tr([html.Td("Moles"), html.Td(moles)])
    assert len(chemical.is_produced_by) + len(chemical.is_consumed_by) == 1
    if chemical.is_produced_by:
        row3 = html.Tr([html.Td("Produced by"), html.Td(chemical.is_produced_by[0])])
    else:
        assert len(chemical.is_consumed_by) == 1
        row3 = html.Tr([html.Td("Consumed by"), html.Td(chemical.is_consumed_by[0])])
    rows = [row0, row1, row2, row3]
    table_body = [html.Tbody(rows)]
    table = dbc.Table(table_body, bordered=True, hover=True, responsive=True, striped=True, color="light")
    return table


def _display_operation_node(selected_node_data: CytoNodeData):
    operation = selected_node_data['operation']
    operation = Operation(**operation)

    rows = [
        html.Tr([html.Td("ID"), html.Td(operation.identifier)]),
        html.Tr([html.Td("notes"), html.Td(operation.annotations['notes'])]),
        html.Tr([html.Td("from reaction"), html.Td(operation.from_reaction)]),
    ]
    table_body = [html.Tbody(rows)]
    node_table = dbc.Table(table_body, bordered=True, hover=True, responsive=True, striped=True, color="light")
    tables = [node_table, ]

    if 'chemical' in operation.annotations:
        chemical = Chemical(**operation.annotations['chemical'])
        tables.append(_display_chemical(chemical))
    return tables


def _display_precedence_edge(selected_edge_data: CytoEdgeData):
    lmin = selected_edge_data['lmin']
    lmax = selected_edge_data['lmax']
    if lmin is None:
        lmin = math.inf
    if lmax is None:
        lmax = math.inf
    if lmin < 1e16:
        lmin = "{:.3f}".format(lmin)
    if lmax < 1e16:
        lmax = "{:.3f}".format(lmax)
    rows = [
        html.Tr([html.Td("min lag"), html.Td(lmin)]),
        html.Tr([html.Td("max lag"), html.Td(lmax)]),
    ]
    table_body = [html.Tbody(rows)]
    node_table = dbc.Table(table_body, bordered=True, hover=True, responsive=True, striped=True, color="light")
    tables = [node_table, ]
    return tables


# callback for displaying info from selection
@app.callback(
    Output(component_information_panel_id, 'children'),
    Input(component_cytoscape_id, 'selectedNodeData'),
    Input(component_cytoscape_id, 'selectedEdgeData'),
    State(component_store_network, 'data'),
)
def displaying_selected_elements(node_data_list: list[dict], edge_data_list: list[dict], network_data: dict):
    if network_data is None:
        return no_update
    if node_data_list is None:
        node_data_list = []
    if edge_data_list is None:
        edge_data_list = []
    network = OperationNetwork(**network_data)
    if not node_data_list and not edge_data_list:
        return _display_network_summary(network)
    elif len(node_data_list) and not edge_data_list:
        node_data = node_data_list[-1]
        return _display_operation_node(node_data)
    elif len(edge_data_list) and not node_data_list:
        edge_data = edge_data_list[-1]
        return _display_precedence_edge(edge_data)

    # TODO do nothing if more than one type is selected, there should be a better solution
    if len(node_data_list) and len(edge_data_list):
        return no_update


# page layout
layout = dbc.Container(
    dbc.Row(
        [
            component_store_network,
            dbc.Col([component_cytoscape], className="col-8"),
            dbc.Col([
                dbc.InputGroup(
                    [
                        dbc.InputGroupText("Reaction Network"),
                        component_selector_available_routes,
                    ]
                ),
                html.Hr(),
                component_information_panel,
            ], className="col-4"),
        ]
    ),
    # style={"width": "calc(100vw - 100px)"}
)
