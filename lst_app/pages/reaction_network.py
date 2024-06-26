import glob
import os

import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash import html, get_app, Input, Output, no_update, dcc, State
from dash import register_page

from libsyn_tools.chem_schema import ReactionNetwork, Chemical, ChemicalReaction
from libsyn_tools.utils import json_load, FilePath, CytoNodeData, is_uuid

"""
reaction network page
- visualize a selected reaction network

# TODO
- quantify the network based on user input
- export/load reaction network, either quantified or unquantified
"""
# page setup and constants
register_page(__name__, path='/reaction_network', description="Reaction")
app = get_app()
cyto.load_extra_layouts()
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PAGE_ID_HEADER = "ReactionNetwork__"  # ID header for components in this page
TARGET_JSON = "reaction_network.json"
routes_folders = sorted(glob.glob(f"{THIS_DIR}/../routes/*"))
routes_folders = [p for p in routes_folders if os.path.isfile(os.path.join(p, TARGET_JSON))]

# selector for available routes
component_selector_available_routes_id = PAGE_ID_HEADER + "component_selector_available_routes"
component_selector_available_routes = dbc.Select(
    id=component_selector_available_routes_id,
    options=[{"label": os.path.basename(p), "value": p} for p in routes_folders],
    placeholder="Select a reaction network...",
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
def load_reaction_network(network_folder: FilePath) -> dict:
    if not network_folder:
        return no_update
    network_path = os.path.join(network_folder, TARGET_JSON)
    return json_load(network_path)


# cytoscape style sheet, the classes are defined in
# `libsyn_tools.chem_schema.network.ReactionNetwork.to_cytoscape_elements`
# TODO too much code here, store it as a separate file
style_sheet_reaction_network = [
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
    {
        'selector': '.RN__ChemicalReaction',
        'style': {
            'width': 800,
            'height': 300,
            'shape': 'rectangle',
            'background-fit': 'contain',
            'background-image': 'data(url)',
            "border-width": "6px",
            "border-color": "black",
            "border-opacity": "1.0",
            "background-color": "white",
        }
    },
    {
        'selector': '.RN__MolecularSMILES',
        'style': {
            'width': 180,
            'height': 100,
            'shape': 'circle',
            'background-fit': 'contain',
            'background-image': 'data(url)',
            "border-width": "6px",
            "border-color": "black",
            "border-opacity": "1.0",
            "background-color": "white",
            # "content": 'data(label)',
            # "text-outline-color": "#77828C"
        }
    },
    {
        'selector': '.RN__MolecularSMILES_starting_material',
        'style': {
            "background-color": "#cff0fa",
        }
    },
    {
        'selector': '.RN__MolecularSMILES_intermediate_material',
        'style': {
            "background-color": "#f0facf",
        }
    },
    {
        'selector': '.RN__MolecularSMILES_targeting_material',
        'style': {
            "background-color": "#ffaba3",
        }
    },
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

# cytoscape canvas for the network
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
    stylesheet=style_sheet_reaction_network,
)


# callback for visualize the network in cytoscape
@app.callback(
    Output(component_cytoscape_id, 'elements'),
    Input(component_store_network_id, 'data'),
)
def visualize_reaction_network(network_data: dict):
    if not network_data:
        return no_update
    network = ReactionNetwork(**network_data)
    elements = network.to_cytoscape_elements(fig_size=80)
    return elements


# information panel
component_information_panel_id = PAGE_ID_HEADER + "component_information_panel"
component_information_panel = html.Div(id=component_information_panel_id)


def _display_network_summary(network: ReactionNetwork) -> dbc.Table:
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


def _display_reaction_node(selected_node_data: CytoNodeData):
    reaction = ChemicalReaction(**selected_node_data['chemcail_reaction'])
    rows = [
        html.Tr([html.Td("Identifier"), html.Td(reaction.identifier)]),
        html.Tr([html.Td("SMILES"), html.Td(reaction.reaction_smiles)]),
        html.Tr([html.Td("Temperature (C)"), html.Td(reaction.temperature)]),
    ]
    if reaction.expected_yields is not None:
        for product_id, expected_yield in reaction.expected_yields.items():
            row = html.Tr(
                [html.Td("Yield of " + reaction.chemical_dictionary[product_id].smiles), html.Td(expected_yield)])
            rows.append(row)
    if reaction.intended_ratios is not None:
        for c_id, ratio in reaction.intended_ratios.items():
            row = html.Tr([html.Td("Intended ratio of " + reaction.chemical_dictionary[c_id].smiles), html.Td(ratio)])
            rows.append(row)

    table_body = [html.Tbody(rows)]
    node_table = dbc.Table(table_body, bordered=True, hover=True, responsive=True, striped=True, color="light")
    return node_table


def _display_molecular_node(selected_node_data: CytoNodeData):
    list_of_chemicals = selected_node_data['list_of_chemicals']
    list_of_chemicals = [Chemical(**d) for d in list_of_chemicals]

    # TODO assuming here all chemicals share the same properties except the quantities
    example_chemical = list_of_chemicals[0]
    rows = [
        html.Tr([html.Td("SMILES"), html.Td(example_chemical.smiles)]),
        html.Tr([html.Td("state of matter"), html.Td(example_chemical.state_of_matter)]),
        html.Tr([html.Td("density"), html.Td(example_chemical.density)]),
    ]
    table_body = [html.Tbody(rows)]
    node_table = dbc.Table(table_body, bordered=True, hover=True, responsive=True, striped=True, color="light")
    tables = [node_table, ]

    for c in list_of_chemicals:
        t = _display_chemical(c)
        tables.append(t)
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
    network = ReactionNetwork(**network_data)
    if not node_data_list and not edge_data_list:
        return _display_network_summary(network)
    if len(node_data_list) and not edge_data_list:
        node_data = node_data_list[-1]
        if not is_uuid(node_data['id']):
            return _display_molecular_node(node_data)
        else:
            return _display_reaction_node(node_data)
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
