import glob
import os
import plotly.graph_objects as go
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash import html, get_app, Input, Output, no_update, dcc, State
from dash import register_page

from libsyn_tools.chem_schema import ReactionNetwork, Chemical, ChemicalReaction, OperationNetwork
from libsyn_tools.opt import SchedulerOutput, SchedulerInput
from libsyn_tools.utils import json_load, FilePath, CytoNodeData, is_uuid

# page setup and constants
register_page(__name__, path='/scheduler', description="Scheduler")
app = get_app()
cyto.load_extra_layouts()
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PAGE_ID_HEADER = "Scheduler__"  # ID header for components in this page
TARGET_JSON0 = "scheduler_output.json"
TARGET_JSON1 = "scheduler_input.json"
TARGET_JSON2 = "operation_network.json"
TARGET_JSON3 = "reaction_network.json"
routes_folders = sorted(glob.glob(f"{THIS_DIR}/../routes/*"))
routes_folders = [p for p in routes_folders if os.path.isfile(os.path.join(p, TARGET_JSON0)) and os.path.isfile(os.path.join(p, TARGET_JSON1))]

# selector for available routes
component_selector_available_routes_id = PAGE_ID_HEADER + "component_selector_available_routes"
component_selector_available_routes = dbc.Select(
    id=component_selector_available_routes_id,
    options=[{"label": os.path.basename(p), "value": p} for p in routes_folders],
    placeholder="Select a reaction network...",
)


# TODO decide if we really need session storage

# cytoscape canvas for the network
component_gantt_id = PAGE_ID_HEADER + "gantt"
component_gantt = dcc.Graph(id=component_gantt_id)

# callback for visualize the network in cytoscape
@app.callback(
    Output(component_gantt_id, 'figure'),
    Input(component_selector_available_routes_id, 'value'),
)
def visualize_scheduler_output(folder: FilePath) -> go.Figure:
    if not folder:
        fig = go.Figure(go.Scatter(x=[], y=[]))
        fig.update_layout(template=None)
        fig.update_xaxes(showgrid=False, showticklabels=False, zeroline=False)
        fig.update_yaxes(showgrid=False, showticklabels=False, zeroline=False)
        return fig
    scheduler_output_path = os.path.join(folder, TARGET_JSON0)
    scheduler_input_path = os.path.join(folder, TARGET_JSON1)
    operation_network_path = os.path.join(folder, TARGET_JSON2)
    reaction_network_path = os.path.join(folder, TARGET_JSON3)

    scheduler_output, scheduler_input, operation_network, reaction_network = json_load(scheduler_output_path), json_load(scheduler_input_path), json_load(operation_network_path), json_load(reaction_network_path)

    scheduler_output = SchedulerOutput(**scheduler_output)
    scheduler_input = SchedulerInput(**scheduler_input)
    operation_network = OperationNetwork(**operation_network)
    reaction_network = ReactionNetwork(**reaction_network)

    fmid2fm = {fm.identifier: fm for fm in scheduler_input.functional_modules}

    fig = go.Figure(layout=go.Layout(margin={"t": 0}))
    for rid, operations in operation_network.operations_by_reaction.items():
        y = []
        x = []
        base = []
        for operation in operations:
            fmid = scheduler_output.assignments[operation.identifier]
            fm = fmid2fm[fmid]
            y.append(fm.name)
            s = scheduler_output.start_times[operation.identifier]
            e = scheduler_output.end_times[operation.identifier]
            base.append(s)
            x.append(e - s)
        bar = go.Bar(base=base, x=x, y=y, orientation='h', name=rid)
        fig.add_trace(bar)
    return fig


# information panel
component_information_panel_id = PAGE_ID_HEADER + "component_information_panel"
component_information_panel = html.Div(id=component_information_panel_id)


# page layout
layout = dbc.Container(
    dbc.Row(
        [
            dbc.Col([
                dbc.InputGroup(
                    [
                        dbc.InputGroupText("Scheduler Instance"),
                        component_selector_available_routes,
                    ]
                ),
                html.Hr(),
                component_information_panel,
            ], className="col-12"),
            dbc.Col([component_gantt], className="col-12"),
        ]
    ),
    # style={"width": "calc(100vw - 100px)"}
)
