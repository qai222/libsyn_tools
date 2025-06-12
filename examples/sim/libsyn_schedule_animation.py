import dash
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, dcc, State

from libsyn_tools.chem_schema import OperationNetwork, OperationType, ReactionNetwork
from libsyn_tools.opt.schema import Solver
from libsyn_tools.utils import json_load, FilePath
import os

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])


SCHEDULER_INSTANCE_FOLDER  = "../../workplace_opt/RUNS/FDA-03-09-0-0"
SOLVER_FILE = os.path.join(SCHEDULER_INSTANCE_FOLDER, "solver_milp.json")
OPERATION_NETWORK_FILE = os.path.join(SCHEDULER_INSTANCE_FOLDER, "operation_network.json")
SOLVER = Solver(**json_load(SOLVER_FILE))
OPERATION_NETWORK = OperationNetwork(**json_load(OPERATION_NETWORK_FILE))
SIMULATION_DF = pd.read_csv("FDA-03-09-0-0_normal_0.1_0.csv")
SIMULATION_DF_ENDS = SIMULATION_DF[SIMULATION_DF['event_type'] == 'ACTION_END']
SIMULATION_DF_STARTS = SIMULATION_DF[SIMULATION_DF['event_type'] == 'ACTION_START']
ACTION_END_TIMES = dict(zip(SIMULATION_DF_ENDS['action_id'], SIMULATION_DF_ENDS['timestamp']))
ACTION_START_TIMES = dict(zip(SIMULATION_DF_STARTS['action_id'], SIMULATION_DF_STARTS['timestamp']))
SIMULATION_END = max(SIMULATION_DF_ENDS['timestamp'])
SLIDER_STEP = 5


def get_gantt_df():
    records = []
    fms = SOLVER.input.functional_modules
    fms = [f.identifier for f in fms]
    for action_id in ACTION_END_TIMES:
        end_time = ACTION_END_TIMES[action_id]
        start_time = ACTION_START_TIMES[action_id]
        operation = OPERATION_NETWORK.operation_dictionary[action_id]
        if operation.type not in (OperationType.Heating, OperationType.ConcentrationAndPurification):
            continue
        assigned_fm = SOLVER.output.assignments[operation.identifier]
        r = {
            "action_id": action_id,
            "start_time": start_time,
            "end_time": end_time,
            "machine_for_this_action": fms.index(assigned_fm),
        }
        records.append(r)
    return pd.DataFrame.from_records(records)


GANTT_DF = get_gantt_df()


@app.callback(
    [Output("time-slider", "value"),
     Output("interval", "disabled")],
    [Input("start-button", "n_clicks"),
     Input("stop-button", "n_clicks"),
     Input("interval", "n_intervals")],
    State("time-slider", "value")
)
def update_slider(start_clicks, stop_clicks, n_intervals, current_value):
    ctx = dash.callback_context
    if not ctx.triggered:
        # No trigger yet – do nothing.
        raise dash.exceptions.PreventUpdate

    # Identify the triggered input.
    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if trigger_id == "start-button":
        # When the start button is clicked, reset slider and enable the interval.
        return 0, False
    elif trigger_id == "stop-button":
        # When the stop button is clicked, disable the interval (stop animation).
        return current_value, True
    elif trigger_id == "interval":
        # When the interval triggers, update the slider value.
        if current_value >= SIMULATION_END:
            # Stop the animation if the slider reaches 100.
            return current_value, True
        new_value = min(current_value + SLIDER_STEP, SIMULATION_END)
        return new_value, False


def get_sim_cyto_stylesheet():
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
        {'selector': '.FINISHED_OPERATION', 'style': {"background-color": "green", }},
        {'selector': '.ONGOING_OPERATION', 'style': {"background-color": "blue", }},
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
    return sheet


@app.callback(
    Output("COMPONENT_CYTOSCAPE", "elements"),
    Input("time-slider", "value"),
    State("COMPONENT_CYTOSCAPE", "elements")
)
def update_operation_network(current_time, elements):
    # elements = OPERATION_NETWORK.to_cytoscape_elements(fig_size=80)
    finished = [k for k, v in ACTION_END_TIMES.items() if v <= current_time]
    ongoing = [k for k, v in ACTION_START_TIMES.items() if v <= current_time and k not in finished]
    for i in range(len(elements)):
        if elements[i]['data']['id'] in finished:
            elements[i]['classes'] = 'FINISHED_OPERATION'
        elif elements[i]['data']['id'] in ongoing:
            elements[i]['classes'] = 'ONGOING_OPERATION'
    return elements


@app.callback(
    Output("gantt-chart", "figure"),
    Input("time-slider", "value"),
)
def update_gantt(simulation_time):
    df = GANTT_DF
    df_active = df[df["start_time"] <= simulation_time].copy()
    fixed_machines = sorted(df['machine_for_this_action'].unique())

    if df_active.empty:
        fig = go.Figure()
        fig.update_layout(
            # title=f"Simulation Time: {simulation_time}",
            xaxis_title="Time (min)",
            yaxis_title="Functional Module",
        )
        return fig
    # Calculate the progress of each action
    # For actions still in progress, elapsed time = simulation_time - start_time (capped at end_time)
    df_active["progress"] = df_active.apply(
        lambda row: min(simulation_time, row["end_time"]) - row["start_time"],
        axis=1
    )

    # Color actions: green if completed, blue if still ongoing
    df_active["color"] = df_active.apply(
        lambda row: "green" if simulation_time >= row["end_time"] else "blue",
        axis=1
    )

    # Build the bar chart.
    # Each bar starts at the action's start_time (using 'base') and extends by the 'progress'
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=df_active["progress"],
        y=df_active["machine_for_this_action"],
        base=df_active["start_time"],
        orientation='h',
        marker_color=df_active["color"],
        hovertext=[
            f"Action ID: {row['action_id']}<br>Start: {row['start_time']}<br>End: {row['end_time']}"
            for _, row in df_active.iterrows()
        ],
        hoverinfo="text"
    ))

    fig.update_layout(
        title=f"Simulation Time: {simulation_time}",
        xaxis_title="Time (min)",
        yaxis_title="Functional Module",
        yaxis=dict(
            tickmode='array',
            tickvals=fixed_machines,
            ticktext=fixed_machines,
            categoryorder="array",
            categoryarray=fixed_machines
        ),
        barmode="overlay",
        # margin=dict(l=20, r=20, t=40, b=20),
        margin=dict(l=0, r=0, t=0, b=0)
    )
    # No need to reverse y-axis if machines are already in desired order.

    return fig


cyto_stylesheet = get_sim_cyto_stylesheet()
cyto.load_extra_layouts()
component_cytoscape = cyto.Cytoscape(
    id="COMPONENT_CYTOSCAPE",
    wheelSensitivity=0.03,
    elements=OPERATION_NETWORK.to_cytoscape_elements(fig_size=80),
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

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(
            dbc.Button("Start", id="start-button", n_clicks=0, color="primary"),
            width="auto"
        ),
        dbc.Col(
            dbc.Button("Pause", id="stop-button", n_clicks=0, color="danger", className="ml-2"),
            width="auto"
        )
    ], className="my-3"),
    dbc.Row([
        dbc.Col(
            dcc.Slider(
                id="time-slider",
                min=0,
                max=SIMULATION_END,
                step=0.1,
                value=0,
                marks=None,
                tooltip={"always_visible": True, "placement": "bottom"}
            ),
            width=12
        )
    ], className="my-3"),
    dbc.Row(
        [
            dbc.Col(
                dcc.Graph(id="gantt-chart", style={'height': f'600px'}),
                width=7
            ),
            dbc.Col(
                component_cytoscape, width=5
            ),
        ]
    ),

    # Interval component triggers every 1000ms (1 second)
    dcc.Interval(id="interval", interval=150, n_intervals=0, disabled=True)
])

# Run the Dash app.
if __name__ == '__main__':
    app.run_server(debug=True)
