import os
import urllib.parse
import xml.etree.ElementTree as ET
from io import BytesIO

import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import networkx as nx
import requests
from dash import Dash, Output, Input, html, ctx, callback

from libsyn_tools.chem_schema import ReactionNetwork
from libsyn_tools.utils import json_load, CytoNodeData, drawing_url, CytoNode, CytoEdge, CytoEdgeData, FilePath

"""
visualize the synthetic routes of a library

download example see 
https://dash.plotly.com/cytoscape/images
"""

cyto.load_extra_layouts()
STYLE_SHEET = [
    {
        # global style for edges
        'selector': 'edge',
        'style': {
            'curve-style': 'unbundled-bezier',
            # 'curve-style': 'bezier',
            'taxi-direction': 'vertical',
            'target-arrow-shape': 'triangle',
            'target-arrow-color': 'black',
            'target-arrow-width': '2',
            "opacity": "0.9",
            "line-color": "black",
            "overlay-padding": "3px",
            'label': "data(reaction_index)",
            "font-size": "40px",
            "arrow-scale": 2,
            'text-margin-y': -30,  # https://github.com/cytoscape/cytoscape.js/issues/1374#issuecomment-267486645
            'text-margin-x': -20,
            # "text-valign": "top",
            # "text-halign": "center",
        }
    },
    {
        'selector': 'node',
        'style': {
            # 'background-image-containment': 'over',
            # 'background-width': 'auto',
            # 'background-height': 'auto',
            'shape': 'round-rectangle',
            # 'background-clip': 'none',
            # 'background-fit': 'contain',
            'background-fit': 'contain',
            'background-image': 'data(url)',
            # "border-width": "6px",
            # "border-color": "black",
            # "border-opacity": "1.0",
            "background-color": "white",
            # "content": 'data(label)',
            # "text-outline-color": "#77828C"
        }
    },
    {
        'selector': '.starting_material',
        'style': {
            "background-color": "#cff0fa",
        }
    },
    {
        'selector': '.intermediate_material',
        'style': {
            "background-color": "#f0facf",
        }
    },
    {
        'selector': '.targeting_material',
        'style': {
            "background-color": "#ffaba3",
        }
    },
]


def get_routes_graph(reaction_network: ReactionNetwork):
    """
    routes graph is a nx digraph that contains only intermediates, targets, and starting materials from leaf reactions

    :param reaction_network:
    :return:
    """
    rng = reaction_network.nx_digraph
    reaction_nodes = [r.identifier for r in reaction_network.chemical_reactions]
    intermediate_nodes = [n for n in rng.nodes if
                          rng.in_degree(n) == 1 and rng.out_degree(n) == 1 and n not in reaction_nodes]
    target_nodes = reaction_network.target_smiles

    reaction_dep = dict()
    for rn1 in reaction_nodes:
        reaction_dep[rn1] = []
    for rn1 in reaction_nodes:
        for rn2 in reaction_nodes:
            if rn1 == rn2:
                continue
            if nx.has_path(rng, rn1, rn2):
                reaction_dep[rn2].append(rn1)
    leaf_reactions = [k for k, v in reaction_dep.items() if len(v) == 0]

    leaf_reactions_main_reactants = []

    for leaf_reaction in leaf_reactions:
        lr = reaction_network.entity_dictionary[leaf_reaction]
        selected_reactant = sorted(lr.reactants, key=lambda r: len(r.smiles))[-1]
        leaf_reactions_main_reactants.append(selected_reactant.smiles)

    for reaction in reaction_network.chemical_reactions:
        rng = nx.contracted_nodes(rng, reaction.product_smiles, reaction.identifier, self_loops=False)

    nodes_to_keep = target_nodes + intermediate_nodes + leaf_reactions_main_reactants

    rng = rng.subgraph(nodes_to_keep)
    return rng, target_nodes, intermediate_nodes, leaf_reactions_main_reactants


def get_routes_vis_app(
        runs_folder: FilePath,
        run_name: str
):
    run_folder = os.path.join(runs_folder, run_name)
    assert os.path.isdir(run_folder)
    reaction_network = os.path.join(run_folder, "reaction_network.json")
    reaction_network = ReactionNetwork(**json_load(reaction_network))
    component_cytoscape = cyto.Cytoscape(
        id='cytoscape-image-export',
        wheelSensitivity=0.01,
        layout={
            'name': 'dagre',
            'nodeDimensionsIncludeLabels': True,
            'animate': True,
            'animationDuration': 1000,
            'rankDir': 'LR',
            'align': 'UL',
        },
        # style={'width': '100%', 'height': '100%'},  # browser exploded when setting height
        className="border-primary border w-100",
        responsive=True,
        stylesheet=STYLE_SHEET,
        elements=get_routes_graph_cyto(reaction_network)
    )

    app = Dash(
        name=__name__,
        title="LibSyn",
        external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME],
    )

    app.layout = html.Div([
        html.Div(children=[
            component_cytoscape
        ]),

        html.Div(children=[
            html.Div('Download graph:'),
            html.Button("as svg", id="btn-get-svg"),
        ])
    ])

    return app


def get_routes_graph_cyto(reaction_network: ReactionNetwork):
    routes_rng, target_nodes, intermediate_nodes, leaf_nodes = get_routes_graph(reaction_network)
    nodes = []
    edges = []

    reaction_ids = sorted([cr.identifier for cr in reaction_network.chemical_reactions])
    reaction_indexer = {
        cr: ii for ii, cr in enumerate(reaction_ids)
    }
    # for rid, ii in reaction_indexer.items():
    #     print(ii, rid, reaction_network.entity_dictionary[rid].reaction_smiles)

    product_to_reaction_index = dict()
    reactant_to_reaction_index = dict()
    for reaction in reaction_network.chemical_reactions:
        product_to_reaction_index[reaction.product_smiles] = reaction_indexer[reaction.identifier]
        for r in reaction.reactants:
            reactant_to_reaction_index[r.smiles] = reaction_indexer[reaction.identifier]

    for idx, n in enumerate(sorted(routes_rng.nodes)):
        node_data = CytoNodeData(id=n, label=n, )

        smi_url = drawing_url(smiles=urllib.parse.quote(n), size=80)
        response = requests.get(smi_url)
        svg_content = BytesIO(response.content)
        tree = ET.parse(svg_content)
        width, height = tree.getroot().attrib["width"], tree.getroot().attrib["height"]

        width = int(width[:-2]) * 2
        height = int(height[:-2]) * 2
        width = f"{width}px"
        height = f"{height}px"

        node_data['url'] = smi_url
        node = CytoNode(data=node_data, classes=f"node{idx}")
        STYLE_SHEET.append(
            {
                'selector': f'.node{idx}',
                'style': {
                    "width": width,  # already in px
                    "height": height,
                },
            },
        )
        if n in target_nodes:
            node['classes'] += " targeting_material"
        elif n in intermediate_nodes:
            node['classes'] += " intermediate_material"
        elif n in leaf_nodes:
            node['classes'] += " starting_material"
        else:
            raise RuntimeError(n)
        nodes.append(node)

    for u, v in routes_rng.edges:
        rid = product_to_reaction_index[v]
        edge_data = CytoEdgeData(id=str((u, v)), source=u, target=v, reaction_index=rid)
        edge = CytoEdge(data=edge_data, classes="")
        edges.append(edge)
    return nodes + edges


@callback(
    Output("cytoscape-image-export", "generateImage"),
    Input("btn-get-svg", "n_clicks"), )
def get_image(get_svg_clicks):
    # File type to output of 'svg, 'png', 'jpg', or 'jpeg' (alias of 'jpg')

    # 'store': Stores the image data in 'imageData' !only jpg/png are supported
    # 'download'`: Downloads the image as a file with all data handling
    # 'both'`: Stores image data and downloads image as file.
    action = 'store'
    ftype = "png"

    if ctx.triggered:
        if ctx.triggered_id == "btn-get-svg":
            action = "download"
            ftype = "svg"
    # QA: somehow I have to first init ftype with "png" when dash init, otherwise I got error:
    # Failed to execute 'readAsDataURL' on 'FileReader': parameter 1 is not of type 'Blob'.

    return {'type': ftype, 'action': action}


if __name__ == '__main__':
    RUNS_FOLDER = "../../workplace_opt/RUNS"
    RUNS_FOLDER = os.path.abspath(RUNS_FOLDER)
    # RUN_NAME = "FDA-03-09-0-0"
    # RUN_NAME = "VS-07-01-1-1"
    RUN_NAME = "VS-04-06-1-1"
    app = get_routes_vis_app(
        runs_folder=RUNS_FOLDER,
        run_name=RUN_NAME
    )
    app.run(port=8064, debug=False)
