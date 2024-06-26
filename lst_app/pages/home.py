import dash_bootstrap_components as dbc
from dash import html, get_app
from dash import register_page
"""
home page for the app
"""

register_page(__name__, path='/', description="Home")
app = get_app()

def get_card(title: str, text: str, link_text: str, link_path: str = "#", width: int = 22):
    """ subpage card """
    card = dbc.Card(
        [
            dbc.CardBody(
                [
                    html.H4(title, className="card-title"),
                    html.P(
                        text,
                        className="card-text",
                    ),
                    dbc.Button(link_text, color="primary", href=link_path),
                ]
            ),
        ],
        style={"min-width": f"{width}rem", "width": f"{width}rem"},
        className="mx-4"
    )
    return card


subpage_cards = []

card_reaction_network = get_card(
    title="Reaction Network",
    text="Given the target compounds and their amounts, visualize the suggested reaction network.",
    link_text="Reaction Network",
    link_path="/reaction_network",
    width=20,
)
subpage_cards.append(card_reaction_network)

card_operation_network = get_card(
    title="Operation Network",
    text="Given a quantified reaction network, visualize the operations that realize the reactions.",
    link_text="Operation Network",
    link_path="/operation_network",
    width=20,
)
subpage_cards.append(card_operation_network)

card_scheduler = get_card(
    title="Scheduler",
    text="The optimal schedule solved as a FJSS given a list of functional modules.",
    link_text="Scheduler",
    link_path="/scheduler",
    width=20,
)
subpage_cards.append(card_scheduler)

#
# card_junior_simulator = get_card(
#     title="Workstation Simulator",
#     text="Given an operation graph, Simulate its operations on the digital twin of a workstation.",
#     link_text="Workstation Simulator",
#     link_path="/workstation_simulator",
#     width=20,
# )

layout = dbc.Row(
    dbc.Col(
        subpage_cards,
        className="justify-content-center col-8 d-flex mx-auto"
    ),
    className="align-items-center",
    style={"min-height": "50vh"}
)

