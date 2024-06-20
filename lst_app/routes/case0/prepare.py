import os

from libsyn_tools.chem_schema.workflow import NetworkWorkflow

workflow = NetworkWorkflow(
    routes_file="routes.json",
    work_folder=os.path.abspath(os.getcwd()),
    scraper_output="scraper_output.json",
    query_askcos=False,
    sample_seed=42,
    sample_n_target=5,
)

workflow.export_reaction_network(dummy_quantify=True)
