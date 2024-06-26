import os

from libsyn_tools.workflow import Workflow

workflow = Workflow(
    routes_file="routes.json",
    work_folder=os.path.abspath(os.getcwd()),
    scraper_output="scraper_output.json",
    query_askcos=False,
    sample_seed=42,
    sample_n_target=5,
)

if __name__ == '__main__':
    workflow.export_reaction_network(dummy_quantify=True)
    workflow.export_operation_network()
    workflow.export_scheduler(max_capacity=1, max_module_number_per_type=1)
