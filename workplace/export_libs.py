import os
import random
from itertools import groupby

import pandas as pd

from libsyn_tools.chem_schema import ReactionNetwork, OperationNetwork
from libsyn_tools.opt.example_routes import ROUTES_PATH_VS, ROUTES_PATH_FDA
from libsyn_tools.utils import FilePath
from libsyn_tools.utils import json_load
from libsyn_tools.workflow import Workflow

"""
export reaction and operation networks of VS-x and FDA-x
"""


def _parse_set(routes_file: FilePath):
    routes = json_load(routes_file)
    routes = {k: routes[k] for k in routes if len(routes[k]['Reactions']) > 1}  # exclude orphans
    targets = sorted(routes.keys())
    return targets


FDA_SET = _parse_set(ROUTES_PATH_FDA)
VS_SET = _parse_set(ROUTES_PATH_VS)


def export_main(up_to_x: int, up_to_y: int, prefix: str, export_to: FilePath):
    if prefix == "FDA":
        targets = FDA_SET
        routes_file = ROUTES_PATH_FDA
    elif prefix == "VS":
        targets = VS_SET
        routes_file = ROUTES_PATH_VS
    else:
        raise ValueError
    assert len(targets) >= up_to_x >= 1
    assert up_to_y >= 1

    rng = random.Random(0)

    for x in range(1, up_to_x + 1):
        for y in range(1, up_to_y + 1):
            wdir = os.path.join(export_to, f"{prefix}-{x:02}-{y:02}")
            os.makedirs(wdir, exist_ok=True)
            targets_x_y = rng.sample(targets, k=len(targets))[:x]
            workflow = Workflow(
                routes_file=routes_file,
                work_folder=wdir,
                scraper_output=None,
            )
            workflow.export_reaction_network(query_askcos=True, dummy_quantify=True, specified_targets=targets_x_y)
            workflow.export_operation_network(rng=rng)


def export_summary_table(prefix: str, export_folder: FilePath, temperature_threshold=30):
    assert prefix in ("FDA", "VS")
    records = []
    for x in range(100):
        for y in range(100):
            lib_dir = os.path.join(export_folder, f"{prefix}-{x:02}-{y:02}")
            if not os.path.exists(lib_dir):
                continue
            reaction_network = ReactionNetwork(**json_load(lib_dir + "/reaction_network.json"))
            operation_network = OperationNetwork(**json_load(lib_dir + "/operation_network.json"))

            # calculate temperature bins
            temperatures = [r.temperature for r in reaction_network.chemical_reactions]
            temperatures = sorted(temperatures)
            groups = []
            for _, g in groupby(temperatures, lambda tt: (tt - 1) // temperature_threshold):
                groups.append(list(g))

            record = {"library_name": os.path.basename(lib_dir), "# temperature bins": len(groups)}
            record.update(reaction_network.summary)
            record.update(operation_network.summary)
            records.append(record)
    return pd.DataFrame.from_records(records)


if __name__ == '__main__':
    EXPORT_LIBS_PATH = os.path.join(os.getcwd(), "LIBS")
    X_MAX = 4
    Y_MAX = 10
    TEMPERATURE_THRESHOLD = 30
    PREFIX = "FDA"

    export_main(up_to_x=X_MAX, up_to_y=Y_MAX, prefix=PREFIX, export_to=EXPORT_LIBS_PATH)
    DF = export_summary_table(prefix=PREFIX, export_folder=EXPORT_LIBS_PATH,
                              temperature_threshold=TEMPERATURE_THRESHOLD)
    DF.to_csv(f"LIBS_{PREFIX}.csv", index=False)
