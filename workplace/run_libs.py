import os
import random
import shutil

from libsyn_tools.utils import FilePath, json_dump
from libsyn_tools.workflow import Workflow, OperationType, FunctionalModule


def _get_fms(module_settings):
    fms = []
    for t in OperationType:
        try:
            capacity, module_number = module_settings[t]
        except KeyError:
            capacity = 1
            module_number = 1
        for i in range(module_number):
            m = FunctionalModule(
                name=f"{t}-{i}",
                can_process=[t, ],
                capacity=capacity,
            )
            fms.append(m)
    return fms


FMS_1 = _get_fms(
    module_settings={
        OperationType.Heating: (1, 2),  # capacity, module_number, default to be (1, 1)
        OperationType.Purification: (1, 2),
        OperationType.TransferLiquid: (1, 2),
        OperationType.TransferSolid: (1, 2),
        OperationType.Concentration: (1, 2),
    },
)

FMS_2 = _get_fms(
    module_settings={
        OperationType.Heating: (3, 2),  # capacity, module_number, default to be (1, 1)
        OperationType.Purification: (1, 2),
        OperationType.TransferLiquid: (1, 2),
        OperationType.TransferSolid: (1, 2),
        OperationType.Concentration: (1, 2),
    },
)

FMS_3 = _get_fms(
    module_settings={
        OperationType.Heating: (3, 3),  # capacity, module_number, default to be (1, 1)
        OperationType.Purification: (1, 1),
        OperationType.TransferLiquid: (1, 2),
        OperationType.TransferSolid: (1, 2),
        OperationType.Concentration: (1, 1),
    },
)

FMS_LIST = [FMS_1, FMS_2, FMS_3]


def run_one(runs_folder: FilePath, libs_folder: FilePath, prefix: str, x: int, y: int, fms_index: int,
            has_work_shifts: bool, temperature_threshold: float, time_limit: int, gb_threads: int):
    fms = FMS_LIST[fms_index]

    lib_dir = os.path.join(libs_folder, f"{prefix}-{x:02}-{y:02}")
    if not os.path.exists(lib_dir):
        return
    run_dir = os.path.join(runs_folder, f"{prefix}-{x:02}-{y:02}-{fms_index}-{int(has_work_shifts)}")
    if os.path.exists(run_dir):
        return

    shutil.copytree(lib_dir, run_dir)

    workflow = Workflow(routes_file="", work_folder=run_dir)

    json_dump([fm.model_dump() for fm in fms], os.path.join(run_dir, workflow.functional_modules_json))
    workflow.export_scheduler_input(rng=random.Random(42), temperature_threshold=temperature_threshold,
                                    dummy_work_shifts=has_work_shifts)
    workflow.export_solver(baseline=False, time_limit=time_limit, gb_threads=gb_threads)
    workflow.export_solver(baseline=True, time_limit=time_limit, gb_threads=gb_threads)


if __name__ == '__main__':
    RUNS_FOLDER = os.path.join(os.getcwd(), "RUNS")
    LIBS_FOLDER = os.path.join(os.getcwd(), "LIBS")
    PREFIX = "FDA"
    TEMPERATURE_THRESHOLD = 30
    TIME_LIMIT = 3600
    GB_THREADS = 48

    X_MAX = 4
    Y_MAX = 10

    for X in range(3, X_MAX + 1):
        for Y in range(1, Y_MAX + 1):
            for IFM in range(len(FMS_LIST)):
                run_one(RUNS_FOLDER, LIBS_FOLDER, PREFIX, X, Y, IFM, has_work_shifts=False,
                        temperature_threshold=TEMPERATURE_THRESHOLD, time_limit=TIME_LIMIT, gb_threads=GB_THREADS)
