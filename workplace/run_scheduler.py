import argparse
import os
import random
import shutil
from loguru import logger
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


FMS_0 = _get_fms(
    module_settings={
        OperationType.Heating: (1, 2),  # capacity, module_number, default to be (1, 1)
        OperationType.TransferLiquid: (1, 1),
        OperationType.MakeSolution: (1, 2),
        OperationType.ConcentrationAndPurification: (1, 2),
    },
)

FMS_1 = _get_fms(
    module_settings={
        OperationType.Heating: (2, 2),  # capacity, module_number, default to be (1, 1)
        OperationType.TransferLiquid: (1, 1),
        OperationType.MakeSolution: (1, 2),
        OperationType.ConcentrationAndPurification: (1, 2),
    },
)

FMS_2 = _get_fms(
    module_settings={
        OperationType.Heating: (3, 2),  # capacity, module_number, default to be (1, 1)
        OperationType.TransferLiquid: (1, 1),
        OperationType.MakeSolution: (1, 2),
        OperationType.ConcentrationAndPurification: (1, 2),
    },
)

FMS_LIST = [
    FMS_0, FMS_1, FMS_2,
]


def run_one(runs_folder: FilePath, libs_folder: FilePath, prefix: str, x: int, y: int, fms_index: int,
            has_work_shifts: bool,
            # temperature_threshold: float,
            time_limit: int, gb_threads: int,
            baseline_only: bool
            ):
    fms = FMS_LIST[fms_index]

    lib_dir = os.path.join(libs_folder, f"{prefix}-{x:02}-{y:02}")
    if not os.path.exists(lib_dir):
        logger.critical(f"Library folder not found: {lib_dir}")
        return
    run_dir = os.path.join(runs_folder, f"{prefix}-{x:02}-{y:02}-{fms_index}-{int(has_work_shifts)}")
    if os.path.exists(run_dir):
        logger.critical(f"Run folder already exists: {lib_dir}")
        return

    shutil.copytree(lib_dir, run_dir)

    workflow = Workflow(routes_file="", work_folder=run_dir)

    json_dump([fm.model_dump() for fm in fms], os.path.join(run_dir, workflow.functional_modules_json))
    workflow.export_scheduler_input(rng=random.Random(42), dummy_work_shifts=has_work_shifts)
    if not baseline_only:
        workflow.export_solver(baseline=False, time_limit=time_limit, gb_threads=gb_threads)
    workflow.export_solver(baseline=True, time_limit=time_limit, gb_threads=gb_threads)


parser = argparse.ArgumentParser(description='library synthesis scheduler examples')
parser.add_argument('--n-target', type=int, help='number of targets')
parser.add_argument('--network-index', type=int, help='index of the generated reaction network')
parser.add_argument('--library', type=str, help='name of the library, VS or FDA')
parser.add_argument('--gb-timelimit', type=int, default=3600, help='time limit for Gurobi')
parser.add_argument('--gb-threads', type=int, default=12, help='threads for Gurobi')
parser.add_argument('--runs-folder', type=str, default="RUNS", help='runs folder path')
parser.add_argument('--libs-folder', type=str, default="LIBS", help='libs folder path')
parser.add_argument('--fms-index', type=int, help='index of the functional module set')
parser.add_argument('--work-shift', action="store_true", help='if include work shifts')
parser.add_argument('--baseline-only', action="store_true", help='if only run the baseline scheduler')

args = parser.parse_args()


def main():
    run_one(
        runs_folder=args.runs_folder,
        libs_folder=args.libs_folder,
        prefix=args.library,
        x=args.n_target,
        y=args.network_index,
        fms_index=args.fms_index,
        has_work_shifts=args.work_shift,
        time_limit=args.gb_timelimit,
        gb_threads=args.gb_threads,
        baseline_only=args.baseline_only,
    )


if __name__ == '__main__':
    main()
