# ### THIS IS THE START OF CONTENT OF tests_sim/sim/test_examples_import.py ###
from __future__ import annotations

import importlib
import sys
from pathlib import Path


def test_quickstart_modules_importable() -> None:
    examples_dir = Path(__file__).resolve().parents[2] / "examples"
    sys.path.insert(0, str(examples_dir))
    try:
        importlib.import_module("quickstart_closed_loop")
        importlib.import_module("quickstart_schedule")
        importlib.import_module("quickstart_semantic_selector")
    finally:
        sys.path.pop(0)


def test_cli_module_importable() -> None:
    importlib.import_module("libsyn_tools.sim.__main__")
# ### THIS IS THE END OF CONTENT OF tests_sim/sim/test_examples_import.py ###
