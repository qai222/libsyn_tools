from __future__ import annotations

import argparse
from pathlib import Path


def _run_quickstart(name: str, out_dir: str | Path | None) -> None:
    if name == "closed_loop":
        from examples.quickstart_closed_loop import run
    elif name == "schedule":
        from examples.quickstart_schedule import run
    elif name == "semantic_selector":
        from examples.quickstart_semantic_selector import run
    else:
        raise ValueError(f"Unknown quickstart: {name}")

    run(out_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description="libsyn_tools.sim CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    quickstart_parser = subparsers.add_parser("run-quickstart", help="Run a quickstart example")
    quickstart_parser.add_argument(
        "name",
        choices=("closed_loop", "schedule", "semantic_selector"),
        help="Quickstart example to run",
    )
    quickstart_parser.add_argument(
        "--out",
        required=True,
        help="Output directory for report artifacts",
    )

    args = parser.parse_args()

    if args.command == "run-quickstart":
        _run_quickstart(args.name, Path(args.out))


if __name__ == "__main__":
    main()
