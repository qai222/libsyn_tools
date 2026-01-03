from __future__ import annotations


def test_sim_adapters_import_does_not_require_optional_solver_deps() -> None:
    # This should import without requiring `gurobipy`.
    import libsyn_tools.sim.adapters as adapters  # noqa: F401
