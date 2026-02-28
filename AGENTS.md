# AGENTS.md

## Scope and intent
This file guides coding agents working in this repository. Keep changes small, test-backed, and aligned with current project priorities.

## Project map
Primary code lives in:

- `libsyn_tools/chem_schema`: shared chemistry and operation network schemas.
- `libsyn_tools/opt`: schedule optimization (published work; treat as locked).
- `libsyn_tools/sim`: simulator and dynamic knowledge-graph execution layer (active development focus).

Supporting areas:

- `tests/`: chem schema, opt, and utils tests.
- `tests_sim/`: simulator-focused tests.
- `docs/sim/EXTENDING.md`: extension workflow for custom simulator operations.
- `examples/`: runnable examples and exploratory scripts.
- `old_sim/`: historical code; do not use as implementation target.

## Priority and change policy
- Prefer changes in `sim` and `chem_schema` when implementing new behavior.
- Assume `opt` is stable/locked unless the user explicitly requests `opt` work.
- If `sim` changes require schema updates, keep `chem_schema` edits minimal and backwards-compatible.
- Do not perform broad refactors unless requested.

## `twa` integration notes
`twa` is a frozen external dependency and the OGM foundation for `libsyn_tools.sim`.
At a high level, `sim` defines its ontology/entities on top of `twa` and executes operations as transactions over the in-memory KG model.

Latest implementation details and integration patterns are tracked in:
- `docs/sim/twa.md`

When `twa` usage in `sim` is changed, update `docs/sim/twa.md` in the same change.

## How to validate changes
Run the smallest relevant test slice first, then expand only if needed.

- For `sim` changes:
  - `PYTHONPATH=. pytest -q tests_sim`
- For `chem_schema` changes:
  - `PYTHONPATH=. pytest -q tests/test_chem_schema.py tests/test_utils.py`
- For `sim`↔`opt` boundary changes (schedule adapter/import behavior):
  - `PYTHONPATH=. pytest -q tests_sim/test_adapters_optional_import.py`
- For `opt` changes (only when explicitly requested):
  - `PYTHONPATH=. pytest -q tests/test_opt.py`
  - Note: this may require a working `gurobipy` setup/license.

## Anti-staleness rules
- Do not encode version pins, license details, or environment assumptions here; read `requirements.txt` and current tests/docs at runtime.
- Do not rely on generated example outputs (`examples/sim/_out_*`) as source of truth.
- Prefer deriving behavior from current tests and module code over historical docs or archived folders.
- If guidance conflicts with code/tests, treat code+tests as authoritative and update this file in the same change.
