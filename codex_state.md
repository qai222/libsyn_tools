# Codex State — sim correctness polish (post-refactor)

Date: 2026-01-09

## Baseline
- `pytest -q tests_sim` passes before starting these tasks.

## Correctness issues selected to fix (simple approach)
1) Selector docstring / contract mismatch
- `_atomic_get_and_lock` doc/docstring indicates "reinsert at front", but implementation uses `store.put` (FIFO append).
- This is misleading and can cause future correctness regressions. Fix: align docs/comments to actual behavior (no behavioral change).

2) Functional properties stored as sets + `next(iter(set))` usage
- Several places treat a set-valued property as functional (e.g., pool_type/capacity) by calling `next(iter(...))`.
- If multiple values exist (data bug), behavior becomes nondeterministic. Fix: explicit helper enforcing cardinality (raise ValueError on 0 or >1) and use it in key spots.

3) Report correctness: makespan and in-progress classification
- `build_report()` previously used only OPERATION_END for completion/makespan.
- If operations ABORT/INTERRUPT, report may misclassify them as "in progress" and makespan may be None.
- Fix: define TERMINAL events = END/ABORT/INTERRUPT and compute in-progress and makespan accordingly; also guard empty history.

4) Deadlock limitation documentation
- Even with stable ordering, multi-resource selection from small pools can deadlock (classic "each holds one, waits for another").
- Non-overengineering fix: document limitation clearly in code comments and/or module docs (no behavior change).

## Semantics decisions (for these tasks)
- These tasks should not change simulator runtime semantics except report outputs and fail-fast validation for bad data (multiple values in functional sets).
- For functional-set validation, raising ValueError is acceptable because multiple values represent invalid state.

## Task log
(append entries)

### Template
- Task N: <title>
  - Status: DONE / IN PROGRESS
  - Files changed:
  - Tests added/updated:
  - Notes:

## Completion checklist
- [ ] Docs: selector reinsertion doc matches behavior; deadlock limitation documented
- [ ] Validation helper added; key functional-set lookups use it; tests cover 0/multi values
- [ ] build_report handles empty history; terminal events; makespan computed from terminal events; tests cover abort/interrupt-only runs
- [ ] `pytest -q tests_sim` passes
