# Codex State — sim correctness hardening (2026-01-09)

## Baseline
- `pytest -q tests_sim` passes before starting these tasks.

## Issues selected as "worth fixing" (simple + robust)
A) Operation/resources tracking
1) Duplicate resources can double-count provenance/utilization:
   - Operation.resources not deduplicated when same resolved resource bound to multiple roles.

B) Transfer preset validation & tolerance correctness
2) TransferMaterialByPortionSize missing source can crash with AttributeError in error message.
3) TransferMaterialByVolume missing/invalid source can crash (no None/type guard).
4) PortionOfMaterial.get_portion_by_volume tolerance allows values that then fail in get_portion.

C) Overlay completeness for subclassed types
5) CurrentVolumeOverlayProvider omits subclasses of MaterialContainer.
6) ChemistryOverlayProvider omits subclasses of PortionOfMaterial.
7) directly_contained_pom_volume undercounts if PortionOfMaterial subclasses exist.
8) get_directly_contained_individuals omits subclass instances when querying base class.

D) Reporting/export correctness
9) Instance history timestamps ignore abort/interrupt (only OPERATION_END).

E) Interrupt classification (optional-but-worth fixing for correctness observability)
10) Selector cancellation returning None can be misclassified as abort rather than interrupt.

## Non-goals (avoid overengineering)
- No global dead confirmation / advanced deadlock prevention.
- No large redesign of graph model.
- Keep semantics stable; only improve correctness/error clarity/coverage.

## Task log
(append entries)

### Template
- Task N: <title>
  - Status: DONE / IN PROGRESS
  - Files changed:
  - Tests added/updated:
  - Notes:

## Completion checklist
- [ ] Operation.resources deduped (order-preserving); utilization history not double-counted
- [ ] Transfer presets: clear, fail-fast errors on missing/invalid source
- [ ] get_portion_by_volume tolerance consistent with get_portion
- [ ] Overlays include subclass instances (containers + POMs); volume computations include subclasses
- [ ] Instance history timestamps include ABORT/INTERRUPT terminal times
- [ ] Selector cancel is logged/classified as interrupt (not abort) when cancelled during acquisition
- [ ] `pytest -q tests_sim` passes
