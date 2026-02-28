SHELL := /bin/bash

CODEX_HOME ?= $(HOME)/.codex
CODEX_SKILLS_DIR ?= $(CODEX_HOME)/skills
LOCAL_SKILLS_DIR ?= docs/skills

SIM_SKILLS := sim-plan sim-review sim-implement
SIM_SELECTOR_OVERLAY_TEST := tests_sim/sim/test_chemistry_overlay.py::test_kg_query_selector_can_filter_by_smiles

.PHONY: sync-skill-sim-plan sync-skill-sim-review sync-skill-sim-implement sync-skills-sim test-sim-selector-overlay-stress _sync-skill

sync-skill-sim-plan:
	@$(MAKE) _sync-skill SKILL=sim-plan

sync-skill-sim-review:
	@$(MAKE) _sync-skill SKILL=sim-review

sync-skill-sim-implement:
	@$(MAKE) _sync-skill SKILL=sim-implement

sync-skills-sim: sync-skill-sim-plan sync-skill-sim-review sync-skill-sim-implement
	@echo "Synced SIM skills: $(SIM_SKILLS)"

test-sim-selector-overlay-stress:
	@set -euo pipefail; \
	n="$${N:-20}"; \
	echo "Running selector/overlay stress check $$n times: $(SIM_SELECTOR_OVERLAY_TEST)"; \
	for i in $$(seq 1 "$$n"); do \
		echo "[$$i/$$n]"; \
		PYTHONPATH=. pytest -q "$(SIM_SELECTOR_OVERLAY_TEST)"; \
	done

_sync-skill:
	@set -euo pipefail; \
	if [[ -z "$(SKILL)" ]]; then \
		echo "SKILL is required"; \
		exit 1; \
	fi; \
	src="$(LOCAL_SKILLS_DIR)/$(SKILL)"; \
	dst="$(CODEX_SKILLS_DIR)/$(SKILL)"; \
	if [[ ! -d "$$src" ]]; then \
		echo "Missing skill source: $$src"; \
		exit 1; \
	fi; \
	mkdir -p "$(CODEX_SKILLS_DIR)"; \
	rm -rf "$$dst"; \
	mkdir -p "$$dst"; \
	cp -a "$$src"/. "$$dst"/; \
	echo "Synced $$src -> $$dst"
