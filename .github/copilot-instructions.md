# Copilot Instructions

## Repository Overview

This is a research monorepo with two layers:

1. **`csdid` — a published Python package** (root-level) implementing Callaway & Sant'Anna (2021) Difference-in-Differences with multiple time periods. This is the primary codebase that gets published to PyPI.
2. **Independent subproject directories** (e.g., `MultiModal/`, `GameLLM/`, `musicgen/`, `data-formulator/`, `causal-demand/`) — each is a self-contained research project with its own stack, dependencies, and AI instruction files (`CLAUDE.md`, `AGENTS.md`, or `.github/copilot-instructions.md`). Always read the subproject's instruction file before working in that directory.

## csdid Package

### Architecture

- `csdid/att_gt.py` — `ATTgt` class: the main entry point. Call `.fit(est_method='dr')` to estimate group-time ATTs, then `.aggte(typec='dynamic'|'calendar'|'group'|'simple')` for aggregation.
- `csdid/attgt_fnc/` — estimation internals: data preprocessing (`preprocess_did.py`), doubly-robust trimming (`drdid_trim.py`), ATT computation (`compute_att_gt.py`).
- `csdid/aggte_fnc/` — aggregation logic for event-study, calendar, and group-level effects.
- `csdid/plots/gplot.py` — visualization using `plotnine`.
- `csdid/utils/` — multiplier bootstrap (`mboot.py`) and helper functions (`bmisc.py`).
- External dependency: `drdid` package (installed from `git+https://github.com/d2cml-ai/DRDID`).

### Build & Test

```bash
# Install dependencies
pip install -r requirements.txt
pip install git+https://github.com/d2cml-ai/DRDID

# Install the package in development mode
pip install -e .

# Run the canonical suite (no R required; missing reference fixtures fail).
pytest csdid/test_csdid

# Run a single test / file
pytest csdid/test_csdid/test_att_gt.py
pytest "csdid/test_csdid/test_r_parity.py::test_attgt_matches_r"

# Compare against clean R source checkouts without updating pinned fixtures.
python csdid/test_csdid/upstream_sync/run.py --did /path/to/did \
  --drdid /path/to/DRDID --output /tmp/new-upstream-comparison
```

The `csdid/test_csdid/` suite is the canonical test set. Behavioral tests are
distinct from direct numerical comparisons. Pinned current-source R fixtures
cover 25 configurations in both Python modes, including full influence functions;
historical R 2.5.0 and JEL references remain separately identified. Missing
fixtures fail; ordinary tests require no R installation. The upstream workflow
generates temporary references from fresh GitHub checkouts, reports changed R
tests and fails on numerical/coverage drift without automatically accepting it.
See `csdid/test_csdid/test_suite_desc.md`, `upstream_sync/README.md` within that
directory, and `FASTER_MODE.md` for scope and reproduction instructions.

### Key Conventions

- Version is defined in `csdid/_version.py` as `__version`; update it there for releases.
- The README is auto-generated: edit `README.qmd` (Quarto), not `readme.md` directly. CI renders and commits the output on push to `main`.
- PyPI publishing is automated via CI on push to `main` using `python setup.py sdist` + Twine.

## Subprojects

Each subproject is independent. Before working in one, check for and read its instruction file:

- `CLAUDE.md` or `AGENTS.md` at the subproject root
- `.github/copilot-instructions.md` within the subproject

These files contain project-specific architecture, commands, and conventions that override anything here.
