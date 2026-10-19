# NAMD-FEP, step by step, in a browser

Interactive notebooks that walk the repository's RBFE workflow end to end —
from a **docked pose or crystal structure** (e.g. the output of
[quantaosun/labodock_binder](https://github.com/quantaosun/labodock_binder))
to a checked ΔΔG — visualising every step: the complex, the sanitised protein,
the ligand pair, the dual-topology hybrid, the run's progress, and the
equilibration/convergence diagnostics behind the final number.

## Launch on Binder

Point mybinder.org at this repository; the build config lives in `binder/`
(repo2docker only discovers the repo root, `binder/`, or `.binder/`), and
`filepath=` opens notebook 00 directly:

```
https://mybinder.org/v2/gh/quantaosun/NAMD-FEP/<branch>?filepath=notebooks/00_binder_setup.ipynb
```

For a local conda env instead, `conda env create -f notebooks/environment.yml`
(an identical copy). Notebook 00 also self-heals: it pip-installs anything the
image lacks.

> Free Binder sessions have **no GPU and ~2 GB RAM**. Everything except the
> NAMD production run works there — preparation, visualisation, sanity checks,
> and the full analysis of the shipped 6I5I example. Notebooks 03–04 detect
> the missing VMD/NAMD binaries, print the exact commands for a GPU box, and
> skip the cells that need them.

## The notebooks

| notebook | step | visualises |
|---|---|---|
| `00_binder_setup.ipynb` | environment check + install | a 3D smoke test |
| `01_complex_visualisation.ipynb` | fetch/upload a complex, `rbfe-prep report` | the protein–ligand complex in 3D |
| `02_prepare_and_sanitise.ipynb` | strip protein, extract + parameterise both ligands, `system.ini`, `rbfe similarity` | sanitised protein, 2D/3D ligand pair overlay |
| `03_hybrid_and_system.ipynb` | `rbfe hybrid` → `build` → `inputs` → `check` | hybrid coloured by its −1/0/+1 B-factor marks; the solvated boxes |
| `04_run_monitor.ipynb` | `rbfe run` (GPU) / `rbfe status` | progress grid, live per-window ΔG preview |
| `05_equilibration_and_analysis.ipynb` | the full audit, interactively | equilibration traces, BAR per window, hysteresis, stationarity, bootstrap CI, trim scan |

## What a notebook does and does not do

The notebooks **drive the same `rbfe` / `rbfe-prep` commands** the RUNBOOK
documents — they are a guided front end, not a reimplementation. Analysis
reuses `rbfe.analysis` (the audit module) directly, so the numbers in a
notebook and `rbfe audit` are the same numbers.

Everything a notebook writes goes to `_tutorial/` (git-ignored). The frozen
`6I5I_DUAL_FEP/` example and the `systems/` configs are only ever read.

Run order is `00 → 05`, but `04`–`05` stand alone for exploring the shipped
6I5I results, and `01`–`02` stand alone for preparing a new complex.
