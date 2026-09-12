# RUNBOOK — operating the NAMD RBFE workflow

Operational guide for the current work. The root `README.md` is the original
2021-era public tutorial (LigParGen + Feprepare) and does **not** describe this
pipeline. `NAMD_RBFE_Guide_zh.md` is the detailed Chinese walkthrough and is kept
in sync with this file.

---

## 1. Repo map — what is current and what is history

| path | status |
|---|---|
| **`6I5I_DUAL_FEP/`** | **CURRENT.** The system that produced the reported ΔΔG. Prep, configs, per-window results, controller, analysis. |
| `toppar/` | **CURRENT.** CHARMM36 parameters. Includes the project-specific `par_water_ions_clean.prm`. |
| `fep_pipeline/` | ⚠️ **NOT WORKING.** Never ran end-to-end; `run` path is broken (see §7). Only `alignment.py` is usable. |
| `T4L_RBFE/` | 📁 historical. NAMD tutorial system, built with the Feprepare/CHARMM-GUI/LigParGen web services. |
| `6I5I_RBFE/` | 📁 historical. 6I5I attempt via acpype + GAFF. |
| `6I5I_FEP/` | 📁 historical. 6I5I via acpype/GAFF + RDKit atom mapping (`atom_map.json`). Source of the ligand PDBs/SDFs still cited in `DDG_preliminary.md`. |
| `6I5I_deck/` | 📁 presentation slides (reveal.js). |
| `*.namd` at root, `NAMD-FEP-tutorial.pdf`, `NAMD-FEP_local.ipynb` | 📁 historical tutorial material. |

Anything marked historical is kept for provenance. **Do not build on it.**

---

## 2. Quick start — the one entry point

Everything runs from `6I5I_DUAL_FEP/`.

```bash
cd 6I5I_DUAL_FEP

# 1. Check the GPU first (a job started without one wastes the setup time)
nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv

# 2. Start / resume the whole job. Idempotent — safe to re-run any number of times.
bash respawn_controller.sh     # detached, survives the box's ~4.5 h job kills
```

`respawn_controller.sh` → `run_checkpointed.sh` → `fep_run.py`. It is
marker-guarded: every λ-window writes a `.done` marker, so re-running resumes at
the interrupted window instead of replaying.

**Watch progress:**
```bash
tail -f controller_run.log          # controller decisions
bash checkpoint_status.sh &         # one status line every 5 min
python3 fep_run.py windows complex md_forward   # JSON: how much is done
```

### ⚠️ `run_all.sh` is NOT the entry point

It is the original non-resumable launcher: it replays every stage from
`nvt_equil` and overwrites `md_forward`/`md_backward`. It now **refuses to run**
when production state exists, and exits 1. Use `run_checkpointed.sh`.
(`run_all.sh --force` overrides — only for a genuinely fresh start.)

---

## 3. The full pipeline, in order

```bash
cd 6I5I_DUAL_FEP

# --- Preparation (only needed once; regenerates hybrid/, ionized.*, *.namd) ---
python3 prepare_hybrid.py      # inputs/{ref,mut}.{rtf,prm,pdb} -> hybrid/hybrid.{rtf,prm,pdb}
python3 build_system.py        # VMD psfgen + solvate + ionize -> complex/, solvent/
python3 write_fep_inputs.py    # ionized.fep B-factors + fep.tcl + the 8 .namd configs

# --- Equilibration (per leg: nvt then npt) ---
# --- Production: one NAMD process per lambda window, forward then backward ---
#     all handled by the controller; do not launch these by hand
```

Each leg is `nvt_equil → npt_equil → md_forward (0→1) → md_backward (1→0)`.
15 λ-windows per direction, 50 000 equilibration + 250 000 production steps each.

> **Regenerating preparation files desyncs the committed results.** `hybrid/*`
> are checked-in *build products* of `prepare_hybrid.py`. If you re-run step 1,
> `hybrid/` no longer matches the run that produced the reported ΔΔG.

---

## 4. Analysis

```bash
cd 6I5I_DUAL_FEP

python3 audit_fep.py              # full audit: parser validation, BAR, hysteresis,
                                  # stationarity, overlap, block-bootstrap error bars
python3 audit_fep.py --skip-boot  # fast
python3 audit_fep.py --scan       # ΔΔG vs equilibration-trim sensitivity

python3 analyze_fep.py bar \
  complex/md_forward_combined.fepout complex/md_backward_combined.fepout \
  solvent/md_forward_combined.fepout solvent/md_backward_combined.fepout
```

**Current result:** ΔΔG(12H → desmethyl) = **−0.106 ± 0.120 kcal/mol**,
95% CI **[−0.355, +0.118]** — i.e. zero within error. See
`DDG_preliminary.md` §7 for the corrections that produced this number and the
known limitations.

> ⚠️ **Every per-window `.fepout` opens with 99 pre-equilibration rows** out of 500
> (`alchEquilSteps 50000` ÷ `alchOutFreq 500`). NAMD's own `dE_avg`/`dG` columns
> exclude them, and so must any analysis. Both scripts above trim by default;
> `--no-trim` reproduces the old (biased) behaviour.
>
> ⚠️ **Do not use `fep_pipeline.analysis`** — deleted 2026-09-12. Its greedy regex
> captured the cumulative `net change until now` instead of the per-window value.

---

## 5. Reassembling a combined fepout

`<stage>_combined.fepout` is the analysis input, built from the per-window files.
Rebuild with:

```bash
python3 fep_run.py assemble complex md_forward    # and md_backward, solvent/*
```

`assemble` is deliberately strict — it refuses unless every window is both
**complete (500 rows)** and **marked `.done`**, because a killed window leaves a
partial `.fepout` that would otherwise be assembled as a silently truncated
window. It preserves NAMD's own `STEPS OF EQUILIBRATION` and
`Free energy change ...` comment lines; those are what make the output
verifiable against NAMD's internal estimator.

---

## 6. Traps that have already bitten

| trap | what happens | guard now |
|---|---|---|
| `--dry` marking windows done | every later real run skips them; the controller declares the leg complete | removed — `--dry` writes no marker |
| assembling a partial window | silent truncation, biased ΔΔG | `assemble` requires the `.done` marker + 500 rows |
| `fep_run.py run <leg> fwd --start N` with no `--from` | window N is seeded from `npt_equil` instead of w(N−1) | pass `--from window_snapshots/md_forward_w05` explicitly |
| running `run_all.sh` on a finished job | ~14 h of duplicated GPU work, clobbered `.fepout` | refuses unless `--force` |
| κ `temperature` in a restart-chained config | `FATAL ERROR: Cannot specify both an initial temperature and a velocity file` | only `nvt_equil` sets it |
| a second `namd3` starting | splits the single GPU — makes everything *slower* | `wait_gpu_free` gate before every launch |

---

## 7. GPU notes (this box)

- **V100, sm_70, 16384 MiB usable — and shared.** Other processes may already hold
  memory. Check `nvidia-smi` before starting, not just that a GPU exists.
- Run NAMD with **`+p1`**, never `+p8`. In GPU-resident mode throughput scales
  *inversely* with PE count (`+p1` 75.8 ns/day vs `+p8` 13.2). `run_all.sh` and
  `run_checkpointed.sh` both set this; don't "optimise" it back.
- One GPU ⇒ **strictly sequential legs**. `wait_gpu_free` enforces this.
- NAMD peaks at ~590 MiB, so memory is not the constraint — latency is.
- Judge throughput from `TIMING:` sec/step, **not** the cumulative `PERFORMANCE:`
  line (a leading minimization drags that average down for hours).

Binary: `/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3`
VMD: `/home/aistudio/vmd-env/bin/vmd` (not on `PATH` — hardcoded in `build_system.py`)

---

## 8. `fep_pipeline/` — why it is marked broken

A 2643-line package advertised as a Feprepare replacement. It has **never run
end-to-end** (no output, no tests, no history). Three independent breaks in `run`:

1. `system_builder.py:60-62` passes `.prm` files and the hybrid **PDB** to
   psfgen's `topology` and emits no `parameters` lines — it cannot build a PSF.
2. `prepare_fep.py:339` references `alchFile ionized_{leg}.fep`; nothing generates it.
3. `fep_file.build_system_fep` — the function that would — is never called.

Salvageable: `alignment.py` (pure-stdlib Kabsch + B-factor marking) and the
`namd_config.py` templates (which still carry the restart-chain `temperature` bug).
Its `analysis.py` was deleted for returning wrong numbers.
