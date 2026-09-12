# NAMD-FEP Project

## Overview
A tutorial repo for Relative Binding Free Energy (RBFE) calculations using NAMD.
Calculates ΔΔG of binding for small molecule ligands against a protein target.

### Docs / files generated (register every new file here)
- `NAMD_RBFE_Guide_zh.md` (repo root, 2026-09-06) — single combined Chinese guide:
  why physics-based RBFE stays the gold standard vs AI · Baidu AI Studio env (NAMD
  3.0.3 build + per-session `setup.sh`) · 6I5I manual per-λ-window run workflow.
  Replaces the earlier `NAMD_AI_STUDIO_SETUP.md` + `6I5I_DUAL_FEP/MANUAL_RUN_WORKFLOW.md`
  (both deleted after merge).
- `6I5I_DUAL_FEP/respawn_controller.sh` (2026-09-06) — self-healing relauncher for
  `run_checkpointed.sh`; survives the box's ~4.5 h job kills (exits when
  `ALL LEGS DONE` or after 2 quick-fails). Launch detached, PID → `controller_run.pid`.
- `6I5I_DUAL_FEP/DDG_preliminary.md` (2026-09-08) — theory + process notes for the ΔΔG:
  thermodynamic cycle, dual-topology/λ mapping, EXP vs BAR estimators, how the preview
  number (+0.09 kcal/mol) was obtained, and the two ligand SMILES (12H N–CH₃ → desmethyl
  N–H on a pyrazole N).
- `6I5I_DUAL_FEP/analyze_fep.py` (2026-09-08) — NAMD .fepout analyzer reading raw
  `FepEnergy` rows directly (the old `fep_pipeline/analysis.py` greps for strings that
  never appear in a NAMD fepout → returns zeros). Two estimators, pure Python:
  `exp` = one-sided Zwanzig from forward data (gives a **preview ΔΔG while backward
  legs still run` — biased), `bar` = two-state Bennett per window from fwd+bwd pairs
  (final number). Run: `python3 analyze_fep.py exp <complex_fwd.fepout> <solvent_fwd.fepout>`
  or `bar <cf> <cb> <sf> <sb>`. Caveat: `dE` column = U(l+dl)−U(l), so backward files
  pair to window i by reversed index. **Two bugs fixed 2026-09-09** (initial `bar` run
  returned a nonsense −894 / per-window −59.616 = −100/β for every window, both legs):
  (1) `bar_window` bisection branch was inverted — `resid` is monotonically *increasing*
  in dG, so `resid(mid)>0` means the root is below mid ⇒ must set `hi=mid`, not `lo=mid`;
  (2) `bar_leg` never negated the backward dE — a backward NAMD run reports U(A)−U(B)
  sampled at B, so the B-side forward work is `−dE_backward`. After both fixes BAR per
  window tracks one-sided forward EXP (overlap good) and the as-is backward EXP ≈ −forward
  EXP (valid reverse). Final: dG_complex=+3.917, dG_solvent=+3.974 ⇒ **ΔΔG = −0.057**.
  ⚠️ **Superseded 2026-09-12** — this number averages the per-window equilibration rows;
  see the corrected ΔΔG below and `audit_fep.py`.
- `6I5I_DUAL_FEP/audit_fep.py` (2026-09-12) — **independent audit** of the RBFE, written to
  cross-check `analyze_fep.py` (separate parser, separate numpy BAR). Two things it adds:
  (1) **discards the per-window equilibration samples** — locates the production start from
  NAMD's own `#<N> STEPS OF EQUILIBRATION AT LAMBDA ... COMPLETED` marker (auto-detected as
  index 99 in all 60 windows here, not hardcoded) rather than averaging all 500 rows;
  (2) **validates the extraction chain against NAMD's own free-energy output**: recomputes
  one-sided EXP per window from the raw `dE` column and compares to NAMD's
  `#Free energy change for lambda window [...] is <dG>` comment (and the `dG` column) —
  agreement to 7e-6 kcal/mol over all 30 forward windows. Also reports hysteresis,
  first-half-vs-second-half stationarity, σ(dE)/kT overlap, and a moving-block bootstrap
  error on ΔΔG. Runs on the *per-window* `.fepout` files (`md_*_wNN.fepout`), bypassing
  `fep_run.py assemble` entirely. Read-only, numpy only.
  `python3 audit_fep.py` (full) · `--skip-boot` (fast) · `--scan` (ΔΔG vs trim) · `--trim <n>`.
  Note the marker-aware parser relies on NAMD's `#... STEPS OF EQUILIBRATION` and
  `#Free energy change` comment lines, which `assemble` strips — so it reads the per-window
  files (and the salvaged single-run `md_forward.fepout` for complex/forward w00–w05), never
  the `*_combined.fepout`.
- `/home/aistudio/work/setup.sh` (edited 2026-09-06) — now also exports `$NAMD`
  (gpuresident binary), `$NAMD_CPU`, `$CUDA_HOME`, PATH for AI Studio ephemeral
  runtimes; `source` it at the start of every session.

## ⏳ 6I5I production run — LIVE + checkpointing (2026-09-06, resumed)

Full dual-topology RBFE (6I5I CLK1 · ligand 12H, N–CH₃ → N–H). GPU-resident namd3, `+p1`,
**strictly sequential legs** (single V100 — never split the GPU between two NAMD jobs).
Controller log: `6I5I_DUAL_FEP/controller_run.log` (controller PID in `controller_run.pid`);
per-stage `{complex,solvent}/<stage>.log`; check alive via `ps aux | grep namd3`,
GPU via `nvidia-smi`.

**Current status (2026-09-09): DONE ✅.** Relaunched the respawn controller 2026-09-09 19:33
(pid 7628 → `controller_run.pid`) after the 09-08 23:28 kill; it re-ran the one unmarked
window (`solvent/md_backward_w08`, ~15 min) then w09–w14 at ~14.5 min/window ≈ 57 ns/day,
logging `ALL LEGS DONE` at **21:14:52** (wrapper exit-0). All four legs assembled
(`*_combined.fepout`, 7500 `FepEnergy` rows / 15 windows each, no truncation).

**Corrected result (2026-09-12) — ΔΔG = −0.106 ± 0.120 kcal/mol, i.e. zero.** The earlier
"−0.057" came from `analyze_fep.py`, which averages **all 500 rows** of each window. Only
**401** are production data: with `alchEquilSteps 50000` at 2 fs = 100 ps and
`alchOutFreq 500` = 1 ps, the first **99 rows of every window (steps 500–49500) are
per-window equilibration** and NAMD's own accumulator resets at `stepInRun ==
alchEquilSteps` (`Controller::outputFepEnergy`), so its `dE_avg`/`dG` columns use the 401
only. Confirmed to 4 decimals: `mean(rows[99:])` reproduces NAMD's `dE_avg` exactly, and
`EXP(rows[99:])` reproduces NAMD's own per-window `dG` to **7e-6 kcal/mol across all 30
forward windows**. With the trim: dG_complex +3.882, dG_solvent +3.988 ⇒ **ΔΔG = −0.106**;
block-bootstrap 95% CI **[−0.355, +0.118]** → indistinguishable from zero. Run
`python3 audit_fep.py` for the full audit (see that file's entry above).

**Open convergence issues (2026-09-12, not yet acted on):**
1. ΔΔG is a near-cancellation of two ~2.3 kcal/mol halves that flip sign at λ≈0.5 (exactly
   `alchElecLambdaStart 0.5`, where the methyl's charge finishes vanishing and the N–H's
   starts appearing): Σ per-window diff ≈ −2.37 (λ<0.5) + 2.26 (λ>0.5). Fragile by
   construction.
2. The **complex leg is not stationary**: Σ[late-half − early-half] = **+0.365** kcal/mol
   (solvent −0.007). Largest single window: w11, +0.258. Still drifting across the 500 ps
   production window ⇒ under-equilibrated/under-sampled.
3. ΔΔG moves ±0.1 kcal/mol with the trim alone (−0.057 at trim 0 → −0.141 at trim 150 →
   +0.067 at trim 300) — a systematic comparable to the statistical error.
4. **Backward w00 restarts from `npt_equil.coor` at λ=1.0**, not from forward w14's endpoint
   — the backward leg's environment is equilibrated around the *methylated* ligand while
   being sampled at the *desmethyl* state. Cheap fix: seed it from `md_forward_w14`.
5. 15 windows × 500 ps = 7.5 ns/leg is short for RBFE (2–5 ns/window is typical).

**State (resumed 2026-09-06 10:14 CST):** the original single-run `complex/md_forward` was
TERM'd at step 1,707,500 (22:52, mid-window-6); complex + solvent legs now run under the fixed
`run_checkpointed.sh`, ONE NAMD invocation per λ-window via `fep_run.py`:
1. `complex/forward` resumed at window 6 from `complex/window_snapshots/md_forward_w06` (a clean
   restart at step 1,500,000). Windows 0–5 are salvaged from the crashed single-run fepout.
2. then `complex/backward`, then `solvent` nvt(✓ done)/npt/fwd/bwd — each marker-guarded.
Remaining ≈ 6.0M (complex) + 7.55M (solvent) steps ≈ **~13 h, ETA ~Mon 00:00 CST** at
measured ~41 ns/day (complex) / ~59 ns/day (solvent). Measured 2026-09-06: complex w06 window
41.6 ns/day (0.00415 s/step) steady.

**Two bugs fixed 2026-09-06 (both would have silently wasted the run):**
1. `fep_run.py` — `FLAGS = ["+p1", "+devices 0"]` passed `"+devices 0"` as ONE argv token to
   `subprocess.run` → Charm RTS couldn't parse it → NAMD died `FATAL ERROR: Unknown command-line
   option +devices 0`. Fixed: `["+p1", "+devices", "0"]` (two tokens). A shell word-splits; a list
   does not.
2. `run_checkpointed.sh` — `PY="python3 fep_run.py"` used as a single command name → "command not
   found". Fixed: `run_fep(){ python3 fep_run.py "$@"; }` (function) at every call site. Also added
   `wait_gpu_free` (blocks while ANY `namd3` runs) before every launch, and run_production now
   aborts (`exit 4`) if a leg's window runs or assemble fail instead of spuriously logging "done".
Do **not** relaunch `run_all.sh` (replays from nvt and is not crash-safe).

**Checkpointing already in place:**
1. NAMD binary restart every 500 steps (~3 s) → per-window `*.coor/.vel/.xsc` always fresh (+`.BAK/.old`).
2. `watch_windows.sh <leg> <stage> <spw>` — copies clean `coor/vel/xsc` into
   `<leg>/window_snapshots/<stage>_wNN.*` at every λ-window boundary, so a killed direction can
   resume at the current window, not window 0. (Per-window runs don't need it — each window keeps
   its own restart.)
3. `fep_run.py` — per-λ-window driver (idempotent, `.md_*_wNN.done` markers): `run <leg> <fwd|bwd>`,
   `assemble <leg> <stage>` (builds canonical `<stage>_combined.fepout`; note NAMD writes only ONE
   `#NEW FEP WINDOW` line per process — windows are delimited by step buckets of 250 000),
   `windows <leg> <stage>` (probe). Verified: templating + step-bucket parser.

`alchEquilSteps 50000` is **per-window** ("equilibration before data collection in the alchemical
window", `SimParameters.C:1267`), and per-window configs set the same `alchLambda/Lambda2/run` as
the single-run `runFEP` — so salvaged single-run windows 0–5 and resumed per-window windows 6–14
are protocol-identical.

Also see: “Gotchas learned (6I5I prep)” below (temperature vs binvelocities fix etc.).


## fep_pipeline — Self-Contained RBFE Preparation Toolkit

`fep_pipeline/` is a Python 3 package that replaces the Feprepare web server,
Maestro, and LigParGen with a single local workflow. It generates everything
needed to run an RBFE simulation with NAMD + VMD + Python.

### ⚠️ Status: NOT WORKING — do not use `run` (audited 2026-09-12)

This package has **never been run end-to-end**. There is no `fep_status.json`, no
output directory, no test, no packaging metadata, and no git history before the
2026-09-12 safety commit. An audit found three independent breaks in the `run`
path, so the command below does not work:

1. `system_builder.py:60-62` passes `.prm` **parameter** files and the hybrid
   **PDB** to psfgen's `topology` command and emits no `parameters` lines at all
   — it cannot produce a PSF. (Compare the working `6I5I_DUAL_FEP/build_system.py:59-60`,
   which correctly does `topology top_all36_prot.rtf` + `topology hybrid.rtf`.)
2. `prepare_fep.py:339` writes `alchFile ionized_{leg}.fep`, but nothing in the
   package ever generates that file — NAMD would die on startup.
3. `fep_file.build_system_fep`, the function that *would* generate it, is defined
   and never called anywhere.

`ligand_param.py` is unreachable code and would crash if called (`structure.save(
..., format="charmm_rtf")` — ParmEd has no such format). `system_builder.py`
re-reads the whole protein PDB per chain, duplicating atoms on multi-chain input.

**What is actually salvageable:** `alignment.py` (pure-stdlib Horn–Kabsch
superposition + B-factor marking) and the `namd_config.py` templates — the latter
still emit `set temp` alongside `binvelocities` in the NPT/production templates,
i.e. the restart-chain bug recorded below.

**`analysis.py` was deleted 2026-09-12.** It returned *wrong numbers*, not zeros:
its greedy regex `re.search(r"Free energy change.*is\s+([-\d.]+)")` captured the
cumulative `net change until now` instead of the per-window value. Measured on
`6I5I_DUAL_FEP/complex/md_forward.fepout`: **−16.821 where the truth is −5.455**
(5 of 6 windows wrong). Use `6I5I_DUAL_FEP/audit_fep.py` instead.

**Key modules (reality, not aspiration):**
| Module | Status |
|--------|--------|
| `alignment.py` | real, usable — pure stdlib Kabsch + hybrid B-factor marking |
| `namd_config.py` | real templates; needs the restart-chain `temperature` fix |
| `system_builder.py` | **broken** — psfgen script is invalid |
| `fep_file.py` | real code, **never called** |
| `ligand_param.py` | dead code, would crash |
| `prepare_fep.py` | CLI is wired; the `run` path is not |
| ~~`analysis.py`~~ | **deleted** — returned double-counted ΔG |

**Usage that does work** (ligand alignment only):
```bash
python -m fep_pipeline align --reference ref.pdb --mobile mut.pdb --output aligned.pdb
```

The 6I5I work never used this package at all — `grep -rn fep_pipeline
6I5I_DUAL_FEP/` returns nothing. For ΔΔG analysis use:
```bash
python3 6I5I_DUAL_FEP/audit_fep.py          # full audit + error bars
python3 6I5I_DUAL_FEP/analyze_fep.py bar <cf> <cb> <sf> <sb>
```

**Requirements:** Python 3.8+, numpy, scipy, VMD (for system building),
NAMD pre-compiled binary (for simulation).

**License:** MIT (this package). FEPrepare (which inspired the architecture) is GPLv3 —
we wrote an independent, clean-room implementation.

## Environment / Hardware

### CPU
- 64 logical cores (`nproc` = 64; 2 sockets × 16 cores × 2 PUs)
- AVX2 support
- ~503 GB RAM

### GPU / Accelerators
- **1× NVIDIA Tesla V100-SXM2-32GB** (compute capability **7.0** / sm_70), driver 525.125.06 (CUDA 12.0 driver)
- **CUDA 11.8 toolkit** at `/usr/local/cuda-11.8` (`nvcc` V11.8.89); `/usr/local/cuda` → `/etc/alternatives/cuda` → cuda-11.8
- NOTE: an earlier CLAUDE.md described 16× Iluvatar BI-V150S with a CUDA 10.2 compat layer — that is **stale/wrong** for the current box. `/usr/local/corex*` does not exist. Treat this as a standard NVIDIA V100 + CUDA 11.8 machine.

### OS / Software
- Ubuntu 20.04.6 LTS
- GCC 9.4.0, G++ 9.4.0
- **No gfortran** installed
- OpenMPI 4.0.3 (`mpicxx`, `mpicc`, `mpirun` available)
- CMake 3.24.2 at `/usr/local/bin/cmake`
- Conda 23.1.0 at `/opt/conda/bin/conda`
- **No sudo access** — must use conda or user-space installs
- Disk: 2.7 TB free at `/home/aistudio/`

### Key Paths
| Resource | Path |
|----------|------|
| CUDA toolkit | `/usr/local/cuda-11.8/` (`nvcc` V11.8.89; `/usr/local/cuda` symlinks here) |
| CUDA libs | `/usr/local/cuda-11.8/lib64/` |
| GPU monitor | `nvidia-smi` (Tesla V100-SXM2-32GB, sm_70) |
| NAMD binary (**use this**) | `/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3` (GPU-resident, 3.9× faster) |
| NAMD binary (old, non-GPU-resident) | `/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++/namd3` |
| Tcl (conda) | `/opt/conda/lib/libtcl8.6.so`, `/opt/conda/lib/tclConfig.sh` |
| FFTW3 runtime | `/usr/lib/x86_64-linux-gnu/libfftw3.so.3` (no dev headers) |
| Conda | `/opt/conda/bin/conda` |

## NAMD Installation — Status: DONE (GPU/CUDA build working)

### Result
NAMD 3.0.3 compiled from source with CUDA, targeting the V100 (sm_70). Binary at
`/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++/namd3`. Smoke test (`./namd3 +p4 src/alanin`)
runs to completion and detects the CUDA device.

### Build recipe (what actually worked)
1. `NAMD_3.0.3_Source.tar.gz` (at `~/`) bundles Charm++ 8.0.0 (`charm-8.0.0.tar` inside) — no separate download.
2. Precompiled TCL + FFTW from `http://www.ks.uiuc.edu/Research/namd/libraries/` (reachable; the *binary* page is 403):
   `fftw-linux-x86_64.tar.gz`, `tcl8.6.13-linux-x86_64.tar.gz`, `tcl8.6.13-linux-x86_64-threaded.tar.gz`.
3. Build Charm++: `cd charm-8.0.0 && ./build charm++ multicore-linux-x86_64 --with-production -j8`
4. Configure: `./config Linux-x86_64-g++ --charm-arch multicore-linux-x86_64 --with-cuda --cuda-prefix /usr/local/cuda-11.8 --cuda-gencode arch=compute_70,code=sm_70`
5. `make -j32`

### Critical gotcha (FFTW non-PIC)
The precompiled `libsfftw.a`/`libsrfftw.a` are FFTW 2.1.x built **without `-fPIC`**, so the final link
fails on Ubuntu 20.04 (`gcc --enable-default-pie`) with `relocation R_X86_64_32 ... can not be used when
making a PIE object`. `-no-pie` did NOT fix it through charmc. Fix: rebuild FFTW 2.1.5 with
`--enable-float CFLAGS="-O3 -fPIC"`, then rename the outputs to `libsfftw.a`/`libsrfftw.a`/`sfftw.h`/`srfftw.h`
(and fix `srfftw.h`'s `#include "fftw.h"` → `#include <sfftw.h>`). Installed in `NAMD_3.0.3_Source/fftw/`;
originals backed up in `~/fftw-orig/`.

### Not needed (contradicts earlier notes)
- **gfortran** — NAMD 3.x has no Fortran code.
- **csh/tcsh** — the `config` script is `#!/bin/bash`.
- **Iluvatar/CoreX** — this box is a standard NVIDIA V100 + CUDA 11.8 (see GPU section).

## Recent Commits (this repo)
- `028cf21` — Fix analyze_fep.py equilibration bias (ΔΔG −0.057 → −0.106); guard run_all.sh (2026-09-12)
- `e4d0adf` — Safety checkpoint: bring untracked operational code under version control (2026-09-12)
- `089aabb` — Fix NAMD FEP configs: alchDecouple off for RBFE, non-zero alchEquilSteps for ParseFEP, shell portability (2026-07-29)
  - ⚠️ The commit subject is misleading: **`alchDecouple off` is NAMD's default**
    (`SimParameters.C:1250-1251`, `&alchDecouple, FALSE`), as are
    `alchElecLambdaStart 0.5` and `alchVdwLambdaEnd 1.0`. Nothing was tuned. See
    `6I5I_DUAL_FEP/DDG_preliminary.md` §7.2 for what the flag actually controls
    versus what the original rationale claimed.
- Previous: README updates, visualization, CHARMM-GUI inputs

## 6I5I_DUAL_FEP — dual-topology RBFE prep (no FEPrepare, no RDKit)

`6I5I_DUAL_FEP/` is a self-contained, from-scratch prep for the **6I5I H3E ("12H")
→ desmethyl** perturbation. Reference ligand carries an **N–CH₃**, the mutant an
**N–H**. Uses NAMD's default **dual topology** (`singleTopology off`): the two
ligands share a common core (present once), and the atoms that differ are BOTH
present and switched on/off through the `alchFile` **B-factor** column:
**B=−1** methyl (vanish at λ=1), **B=0** common, **B=+1** N–H (appear at λ=1).

Three plain-Python scripts (numpy only, no FEPrepare/RDKit):
1. `prepare_hybrid.py` — maps the ligands chemically (identifies methyl/N–H by
   local environment, NOT graph-isomorphism which fails on the symmetric ring),
   places the N–H by local N-frame Kabsch, writes `hybrid.rtf/prm/pdb`.
2. `build_system.py` — VMD psfgen (protein split at the 412→416 gap into P1/P2
   + ligand) → solvate → autoionize, for complex + solvent.
3. `write_fep_inputs.py` — marks `ionized.fep` B-factors, writes `fep.tcl`
   (16-λ uneven schedule) + `nvt/npt/md_forward/md_backward.namd`.

Workflow: `python3 prepare_hybrid.py && python3 build_system.py &&
python3 write_fep_inputs.py`, then `bash run_all.sh`.

### Gotchas learned (6I5I prep)
- **Restart-chained `.namd` stages must NOT set `temperature`.** `npt_equil` /
  `md_forward` / `md_backward` read a prior stage's `binvelocities`; adding
  `temperature $temp` in the same config is fatal: `FATAL ERROR: Cannot specify
  both an initial temperature and a velocity file`. Only `nvt_equil` (which
  seeds via `reinitvels $temp`) sets `temperature`. This was a bug in
  `write_fep_inputs.py` (emitted `temperature` in the common header for every
  stage) — fixed 2026-09-05; keep it out of the shared header.
- Water/ions params: use `toppar/par_water_ions_clean.prm` (NAMD cannot parse
  `toppar_water_ions.str`, a topology+param stream file).
- acpype CHARMM rtf uses `IMPH` (not `IMPR`) for impropers.
- CHARMM prm `X` wildcard must stay uppercase (don't `.lower()` it).
- The N–H improper must mirror the mutant's `C H N N1` = types `cc hn na nc`.
- `runFEP <l0> <l1> <dl>` with `dl=0` divides by zero — guard with `dl>0`.
- Protein chain A has a disordered gap 412→416 (two segments P1/P2).

## ✅ NAMD GPU-resident — SOLVED (2026-08-25)

The V100 slowness was a **missing config flag**, not a bad Charm++ arch.

`GPUresident`/`CUDASOAintegrate` is gated at compile time by `NODEGROUP_FORCE_REGISTER`
(`src/SimParameters.C:5008`, inside an `#ifndef`). That macro is emitted **only** by
`./config --with-single-node-cuda` (`config:1052-1053`). Plain `--with-cuda` never
defines it — hence `FATAL ERROR: GPUresident not supported on regular multicore builds`.

An earlier note in this file suggested rebuilding Charm++ as `verbs`/`mpi`. **That was
wrong.** `multicore-linux-x86_64` is the *correct* arch for single-node GPU-resident;
the existing Charm++ build was reused untouched.

### GPU-resident build
```bash
./config Linux-x86_64-g++.gpuresident \
  --charm-arch multicore-linux-x86_64 \
  --with-single-node-cuda \
  --cuda-prefix /usr/local/cuda-11.8 \
  --cuda-gencode arch=compute_70,code=sm_70
make -j32   # ~13 min on 64 cores
```
**Use this binary for production:**
`/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3`
(the old `Linux-x86_64-g++/namd3` is kept for comparison).

Add `GPUresident on` to configs (`CUDASOAintegrate` is the deprecated spelling).

### 🔴 CRITICAL: run with `+p1`, NOT `+p8`

In GPU-resident mode virtually all work runs on the GPU, so **extra CPU threads only
add synchronization contention**. Throughput is almost exactly *inversely* proportional
to PE count. Measured (64651-atom complex, plain MD, same config):

| threads | ns/day |
|---------|--------|
| **`+p1`** | **75.79** ← use this |
| `+p2` | 59.54 |
| `+p4` | 27.22 |
| `+p8` | 13.22 |
| `+p16` | 6.29 |
| `+p32` | 2.70 |

Running `+p8` (the old habit, correct for the *non*-GPU-resident build) costs ~5.7×.
Counter-intuitive but reproducible: **one PE, one GPU.**

`stepspercycle 400` / `pairlistsPerCycle 40` / `margin 4.0` made **no** difference
(75.66 vs 75.79 ns/day at `+p1` — noise). Don't bother; PE count is the whole story.

### Measured (64651-atom 6I5I complex, identical configs)
| run | old build `+p8` | GPU-resident `+p8` | GPU-resident `+p1` | total speedup |
|-----|-----------------|--------------------|--------------------|---------------|
| plain MD | 3.57 ns/day | 14.08 ns/day | **75.79 ns/day** | **21×** |
| FEP (1 window) | 2.26 ns/day | 8.73 ns/day | **49.65 ns/day** | **22×** |

Against the originally recorded baseline (~1.7 / ~1.2 ns/day) this is ~45× / ~41×.

### Production speeds (measured on the *generated* configs, not synthetic benchmarks)

Steady-state, `alch on`, minimization excluded, `TIMING:` sec/step:

| system | atoms | rate | sec/step |
|--------|-------|------|----------|
| complex (`6I5I_DUAL_FEP/complex`) | 64,651 | **47 ns/day** | 0.00367 |
| solvent (`6I5I_DUAL_FEP/solvent`) | 5,013 | **59 ns/day** | 0.00293 |

The solvent leg is only 13× smaller but latency-bound on the V100, so it's not 13×
faster — it saturates around ~59 ns/day regardless of size at `+p1`.

Per leg: 15 windows × 250,000 steps = 3.75M steps per direction (7.5 ns) + ~100k steps
equilibration. Projected wall-clock (sequential, `run_all.sh`):

| stage | complex | solvent |
|-------|---------|---------|
| equil (nvt+npt) | ~6 min | ~5 min |
| FEP forward | ~3.8 h | ~3.1 h |
| FEP backward | ~3.8 h | ~3.1 h |
| **leg total** | **~7.7 h** | **~6.2 h** |

**Full job (complex + solvent) ≈ 14–16 hours**, vs ~25–30 days at the original
baseline. The small system is not faster per-atom; it's latency-limited.

### GPU memory is NOT a constraint
NAMD peaks at **~590 MiB** on this system. The card reports 16384 MiB available
(a 32 GB V100 SKU, partitioned). **Moving to a larger-memory GPU would change nothing** —
the limit is compute/latency, not capacity. Only a faster GPU (A100/H100) or a smaller
PE count helps.

**Accuracy validated:** single-point energy on identical coordinates agrees to ~8
significant figures (POTENTIAL −199582.4277 vs −199582.4190; rel. diff 4e-8, i.e.
float32 rounding). The speedup is not an accuracy tradeoff.

### Caveats
- **Don't misread `PERFORMANCE:` lines.** That number is a *cumulative* average from
  step 0, so a leading `minimize 5000` drags it down for a long time — `nvt_equil`
  shows ~3.9 ns/day during minimization and ~23 ns/day cumulative well into MD, while
  the *instantaneous* rate (`TIMING:` sec/step) is **49.8 ns/day**. Judge throughput
  from `TIMING:` sec/step, not the cumulative average.
- **`GPUresident` must be set in exactly one place.** Setting it both in the config and
  as `--GPUresident on` on the command line is fatal:
  `ERROR: Multiple definitions of 'GPUresident'`. The 6I5I configs carry it *in the
  config* (so `run_all.sh` passes only `+p1 +devices 0`); `fep_pipeline` instead passes
  it *on the command line* (so the `--no-gpu` CPU path still works). Don't mix.
- Output frequency is **not** a bottleneck here — everything at 500
  (`restartfreq`/`dcdfreq`/`XSTFreq`/`alchOutFreq`) measured identically to 5000.
- Logs `Disabling lonepair support due to incompatability with GPU-resident`.
  Harmless **for this system** — `ionized.psf` has no `NUMLP` section (0 LP atoms).
  Re-check if a future ligand uses CGenFF halogen lone pairs.
- Raw `ionized.pdb` has severe clashes (VDW ~2.6e6) — it must be minimized before any
  restart-free start, or RATTLE fails on atom 2519.
