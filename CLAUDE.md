# NAMD-FEP Project

## Overview
A tutorial repo for Relative Binding Free Energy (RBFE) calculations using NAMD.
Calculates ΔΔG of binding for small molecule ligands against a protein target.

## fep_pipeline — Self-Contained RBFE Preparation Toolkit

`fep_pipeline/` is a Python 3 package that replaces the Feprepare web server,
Maestro, and LigParGen with a single local workflow. It generates everything
needed to run an RBFE simulation with NAMD + VMD + Python.

**Key modules:**
| Module | Purpose |
|--------|---------|
| `alignment.py` | Ligand alignment (Kabsch), hybrid PDB, FEP file generation |
| `namd_config.py` | NAMD config generation, lambda schedules, GPU command builder |
| `system_builder.py` | VMD Tcl script generation for psfgen/solvate/ionize |
| `fep_file.py` | FEP file builder (B-factor column marking λ-dependence) |
| `analysis.py` | Parse .fepout files, BAR analysis, ΔΔG computation |
| `ligand_param.py` | Ligand parameterization (OpenMM/GAFF or Antechamber) |
| `prepare_fep.py` | Main orchestrator + CLI |

**Usage:**
```bash
# Full pipeline (needs VMD installed)
python -m fep_pipeline run \
    --protein protein.pdb \
    --ligand-a ref.pdb --ligand-b mut.pdb \
    --workdir ./fep_run

# Align ligands only
python -m fep_pipeline align --reference ref.pdb --mobile mut.pdb --output aligned.pdb

# Hybrid PDB + FEP file only
python -m fep_pipeline hybrid --ligand-a ref.pdb --ligand-b mut.pdb \
    --output hybrid.pdb --fep-output ligand.fep

# Analyze results
python -m fep_pipeline.analysis complex/md_forward.fepout complex/md_backward.fepout \
    --solvent-forward solvent/md_forward.fepout \
    --solvent-backward solvent/md_backward.fepout
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
| NAMD binary | `/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++/namd3` |
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
- `089aabb` — Fix NAMD FEP configs: alchDecouple off for RBFE, non-zero alchEquilSteps for ParseFEP, shell portability (2026-07-29)
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
- Water/ions params: use `toppar/par_water_ions_clean.prm` (NAMD cannot parse
  `toppar_water_ions.str`, a topology+param stream file).
- acpype CHARMM rtf uses `IMPH` (not `IMPR`) for impropers.
- CHARMM prm `X` wildcard must stay uppercase (don't `.lower()` it).
- The N–H improper must mirror the mutant's `C H N N1` = types `cc hn na nc`.
- `runFEP <l0> <l1> <dl>` with `dl=0` divides by zero — guard with `dl>0`.
- Protein chain A has a disordered gap 412→416 (two segments P1/P2).

## ⚠️ NAMD GPU offload is BROKEN/slow on this box (unresolved)

The V100 is **not** actually accelerating NAMD3. `namd3 +p8 +devices 0` binds the
CUDA device but the GPU sits at **~2% utilization / 534 MB**, and the 64651-atom
complex runs at **~1.7 ns/day (plain MD) / ~1.2 ns/day (FEP)** — ~15–25× too slow
(a V100 should do 20–40 ns/day here). CPU is ~9% busy, so it's a serialization /
offload problem, not compute. Clues:
- `FATAL ERROR: GPUresident not supported on regular multicore builds` when using
  `--CUDASOAintegrate on`.
- `Warning: Always using force tables for GPU nonbonded kernel due to unsupported
  config parameters` (ruled out NBFIX and switching as causes — both still ~1.7 ns/day).

**TODO for next session:** diagnose the GPU offload. Likely needs a rebuild of NAMD
with proper GPU-resident support (the current build is `multicore-linux-x86_64`
Charm++ + `--with-cuda`, which apparently lacks effective GPU offload), or a
different Charm++ arch (e.g. `verbs`/`mpi`) with `--with-cuda`. Until fixed, a
16-window production run is impractical (~25–30 days).
