# 6I5I dual-topology RBFE — preparation (from scratch, no FEPrepare)

Relative binding free energy for the **H3E (12H) → desmethyl** perturbation:
reference ligand carries an **N–CH₃**, the mutant an **N–H**.  In NAMD's default
**dual topology** scheme (`singleTopology off`) the two ligands share a common
core (present once) while the atoms that differ are both present and are
switched on/off through the `alchFile` **B-factor** column:

| B-factor | meaning                                  |
|----------|------------------------------------------|
| −1       | reference-only (the methyl C12 + H7/H8/H9) → vanishes at λ=1 |
|  0       | common atoms (identical in both states) → always on |
| +1       | mutant-only (the N–H, H17) → appears at λ=1 |

## What's here

```
6I5I_DUAL_FEP/
├── inputs/           raw inputs (GAFF/acpype CHARMM rtf/prm/pdb + protein.pdb)
├── hybrid/           the merged dual-topology ligand
│   ├── hybrid.pdb    45 atoms, B-factor column marks λ-dependence
│   ├── hybrid.rtf    merged CHARMM topology (common + methyl + N-H)
│   └── hybrid.prm    merged parameters (ref params + the N-H terms)
├── complex/          protein + ligand, solvated + ionized
│   ├── ionized.psf / ionized.pdb / ionized.fep
│   ├── nvt_equil.namd / npt_equil.namd / md_forward.namd / md_backward.namd
├── solvent/          ligand alone, solvated + ionized  (same files)
├── fep.tcl           runFEP / runFEPmin procs (sourced as ../fep.tcl)
├── prepare_hybrid.py builds the hybrid ligand (atom mapping + B-factors + params)
├── build_system.py   VMD psfgen + solvate + ionize (complex + solvent)
└── write_fep_inputs.py  .fep files + fep.tcl + NAMD configs
```

## Workflow (regenerate everything)

```bash
cd 6I5I_DUAL_FEP
python3 prepare_hybrid.py      # 1. hybrid ligand (B-factors + topology + params)
python3 build_system.py        # 2. VMD: psfgen -> solvate -> ionize
python3 write_fep_inputs.py    # 3. ionized.fep + fep.tcl + *.namd
```

## Run

NAMD3 binary: `/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++/namd3`

```bash
# each leg, in order (nvt -> npt -> forward -> backward)
cd complex
namd3 +p8 +setcpuaffinity --CUDASOAintegrate on +devices 0 nvt_equil.namd
namd3 +p8 +setcpuaffinity --CUDASOAintegrate on +devices 0 npt_equil.namd
namd3 +p8 +setcpuaffinity --CUDASOAintegrate on +devices 0 md_forward.namd
namd3 +p8 +setcpuaffinity --CUDASOAintegrate on +devices 0 md_backward.namd
cd ../solvent   # repeat the same four
```

Each leg runs 16 λ windows through `runFEP` (uneven spacing near the endpoints
is recommended for the electrostatics decoupling; the current `md_*.namd` uses
even 1/15 spacing — adjust `runFEP` args if you want the uneven schedule).

## Analysis

```bash
# in VMD: Extensions -> Analysis -> FEP (ParseFEP)
# forward + backward .fepout for complex, then for solvent
# ΔΔG = ΔG(complex) − ΔG(solvent)
```

## Notes / assumptions

- The two ligands differ only by the methyl↔H on the shared ring nitrogen
  (verified by atom-type matching: `na↔na`, `nc↔nc`; 40 common / 4 vanish / 1 appear).
- Common atoms use the **reference** conformation, types and charges; the N–H is
  placed by a local N-frame rigid alignment (N–H bond = 1.01 Å).
- `par_water_ions_clean.prm` (not `toppar_water_ions.str`) is used for
  water/ions — NAMD cannot parse the `.str` stream file directly.
- One cross angle (`c3 na hn`) has no source parameter (the methyl and N–H are
  never simultaneously present); it is filled with an approximate value.
- Protein is chain A with a disordered gap 412→416, split into segments P1/P2.
