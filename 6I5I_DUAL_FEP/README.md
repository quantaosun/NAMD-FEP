# 6I5I dual-topology RBFE — worked example

The concrete system the [root README](../README.md) walks through. Read that
first: it has the general workflow, the prerequisites, the `+p1` rule, the
commands and the analysis. This file records what is specific to **6I5I**.

## The perturbation

PDB **6I5I**, CLK1 kinase, ligand **H3E ("12H")**. The reference ligand carries
an **N–CH₃** on a pyrazole nitrogen; the mutant is the **desmethyl** analogue
with an **N–H** in its place.

In NAMD's default **dual topology** scheme (`singleTopology off`) the two
ligands share a common core (present once) while the atoms that differ are both
present and are switched on/off through the `alchFile` **B-factor** column:

| B-factor | meaning                                                      |
|---------:|--------------------------------------------------------------|
| −1       | reference-only (the methyl C12 + H7/H8/H9) → vanishes at λ = 1 |
|  0       | common atoms (identical in both states) → always on           |
| +1       | mutant-only (the N–H, H17) → appears at λ = 1                 |

Atom mapping verified by atom-type matching (`na↔na`, `nc↔nc`):
**40 common / 4 vanish / 1 appear**.

## Layout

```
6I5I_DUAL_FEP/
├── inputs/              raw inputs
│   ├── protein.pdb        chain A, disordered gap 412→416, no hydrogens
│   ├── ref.{mol2,pdb,rtf,prm}   reference ligand (N–CH₃), acpype/GAFF→CHARMM
│   └── mut.{mol2,pdb,rtf,prm}   mutant ligand (N–H)
├── hybrid/              the merged dual-topology ligand (generated)
│   ├── hybrid.pdb        45 atoms; B-factor column marks λ-dependence
│   ├── hybrid.rtf        merged CHARMM topology (common + methyl + N–H)
│   └── hybrid.prm        merged parameters (ref params + the N–H terms)
├── complex/             protein + ligand, solvated + ionized
├── solvent/             ligand alone, solvated + ionized
├── fep.tcl              runFEP / runFEPmin procs (sourced as ../fep.tcl)
├── prepare_hybrid.py    ① hybrid ligand (mapping + B-factors + params)
├── build_system.py      ② VMD psfgen → solvate → ionize (complex + solvent)
├── write_fep_inputs.py  ③ ionized.fep + fep.tcl + *.namd
├── fep_run.py           per-λ-window NAMD driver
├── run_checkpointed.sh  crash-safe controller — the way to run
├── respawn_controller.sh, watch_windows.sh, checkpoint_status.sh
├── run_all.sh           original one-shot launcher (fresh runs only)
├── analyze_fep.py       EXP / BAR estimators
└── audit_fep.py         independent audit + error bars — what to trust
```

## Regenerating from scratch

```bash
cd 6I5I_DUAL_FEP
python3 prepare_hybrid.py      # 1. hybrid ligand
python3 build_system.py        # 2. VMD: psfgen → solvate → ionize
python3 write_fep_inputs.py    # 3. ionized.fep + fep.tcl + *.namd
```

Then run and analyse as described in the [root README §3–§4](../README.md#3-running-the-production-job).

## Result

```
dG_complex = +3.882 kcal/mol
dG_solvent = +3.988 kcal/mol
─────────────────────────────────
ΔΔG_bind   = −0.106 kcal/mol
95% CI (moving-block bootstrap) = [−0.355, +0.118]
```

**Indistinguishable from zero.** The [root README §5](../README.md#5-worked-example--6i5i-clk1)
explains why this number is fragile — the near-cancellation of two ~2.3 kcal/mol
halves that flip sign at λ ≈ 0.5, the non-stationary complex leg, and the
sensitivity to the equilibration trim. Do not quote it as a binding difference.

## Notes and assumptions

- Common atoms use the **reference** conformation, types and charges; the N–H is
  placed by a local N-frame rigid alignment (N–H bond = 1.01 Å).
- `par_water_ions_clean.prm` (not `toppar_water_ions.str`) is used for water and
  ions — NAMD cannot parse the `.str` stream file directly.
- One cross angle (`c3 na hn`) has no source parameter, because the methyl and
  the N–H are never simultaneously present. It is filled with an approximate
  value.
- The protein is chain A with a disordered gap 412→416, split into segments
  P1/P2 by `build_system.py`.
- The λ schedule in `fep.tcl` is the **uneven** 16-value one quoted in the root
  README, not an even 1/15 spacing.

Full theory and process notes: [DDG_preliminary.md](DDG_preliminary.md).
