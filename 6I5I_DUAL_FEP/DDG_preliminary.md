# 6I5I ΔΔG — preliminary analysis notes

**Date:** 2026-09-08 · **Status:** PRELIMINARY — one-sided EXP preview; final BAR pending
the solvent backward leg (~01:00 CST 2026-09-09).

---

## 1. The modification (before → after)

Relative binding free energy between two ligands of **CLK1** (PDB entry 6I5I):

| | state | name | SMILES (canonical) | formula | MW |
|---|---|---|---|---|---|
| **Before** (λ = 0) | REF | **12H** (PDB ligand **H3E**) | `Cn1cc(-c2ccc3occ(-c4cnn(Cc5ccccc5)c4)c3n2)cn1` | C₂₁H₁₇N₅O | 355.40 |
| **After**  (λ = 1) | MUT | **desmethyl** | `c1ccc(Cn2cc(-c3coc4ccc(-c5cn[nH]c5)nc34)cn2)cc1` | C₂₀H₁₅N₅O | 341.37 |

**The change is a single demethylation of a pyrazole-ring nitrogen:**
`aromatic N–CH₃  →  N–H`  (loss of CH₂ = −14 Da, 27 → 26 heavy atoms).

Verified with RDKit: the maximum common substructure of the two ligands is the whole
mutant (26 heavy atoms), i.e. **nothing differs except that one methyl group**.

> SMILES in `6I5I_FEP/1_ligand/ref_crystal.sdf` (REF) and `.../mut_aligned.sdf` (MUT),
> generated from those structures with RDKit 2026.03.4. Both neutral overall.

### Alchemical meaning
In the dual-topology setup the two end states are both "present" but switched by the
B-factor column of the FEP file: **B = −1** marks the methyl (it *vanishes* as λ→1),
**B = 0** the common core, **B = +1** the N–H (it *appears* as λ→1). So

- **λ = 0** ⇔ REF (12H / N–CH₃) fully interacting
- **λ = 1** ⇔ MUT (desmethyl / N–H) fully interacting
- the **forward** run travels 0 → 1, i.e. REF → MUT in the given environment.

---

## 2. Theory: why ΔΔG = ΔG_complex − ΔG_solvent

Absolute binding free energies are expensive. RBFE instead closes a **thermodynamic
cycle** and measures only *relative* changes:

```
        complex·REF ──ΔG_bind(REF)──▶ protein + REF(solv)
            │                            │
 ΔG_complex │  (mutate ligand           │ ΔG_solvent
  (in site) │   in the binding site)    │  (in water)
            ▼                            ▼
        complex·MUT ──ΔG_bind(MUT)──▶ protein + MUT(solv)
```

By path-independence of free energy,

```
ΔΔG = ΔG_bind(MUT) − ΔG_bind(REF) = ΔG_complex − ΔG_solvent
```

where the two **alchemical legs** mutate the ligand (REF → MUT) either bound to the
protein (`complex`) or free in water (`solvent`). **ΔΔG < 0 ⇒ the mutant binds more
tightly** (more negative ΔG_bind). Environment-specific (protein) interactions cancel
to zero in this difference; only the *change in ligand–environment coupling upon
demethylation* survives.

### Why two "physical states"? (forward/backward)
ΔG of a leg is computed by *sampling*, which is stochastic and biased if done only one
way. The rigorous estimator, **BAR**, needs the difference-energy ΔU sampled from **both**
end states of every λ-window — hence each leg is run twice: `md_forward` (0→1) and
`md_backward` (1→0). The two directions also give the **hysteresis** = |ΔG_fwd + ΔG_bwd|,
a convergence/quality check (0 would be perfect reversibility).

---

## 3. Dual-topology FEP in NAMD

- `singleTopology off` (NAMD default) — both ligands coexist on a shared core; the atoms
  that differ are BOTH present and turned on/off by the `alchFile` B-factor column.
- `alchDecouple off` — appropriate for RBFE (electrostatics of the changed group are
  switched together with its vdW, preserving net charge evolution along λ).
- 15 λ-windows on an **uneven** schedule (denser near the end states), one NAMD process
  per window, each window = 50 000 equilibration steps + 250 000 production steps.
- At every 500th step NAMD evaluates the hybrid Hamiltonian at both λ and λ+dl and writes
  a `FepEnergy:` row whose column 6 is the instantaneous **dE = U(λ+dl) − U(λ)** — the raw
  material for every estimator below. Per window: 500 independent-of-output samples.

---

## 4. Estimators

### One-sided exponential averaging (Zwanzig / EXP) — *used for the preview*
For a window spanning λ → λ+dl, sampling at the base λ:

```
ΔG_i = −(1/β) · ln ⟨ exp(−β · dE) ⟩ ,        β = 1/(k_B T)
```
Implement with the log-sum-exp trick for numerical stability. Sum over windows →
total leg ΔG. **Uses forward data only**, which is exactly why the preview is possible
while solvent backward is still running. Caveat: EXP is a *biased* (and noisier)
estimator — large, poorly-sampled dE tails bias it; the bias partly cancels in the
complex−solvent difference but is not guaranteed to.

### Two-sided Bennett (BAR) — *final number, pending*
Solves for the window ΔG such that the acceptance ratio of forward and backward ΔU
samples is balanced (implemented in pure Python in `analyze_fep.py`, no pymbar).
Unbiased to first order and the standard final analysis.

---

## 5. How the current ΔΔG was obtained (process)

1. **Production sampling** — four legs, marker-guarded, GPU-resident NAMD3 `+p1` on the
   V100: `complex` fwd/bwd, `solvent` fwd/bwd. 15 windows × 250 000 steps each.
   (Solvent backward was still running when this preview was made; complex fwd/bwd and
   solvent fwd were complete.)
2. **Assemble** — `fep_run.py assemble` builds one canonical `<stage>_combined.fepout`
   per leg (500 FepEnergy rows/window, windows delimited by the λ schedule).
3. **Parse** — `analyze_fep.py exp` reads the **dE column** of each forward combined
   fepout, does per-window EXP, and sums.
4. **Result (forward-only EXP preview, T = 300 K, k_BT = 0.596 kcal/mol):**

| leg (0→1) | per-window EXP sum |
|---|---|
| `complex/md_forward_combined.fepout` | **+4.061** kcal/mol |
| `solvent/md_forward_combined.fepout` | **+3.969** kcal/mol |
| **ΔΔG = complex − solvent** | **≈ +0.09 kcal/mol** |

    → the desmethyl mutant and 12H bind **essentially equally** (ΔΔG ≈ 0.1 kcal/mol).

5. **Convergence cross-check (complex leg already has both directions):**
   EXP over the complex *backward* file = **−3.930** kcal/mol vs forward **+4.061** ⇒
   one-sided **hysteresis ≈ 0.13 kcal/mol** — the leg is well converged, so the near-zero
   ΔΔG is not an artifact of one direction.

6. **Pending:** when solvent backward completes, run the unbiased **BAR** over both legs
   (2-state Bennett per window) and report ΔΔG with the hysteresis for both legs.

Command (once all legs are done):
```bash
python3 analyze_fep.py bar \
  complex/md_forward_combined.fepout  complex/md_backward_combined.fepout \
  solvent/md_forward_combined.fepout  solvent/md_backward_combined.fepout
```

---

## 6. Caveats (preliminary)

- Numbers above are **one-sided EXP**, a biased estimator. Treat as a preview/ordering,
  not the publication value.
- No error bars yet (EXP standard error needs block analysis; BAR gives none directly —
  bootstrap/block would follow).
- dG column printed by NAMD is *not* used; we recompute from raw dE.
- Readout is for the current, mostly-complete trajectory; if solvent backward's windows
  disagree with forward, the EXP preview and final BAR will differ visibly.

---

## 7. Corrections after the 2026-09-12 audit (`audit_fep.py`)

The sections above are kept as written (they document the process). Three things in them
are now known to be wrong or unsupported.

### 7.1 The reported ΔΔG was computed over the wrong samples

`analyze_fep.py` averages **all 500** `FepEnergy` rows of every window. Only **401** are
production data. With `alchEquilSteps 50000` at 2 fs (= 100 ps) and `alchOutFreq 500`
(= 1 ps), the first **99 rows of every window** (step 500 … 49500) are the per-window
equilibration that `alchEquilSteps` exists to discard, and NAMD's own accumulator resets at
`stepInRun == alchEquilSteps` (`Controller::outputFepEnergy` in `src/Controller.C`), so the
`dE_avg` and `dG` columns are built from the 401 production rows only.

This is checkable, not a judgement call: `mean(rows[99:])` reproduces NAMD's own final
`dE_avg` to 4 decimals, and one-sided EXP over `rows[99:]` reproduces NAMD's own per-window
`dG` (the `#Free energy change for lambda window [...] is <dG>` comment) to **7e-6 kcal/mol
across all 30 forward windows**.

| | dG_complex | dG_solvent | **ΔΔG** |
|---|---|---|---|
| all 500 rows (as reported in §5) | +3.917 | +3.974 | **−0.057** |
| 401 production rows (correct) | +3.882 | +3.988 | **−0.106** |

So the headline number moves by ~0.05 kcal/mol, and — more importantly — it was never more
precise than that. Block bootstrap (25 ps blocks, 400 resamples) on the corrected number:
**ΔΔG = −0.106 ± 0.120 kcal/mol, 95% CI [−0.355, +0.118]** → consistent with zero. The
qualitative conclusion of §5 ("mutant and 12H bind essentially equally") survives; the
precision implied by "−0.057" does not.

### 7.2 `alchDecouple off` was chosen for a reason that does not apply to it

§3 justifies `alchDecouple off` as "electrostatics of the changed group are switched together
with its vdW, preserving net charge evolution along λ". That describes **`alchElecLambdaStart`
/ `alchVdwLambdaEnd`**, not `alchDecouple`. `alchDecouple` is a separate flag whose ON branch
is described in the source (`src/ComputeNonbondedCUDAExcl.C:102-107`) as changing how PME
handles the alchemical partition — extra grids so that the alchemical atoms' reciprocal-space
interaction with the full system is switched off while a grid containing only alchemical atoms
is switched on.

What is certain: **`off` is NAMD's default** (`SimParameters.C:1250-1251`,
`&alchDecouple, FALSE`), as are `alchElecLambdaStart 0.5` and `alchVdwLambdaEnd 1.0`
(`SimParameters.C:1233-1239`). So the alchemical protocol here is unmodified NAMD defaults —
nothing was actually tuned, and the config comment records a rationale for something it does
not do. Note also that `alchElecLambdaStart (0.5) < alchVdwLambdaEnd (1.0)`, i.e. the
appearing N–H's charge starts turning on before its vdW is fully on; the soft-core
(`alchVdWShiftCoeff 5.0`, a real option — `ComputeNonbondedUtil.h:390`, applied as
`alchVdwShiftCoeff*(1-vdwLambda)` at `ComputeNonbondedBase.h:408`) is what keeps that from
being an endpoint catastrophe. Whether `off` is the better choice for this specific
methyl→H perturbation is **untested here** — it mostly cancels in the complex−solvent
difference, but "mostly cancels" is an assumption, not a measurement.

### 7.3 §6's "no error bars yet" — now there are, and they dominate

| diagnostic | complex | solvent |
|---|---|---|
| EXP hysteresis (fwd + bwd) | +0.059 | −0.092 |
| Σ(second half − first half) per window | **+0.365** | −0.007 |

The complex leg is **not stationary**: it is still drifting by +0.37 kcal/mol between the
first and second half of its own production windows (largest single window: w11, +0.258).
The solvent leg is flat. Two things plausibly cause this and neither is a code bug:
**backward w00 restarts from `npt_equil.coor` at λ = 1.0** instead of from forward w14's
endpoint (so its environment is equilibrated around the methylated ligand while being sampled
at the desmethyl state), and 500 ps per window is short for a protein site.

### 7.4 A structural fragility worth knowing about

ΔΔG is a near-cancellation of two large opposing halves. Summing the per-window
`dG_complex − dG_solvent` differences, the λ < 0.5 half contributes ≈ **−2.37** kcal/mol and
the λ > 0.5 half ≈ **+2.26**, giving ≈ −0.1. The sign flips at λ ≈ 0.5, which is exactly
`alchElecLambdaStart` — the crossover where the vanishing methyl's charge finishes going off
and the appearing N–H's starts coming on. So the answer is the residue of two ~2.3 kcal/mol
terms that nearly cancel; any systematic error in either half lands almost undiluted on ΔΔG.
That, not the estimator, is why this quantity is hard.

**Practical implication.** The honest error bar is ≈ ±0.2 kcal/mol, not ±0.06. Improving it
means more sampling per window (2–5 ns), more windows around λ 0.5–0.8, and seeding the
backward leg from the forward endpoint — not a different estimator.

### 7.5 Known limitation of the existing result: the hybrid charge cycle did not close

`hybrid/hybrid.rtf` sums to **+0.174789** while both `inputs/ref.rtf` and `inputs/mut.rtf`
sum to 0. The shared core keeps the *reference* charges, but ref and mut disagree on the
anchor nitrogen (ref `N4` = −0.235229, mut `N` = −0.245202). `prepare_hybrid.py` copied the
appearing H's charge verbatim out of `mut.rtf` (+0.174786), so nothing absorbed that
difference and:

| state | net charge |
|---|---|
| λ = 0 (core + methyl) | +0.000003 ✓ |
| λ = 1 (core + N–H) | **+0.0162** ✗ |

**So the λ = 1 endpoint is not exactly the neutral desmethyl ligand** — it carries +0.0162 e.
The magnitude is small and the same defect is present in both the complex and solvent legs,
so it largely cancels in ΔΔG; it does not change the qualitative conclusion in §7.1. It is
recorded here rather than silently fixed because fixing it changes the topology and would
require re-running both legs (~14 h) to produce a paired number.

`prepare_hybrid.py` **is now fixed** for future systems: it computes the required
appearing-group total as `mut_total − core_total`, applies the correction, and raises if the
endpoints do not close. Re-running it on the 6I5I inputs now gives λ=0 +0.000003 / λ=1
+0.000000, and changes **exactly two lines** versus the committed file —
`ATOM H17 0.174786 → 0.158610` and `RESI UNL 0.000 → 0.158613`. `hybrid.prm` and
`hybrid.pdb` are byte-identical.

⚠️ **The committed `hybrid/` was deliberately NOT regenerated.** It is the topology the
production run actually used, so regenerating it would desync the repo from the number in
§7.1. If you ever re-run 6I5I, regenerate first and note that the new result will differ
slightly from −0.106 for this reason as well as for sampling.
