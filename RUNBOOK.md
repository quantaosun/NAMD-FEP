# RUNBOOK — relative binding free energy, one complex at a time

The commands are the **same six for every system**. Everything that differs
between one protein-ligand complex and the next lives in that system's
`system.ini` — not in the code, and not in a directory you copy and hand-edit.

```
rbfe hybrid → rbfe build → rbfe inputs → rbfe check → rbfe run → rbfe audit
```

`README.md` is the theory and the background. This file is the procedure.

---

## 0. What this costs

Per system, on one V100 with `+p1`: **~15–31 h**, depending on the number of
windows and the size of the box. The two worked systems bracket it — 6I5I is
15 ns/leg and took ~14–16 h; 4YLJ is 15 ns/leg on a larger system at 39 ns/day
and projects ~31 h. `rbfe check` measures *your* system rather than guessing.

It is resumable at window granularity. A crash costs at most the window in
flight, so a machine that reboots every few hours is survivable (see §5).

---

## 1. One-time setup

```bash
cd /path/to/NAMD-FEP
export PATH="$PWD/bin:$PATH"          # once per shell; put it in your profile
```

Then check the three things the pipeline actually needs:

```bash
"$NAMD" --version                  # NAMD 3.x, GPU-resident build
"$VMD"  -dispdev text -e /dev/null # VMD with psfgen + solvate + autoionize
python3 -c "import numpy; print(numpy.__version__)"
```

**The pipeline's dependencies are Python 3.8+ and numpy — nothing else.** No
FEPrepare, no Maestro, no LigParGen, no notebook. The config format is INI
precisely so this stays true (`tomllib` needs 3.11; YAML would add a package).

Two optional extras, both degrading rather than failing when absent:

| for | needs | without it |
|---|---|---|
| `rbfe-prep` (§2 — making inputs from a PDB) | `acpype` and `obabel` | you prepare ligands by hand |
| `rbfe similarity` and `rbfe-prep ligand`'s formula check | `rdkit` | the screen says `similarity NOT checked` on stderr, never a silent pass |

Binary paths: `system.ini` defaults them to `$NAMD` / `$VMD` if set, and to the
paths on this box otherwise. Nothing else in the config is machine-specific.

---

## 2. Prepare the inputs

Before step ①, a system needs:

| file | what |
|---|---|
| `inputs/protein.pdb` | the protein. Hydrogens optional — psfgen adds them. One chain. |
| `inputs/ref.{pdb,rtf,prm,mol2}` | the reference ligand, **all-atom**, with CHARMM parameters |
| `inputs/mut.{pdb,rtf,prm,mol2}` | the mutant ligand, same |

`rbfe-prep` makes them from a PDB entry. It does the file handling; it does
**not** choose the mutation, because that is a chemistry decision.

```bash
rbfe-prep fetch   3HTB                    # or use a PDB you already have
rbfe-prep report  3HTB.pdb                # read this before choosing anything
rbfe-prep strip   3HTB.pdb -o systems/7abc/inputs/protein.pdb --chain A
rbfe-prep ligand  3HTB.pdb JZ4 -o systems/7abc/inputs/mut.sdf
rbfe-prep ligand  3HTB.pdb JZ4 -o systems/7abc/inputs/ref.sdf --drop C4
rbfe-prep params  systems/7abc/inputs/ref.sdf --name ref -o systems/7abc/inputs/
rbfe-prep params  systems/7abc/inputs/mut.sdf --name mut -o systems/7abc/inputs/
rbfe-prep scaffold systems/7abc --chain A --resname UNL
```

`report` is read-only and tells you which chains and HETATM species exist,
whether the protein has numbering gaps, how many altlocs, and whether any
modified residue will not build. `--drop` deletes atoms before the hydrogens are
added, which is how you derive the smaller ligand from the crystal one while
keeping the crystallographic pose.

**Check what it printed, every time.** `ligand` prints the formula and SMILES it
perceived, and `params` prints the ligand's total charge. These are perceived,
not known, and acpype's charges are computed from whatever was perceived. If the
SMILES is wrong, pass `--smiles` with the intended molecule and the crystal
frame is kept. If the charge is not within ~0.02 of an integer, `rbfe hybrid`
will refuse the ligand.

Requirements: `acpype` and `obabel` (AmberTools' `antechamber`/`sqm` come bundled
inside acpype). On this box they are in `~/external-libraries/bin`, which is not
on PATH; `rbfe-prep` finds them and sets PATH for the subprocess. Override with
`$ACPYPE` / `$OBABEL`.

**All-atom is not optional.** The hybrid builder matches the two ligands by atom
name (element swap) or by chemical environment (atom addition); a
heavy-atom-only ligand has no hydrogens for it to match or place. Note that
acpype will *accept* a bare `.pdb` and silently produce a heavy-atom-only
parameterisation — always go via the SDF from `rbfe-prep ligand`.

Then edit `systems/7abc/system.ini`. `scaffold` fills in everything except the
mutation, which it leaves blank so that `rbfe` refuses to guess. Only these need
thought:

```ini
[system]
name        = 7abc

[ligand]
dir     = inputs
ref     = ref
mut     = mut
resname = UNL          # the residue name your ligand .rtf declares

[protein]
pdb   = inputs/protein.pdb
chain = A              # '' takes every chain

[mutation]
strategy = element_swap        # or atom_addition, or mcs
vanish   = I1                  # atoms only the REFERENCE has
appear   = Br1                 # atoms only the MUTANT has
```

Appendix A lists every other key; all of them have working defaults.

### Choosing the mutation strategy

| strategy | when | `vanish` / `appear` |
|---|---|---|
| `element_swap` | one atom becomes a different element at the same position (I → Br, Cl → F). Nothing is added or removed. | atom names, identical in both ligands |
| `atom_addition` | the mutant has an atom the reference lacks (N–CH₃ → N–H). The new atom's position is **computed** by a rigid alignment. | `vanish` = the whole leaving group; `appear` = the mutant's atom that the new one replaces in kind |
| `mcs` | anything the other two do not cover: a group added *and* removed, a methyl → ethyl, a ring edit, or two ligands from **different parameterisation runs** that share no atom names. | **derived, and must not be declared** |

`mcs` finds the maximum common substructure (RDKit) and reads the mapping off
it, so it needs no declaration at all:

```ini
[mutation]
strategy      = mcs
missing_terms = synth_from_geometry
```

Placement is still a local rigid frame at each attachment point, so an appearing
atom inherits the mutant's own bond length and angle. Verified against both
worked systems: `mcs` reproduces `element_swap`'s 4YLJ hybrid and
`atom_addition`'s 6I5I hybrid **byte for byte**, including the placed coordinate.

> **`mcs` is a fallback, not a replacement.** The MCS fixes *which* atoms are
> common but not *which maps to which*: a symmetric core admits several equally
> large correspondences, and the wrong one still builds a valid-looking hybrid of
> the wrong mutation. On 6I5I the first correspondence RDKit returns swaps the
> pyrazole's two nitrogens — moving the methyl rather than removing it — so the
> strategy scores every correspondence by CHARMM type agreement and takes the
> best. That resolves 6I5I correctly (34 of 40 atom pairs agree), but it is a
> heuristic. **Read the `atom mapping` block it prints.** If the derived
> `vanish`/`appear` are not what you meant, do not tune it — write the
> declaration and use `element_swap` or `atom_addition`.

**An unrecognised strategy is a hard error.** The engine will not guess at
chemistry — that is the entire point of naming one.

---

## 3. The commands

Run them from inside the system directory. Each one discovers `system.ini` by
walking up, so no path arguments are needed.

### `rbfe similarity` — is this pair suitable for RBFE at all?

```bash
rbfe similarity
```

```
ligand similarity    : ECFP4 0.714 (threshold 0.600) -- ok
                       ECFP6 0.631   FCFP4 1.000   MACCS 0.938
```

Scores the two ligands and warns when they are too different for a *relative*
free energy to mean anything. `rbfe hybrid` runs the same screen before it
builds, so this command exists to answer the question **before** you invest in
inputs. It never fails and never changes an exit code — see `[similarity]` in
Appendix A for the threshold and what a low score implies.

### ① `rbfe hybrid` — build the dual-topology hybrid ligand

```bash
cd systems/7abc
rbfe hybrid
```

Writes `hybrid/hybrid.{pdb,rtf,prm}` and prints a report. **Read the report.**
Four lines matter:

```
atom mapping         : 29 common, 1 vanish ['I1'], 1 appear ['Br1']
charge closure       : ... lambda=0 -0.001000   lambda=1 -0.001000
parameters           : 34 bonds, 56 angles, 90 dihedrals, 18 impropers
                       1 term(s) with no source parameter (1 alchemical)
```

- **atom mapping** — did it find what you meant? A wrong `vanish`/`appear` shows
  up here as the wrong atom names. Under `mcs` this is the *only* check on the
  mapping, because nothing was declared — read it every time.
- **charge closure** — `lambda=0` must equal the reference's net charge and
  `lambda=1` the mutant's. If they do not, the build stops: a hybrid whose
  endpoints are not the physical molecules gives a meaningless ΔΔG.
- **terms with no source parameter** — every term that exists in *neither*
  endpoint state, and what was done about it. For an element swap this is the
  angle between the two halogens, a species no real ligand contains. Each line
  says `[mirrored]`, `[synthesised]`, or names the ambiguity it resolved.

`hybrid.pdb`'s B-factor column is the contract for everything downstream:

| B | meaning |
|---|---|
| `-1` | vanishing (reference only) |
| `0` | common |
| `+1` | appearing (mutant only) |

Open it in VMD and confirm the markers are where you expect **before** spending
GPU time. This is the single highest-value 60 seconds in the whole procedure.

### ② `rbfe build` — psfgen → solvate → ionize

```bash
rbfe build
```

Builds `complex/` (protein + hybrid) and `solvent/` (hybrid alone). Writes
`<leg>/{psf,pdb}` and `<leg>/box.json`.

The psfgen segment split is **derived** from the protein's numbering breaks, so
a structure with a disordered gap works without hand-editing anything. The cell
that VMD produces is written to `box.json` and read back in step ③ — the cell is
never a constant in a config.

Takes minutes. `--legs complex` builds only one leg.

### ③ `rbfe inputs` — the NAMD runtime files

```bash
rbfe inputs
```

Writes, once: `fep.tcl` (generated from `[run] lambdas`). Per leg:
`ionized.fep` (the `alchCol B` marking) and the four `.namd` configs.

Two checks run here and both are load-bearing:

- **the marked atoms are checked against the built system.** psfgen upper-cases
  atom names, so an element-swap ligand comes back as `BR1` where the hybrid
  wrote `Br1`. Matching those by exact string silently finds no appearing atom
  and writes an `alchFile` that mutates in one direction only — the run looks
  completely normal and the ΔΔG is meaningless.
- **the cell comes from `box.json`**, and the step refuses to run without it.

### ④ `rbfe check` — before burning GPU time

```bash
rbfe check                 # parse every config with --dryrun
rbfe check --smoke         # ... and time a short real run
```

`--dryrun` makes NAMD read and validate all 8 configs without integrating, so a
mistake costs seconds. `--smoke` then runs a short `nvt_equil` and reports
sec/step and ns/day, from NAMD's `TIMING:` line rather than the cumulative
`PERFORMANCE:` line (which is an average over the whole run and reads high
early). It projects the full job and warns if it exceeds `[run] max_hours`.

### ⑤ `rbfe run` — production

```bash
rbfe run                          # every leg and direction, in order
rbfe run complex forward          # one direction
rbfe run --respawn                # start the self-healing wrapper
rbfe run --watchdog 300           # start the read-only status watchdog
```

The leg order is fixed: **complex forward → complex backward → solvent forward →
solvent backward**. Both directions of the complex complete before the solvent
starts, because a crash before the cheap leg costs least.

Runs one NAMD process **per lambda window**, each with its own restart and a
`.done` marker, so re-running the same command after any interruption skips what
is finished. `rbfe status` shows where you are:

```
4ylj: 15 windows per direction
  complex  forward   [####...........] 4/15  <- next: w04
  complex  backward  [...............] 0/15  <- next: w00
  solvent  forward   [...............] 0/15  <- next: w00
  solvent  backward  [...............] 0/15  <- next: w00
```

Completion is judged by the **markers**, never by whether a `.fepout` exists — a
killed window leaves a partial file behind, and a file-exists test would declare
the leg complete and assemble a truncated one.

### ⑥ `rbfe audit` — the analysis

```bash
rbfe audit                 # full audit + bootstrap error bars
rbfe audit --skip-boot     # fast, no error bars
rbfe audit --scan          # ΔΔG sensitivity to the equilibration trim
```

`rbfe audit` is the estimator to trust. It uses a different parser and a
different BAR solver from the quick path, and it validates the whole dE
extraction against NAMD's own free-energy output. Sections:

| | |
|---|---|
| A / B | per-leg and net ΔΔG |
| C | **hysteresis** — forward vs backward. Should be inside the error. |
| D | bootstrap error bars on ΔΔG |
| E | **stationarity** — first vs second half of each window. Should be ~0. |
| F | overlap — spread of dE within a window, in kT |

Report the honest reading. A CI that spans zero is a null result, and §6 below
is what a healthy one looks like.

---

## 4. Reproducing the frozen systems

```bash
rbfe verify                          # both cases, tier 1, seconds
rbfe verify --case 4ylj --deep       # also re-run VMD (minutes; stop the job first)
```

This regenerates each frozen system's artefacts into a scratch directory and
compares them, and it **asserts that the frozen trees were not modified** by
hashing every tracked file before and after. Any difference that is not listed
in `verify/expected/deltas.json` fails.

Current state: **30 checks, 0 failures.** Byte-identical: `hybrid.pdb` for 6I5I,
`hybrid.{rtf,prm}` for 4YLJ, all `ionized.fep`, and all 12 `.namd` (6I5I's are
byte-identical including their relative parameter paths; 4YLJ's differ only in
the path to *this* build's `hybrid.prm`, which is the scratch location).

The deltas that remain are enumerated, with reasons, in
`verify/expected/deltas.json`. The two that matter:

- **4YLJ's `hybrid.pdb` differs on one line.** `element()` used to take the
  first letter of the GAFF type, so `br` became **boron**. It now derives the
  element from the MASS. The wrong value propagated into `complex.pdb` and
  `ionized.pdb`, so this is a correction.
- **6I5I's `hybrid.rtf` carries the corrected charges.** The frozen file predates
  6I5I's own charge-closure fix: its λ=1 state sums to **+0.0162 e**, a state
  that is supposed to be a neutral ligand. Its source comment describes exactly
  this bug; the file on disk is older than the fix.

---

## 5. Running on a machine that reboots

If the box kills long jobs (this one does), run the wrapper instead of `rbfe run`:

```bash
setsid nohup rbfe run --respawn > controller_run.log 2>&1 &
setsid nohup rbfe run --watchdog 300 &
tail -f controller_run.log
```

`--respawn` relaunches the controller whenever it exits before the job is done,
and gives up only if it exits twice in a row within 120 s (a real bug rather
than an eviction). `--watchdog` appends one line every 5 minutes to
`checkpoint_status.log` with the active window, its step and its restart
freshness — the only forensic record of where a kill landed.

Two caveats worth knowing:

- **Nothing keeps the watchdog alive.** A container recycle kills it along with
  everything else; the wrapper restarts the controller, not the watchdog.
- **A killed window replays from step 0.** It restarts from the previous
  window's endpoint, not from its own partial restart, so a kill 20 minutes into
  a 37-minute window still costs the whole window. That is the intended trade.

---

## 6. What a healthy result looks like

The 6I5I worked example, with its warts:

```
dG_complex +3.8816   dG_solvent +3.9878   ΔΔG −0.1062 kcal/mol
bootstrap SD 0.1199   95% CI [−0.3552, +0.1182]
```

**Reversible but not converged**, and the split is the whole story:

- ✅ hysteresis **+0.0587** (complex) / **−0.0920** (solvent) — both inside the
  ±0.120 error, randomly signed per window and largely cancelling.
- ❌ stationarity **+0.3653** (complex) / **−0.0065** (solvent) — **3× the error,
  3.4× the answer**, and *systematic*, so it accumulates instead of cancelling.
  Window 11 alone contributes +0.258. The complex leg is under-equilibrated.

The ΔΔG is a near-cancellation of two ~2.3 kcal/mol halves that flip sign at
λ≈0.5. It is a null result, and it should be reported as one.

---

## 7. Adding a mutation the engine does not know

**Try `strategy = mcs` first.** It derives the mapping from the two structures
and covers arbitrary substitutions — a group added *and* removed, a ring edit, a
methyl → ethyl, ligands parameterised independently. Only write a new strategy
when `mcs` refuses, or when it picks a mapping you cannot accept (a symmetric
core; §2 explains how to recognise that).

If `rbfe hybrid` stops with *unknown [mutation] strategy*, the answer is to add
a strategy, **not** to edit the builder. A strategy is small: it says which
atoms vanish, which appear, what they are called, and where the appearing ones
go. Parsing, connectivity enumeration, charge closure, parameter sourcing and
emission are shared.

1. Write `rbfe/mutations/<name>.py` with `name`, `required_keys`,
   `optional_keys`, `key_types`, `from_spec()` and `map()`.
2. Decorate the class with `@register`; import it in
   `rbfe/mutations/__init__.py`.
3. Set `Mapping.share_names = False` if the two ligands do **not** name their
   atoms consistently — otherwise the mutant's bond list is merged into the
   hybrid and invents connectivity that does not exist.
4. Add it to the table in §2.

Anything the strategy cannot derive, it must **refuse** with a message naming
the atom, not guess. `atom_addition` is worked example: it derives the anchor,
the local frame and the new atom's name, and demands `frame_atoms` /
`anchor` / `new_atom_name` from the config when the derivation is ambiguous.

---

## Appendix A — `system.ini` reference

Every key, its type and its default. Comments must be on their own line
(configparser does not strip inline `#`, which is deliberate — paths may
contain one). Unknown keys are fatal, with a "did you mean" suggestion.

### `[system]`
| key | type | default | |
|---|---|---|---|
| `name` | str | **required** | used in the audit banner |
| `description` | str | `''` | free text |

### `[ligand]`
| key | type | default | |
|---|---|---|---|
| `dir` | path | `inputs` | relative to the system root |
| `ref` / `mut` | str | **required** | stem: `<dir>/<stem>.{rtf,prm,pdb}` |
| `resname` | str | **required** | the residue name the hybrid is written under |

### `[protein]`
| key | type | default | |
|---|---|---|---|
| `pdb` | path | **required** | |
| `chain` | str | **required** | `''` takes every chain |
| `segment_prefix` | str | `P` | psfgen segment names `P1`, `P2`, … |
| `first` / `last` | str | `NTER` / `CTER` | psfgen terminus types |

### `[mutation]`
| key | type | default | |
|---|---|---|---|
| `strategy` | str | **required** | `element_swap` \| `atom_addition` \| `mcs` |
| `vanish` / `appear` | list | `''` | atom names; see §2. **Not accepted by `mcs`**, which derives them |
| `placement` | str | `native` | `native` (element_swap) \| `local_frame` (atom_addition, mcs) |
| `missing_terms` | str | `fail` | `fail` \| `mirror` \| `mirror_then_synth` \| `synth_from_geometry` |
| `anchor` | str | derived | atom_addition: the atom the leaving group hangs off |
| `frame_atoms` | list | derived | atom_addition: the 2 atoms defining the alignment frame |
| `new_atom_name` | str | derived | atom_addition, mcs: the hybrid name (next free `H<n>`) |
| `improper` | list | derived | atom_addition: 4 atoms for the planarity improper |
| `mcs_timeout` | int | `60` | mcs: seconds before the substructure search gives up |

Only the keys the named strategy declares are accepted; a key from a different
strategy is rejected like any other unknown key. A `mcs_timeout` expiry is a
hard failure, not a warning: a truncated search returns a *smaller* common core,
which would silently enlarge the alchemical transformation.

### `[similarity]`

Screens the ligand pair before it is built, because RBFE is a *relative* method
and is only meaningful between ligands similar enough to share configuration
space. It **warns and continues** — it never fails a build and never changes an
exit code.

| key | type | default | |
|---|---|---|---|
| `threshold` | float | `0.60` | ECFP4 Tanimoto below which the pair is queried |
| `enabled` | bool | `yes` | `no` skips the screen entirely |

The two worked systems sit just above the threshold — 4YLJ (I→Br) 0.714, 6I5I
(N–CH₃→N–H) 0.719 — while a cross-family pair scores 0.105–0.117. `rbfe
similarity` reports the same number without building anything.

RDKit is an **optional** dependency: it is imported lazily, so the rest of the
pipeline runs without it. When it is absent, or a ligand will not build, the
screen says `similarity NOT checked` on stderr rather than passing quietly.

> **A low score is a stop-and-think, not a refusal.** Past this point the two end
> states are effectively different molecules and the resulting ΔΔG is not a
> relative binding free energy. **ABFE** is the appropriate method for a change
> that large, and it is out of scope for this repo.

### `[hybrid]`
| key | type | default | |
|---|---|---|---|
| `title` | str | `dual-topology hybrid ligand` | line 1 of the `.rtf`/`.prm` |
| `type_prefix` | str | `L` | `raw` reproduces a legacy build exactly. **See the warning below.** |
| `rtf_order` | str | `sorted` | `sorted` \| `source_then_new` |
| `prm_order` | str | `source_then_new` | `.rtf` and `.prm` are ordered independently |

> **`type_prefix = raw` is a compatibility mode, not a default.** NAMD matches
> atom types case-insensitively and takes the *last* duplicate, and
> `hybrid.prm` is listed last — so a ligand's GAFF types (`ca`, `cc`, `cd`,
> `ha`, `os`) silently override CHARMM36's **protein** types of the same name.
> 6I5I's published build does exactly this: 7 `DUPLICATE ... ENTRY` warnings per
> process, and the aromatic sidechains of PHE/TYR/TRP/HIS simulated with the
> ligand's parameters. The `L` prefix keeps the ligand in its own namespace and
> is what new systems should use.

### `[build]`
| key | type | default | |
|---|---|---|---|
| `padding` | float | `15.0` | `solvate -t`, Angstrom |
| `neutralize` | bool | `yes` | |
| `salt` | float | `0.0` | molar; `0` omits the flag |
| `pdbalias` | lines | ILE CD1/SER HG/HIS HSD | one `pdbalias` argument per line |

### `[topology]`
| key | type | default | |
|---|---|---|---|
| `prot_rtf` | path | **required** | `{root}` expands to the system dir |
| `params` | lines | **required** | emitted into every `.namd`; `{hybrid}` expands to the built `hybrid.prm` |

A `params` entry containing `{root}` is resolved to an absolute path, so the
config works at any depth. An entry **without** `{root}` is emitted verbatim —
that is the escape hatch for reproducing a legacy config whose paths were
relative.

### `[run]`
| key | type | default | |
|---|---|---|---|
| `legs` | list | `complex solvent` | |
| `lambdas` | floats | **required** | must run 0.0 → 1.0, strictly ascending |
| `steps_per_window` | int | **required** | 500000 = 1 ns at 2 fs |
| `equil_steps` | int | `50000` | nvt and npt |
| `alch_equil_steps` | int | `50000` | per window |
| `outfreq` | int | `500` | `steps_per_window` must divide by it |
| `min_steps` | int | `5000` | minimisation before nvt |
| `timestep` | float | `2.0` | |
| `temperature` | float | `300.0` | |
| `cutoff` / `switchdist` / `pairlistdist` | float | `12.0` / `10.0` / `14.0` | |
| `elec_lambda_start` / `vdw_lambda_end` | float | `0.5` / `1.0` | |
| `vdw_shift_coeff` | float | `5.0` | |
| `decouple` | str | `off` | |
| `max_hours` | float | `0.0` | `0` disables the wall-clock warning |

### `[binaries]`
| key | type | default | |
|---|---|---|---|
| `namd` | str | **required** | `${NAMD:-...}` supported |
| `vmd` | str | `vmd` | |
| `flags` | list | `+p1 +devices 0` | separate tokens |
| `gpu_resident` | str | `on` | |

---

## Appendix B — traps that have already bitten

Each of these produced a plausible-looking wrong answer at least once.

1. **`+p1`, never `+p8`.** In GPU-resident mode throughput is *inversely*
   proportional to PE count: 75.8 ns/day at `+p1` vs 13.2 at `+p8`. `+devices`
   and `0` must be separate argv tokens.
2. **Never set `temperature` on a stage that reads `binvelocities`.** NAMD dies
   with *"Cannot specify both an initial temperature and a velocity file."* Only
   `nvt_equil` sets it.
3. **The first 99 of every window's 500 rows are equilibration.** `alchEquilSteps
   50000` at `alchOutFreq 500`. NAMD's own `dE_avg`/`dG` columns exclude them;
   any analysis must too. This is derived, not copy-pasted, and the trim is
   printed by `rbfe audit`.
4. **A backward `.fepout` pairs to window *i* by reversed index** and its `dE`
   must be **negated** — a backward run samples B and reports `U(A)−U(B)`.
5. **psfgen upper-cases atom names.** `Br1` becomes `BR1`. Match case-insensitively
   and verify the marked atoms against the built system.
6. **`GPUresident` in exactly one place** — the config *or* the command line,
   never both.
7. **`par_water_ions_clean.prm`, not the `.str`.**
8. **Scraping the "Free energy change" text yields the cumulative
   `net change until now`** — which for 6I5I reads −16.821 against a true −5.455,
   and looks like a normal result. Parse the per-window columns, not the prose.

---

## Appendix C — when a step refuses to run

Every refusal names the thing that is wrong. The common ones:

| message | what it means |
|---|---|
| `unknown key 'x' in [section]` | typo; the message suggests the nearest legal key |
| `unknown [mutation] strategy` | the strategy is not implemented — see §7 |
| `no {config} found in ... or any parent` | run from inside a system directory |
| `[mutation] vanish names atom(s) not in ref` | the atom name is not in the ligand `.rtf` |
| `charge closure failed` | the ligand charges do not sum to the right net charge |
| `no path parameter for ...` | a term exists in neither ligand; set `missing_terms` |
| `cannot derive the local frame` | ambiguous — set `frame_atoms` explicitly |
| `carries no alchemical atom named` | psfgen renamed it, or the hybrid is wrong |
| `box.json missing` | run `rbfe build` first |
| `fep.tcl does not match [run] lambdas` | re-run `rbfe inputs` |
