"""Is this ligand pair suitable for RBFE at all?

Relative binding free energy is a *relative* method: it computes the difference
between two end states by alchemically mutating one into the other, and that is
only meaningful when the two ligands are similar enough that their configurations
overlap.  A pair that is too different still builds a perfectly valid hybrid and
still produces a converged-looking number -- the number is simply not the thing
the method assumes it is.  Nothing else in the pipeline notices, so this module is
where the question gets asked.

Where the molecule comes from
-----------------------------
From the engine's **own parsed `Ligand`** (`rbfe.charmm.load_ligand`), not from a
separate structure file: elements from the `.rtf` MASS line, bonds from the `.rtf`,
coordinates from the `.pdb`.  Two reasons, both empirical:

* The `.mol2` sitting beside the inputs cannot be read by RDKit at all.  GAFF
  SYBYL types are not element symbols, and `ca` (aromatic carbon) fails the
  element lookup -- reading as *calcium*.  `MolFromMol2File` returns None; a
  hand-patched parse gave 6I5I a Tanimoto of 0.333, a false rejection of a
  textbook RBFE pair.  `.sdf` is no better as a general source: only 4YLJ has one.
* Building from the topology means the fingerprint describes the same molecule the
  simulation does.  It also sidesteps a bad element column: `4YLJ/inputs/mut.pdb`
  writes `I` in columns 77-78 for the atom named `Br1`, and because the mass is the
  authority here, the bromine still comes out bromine.

Which fingerprint
-----------------
ECFP4 (Morgan, radius 2, 2048 bits) is the verdict.  **Fast filtered fingerprints
are deliberately reported but never scored**: FCFP4 rates the 4YLJ I->Br mutation
at 1.000, because iodine and bromine are the same pharmacophore feature there.  A
metric that cannot see the only difference between the two ligands is the wrong
one to gate on.

Calibration: both worked systems sit just above the default threshold -- 4YLJ
(I->Br) 0.714, 6I5I (N-CH3 -> N-H) 0.719 -- while a cross-family pair scores
0.105-0.117.  So 0.60 accepts both known-good pairs with margin without being
vacuous.  It is a heuristic, not a law: it screens, it does not decide.

Failing closed
--------------
RDKit is an optional dependency.  Importing it at module scope would break
`rbfe.config` -- and therefore every subcommand -- on a machine that does not have
it, so every import here is lazy.  When RDKit is missing, or a ligand will not
build, the check reports that it could NOT run, on **stderr**, and returns None.
It never reports a pass it did not perform.
"""

from __future__ import annotations

import sys

from rbfe.model import Ligand

# ECFP4 above this is taken as "similar enough for RBFE".  See the calibration
# note in the module docstring before moving it.
DEFAULT_THRESHOLD = 0.60

# Morgan radius and bit count for the scored fingerprint.
RADIUS = 2
NBITS = 2048


class Unavailable(Exception):
    """RDKit is absent, or a ligand could not be turned into a molecule."""


def available() -> bool:
    """Whether RDKit can be imported.  Cheap enough to call per run."""
    try:
        import rdkit  # noqa: F401
    except ImportError:
        return False
    return True


def build_mol(ligand: Ligand, keep_h: bool = False):
    """An RDKit molecule for `ligand`, built from its parsed topology.

    Every atom carries a `charmm_name` property, so a result can be mapped back
    onto the hybrid's names.  Raises `Unavailable` rather than returning None --
    a caller that cannot fingerprint must say so, not score zero.

    Bond orders are not in the `.rtf`, so they are inferred from the 3D geometry
    by `rdDetermineBonds`; the ligand's own charge total pins the formal charges.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    from rbfe.charmm import element

    rw = Chem.RWMol()
    index: dict[str, int] = {}
    for name in ligand.names:
        atom = Chem.Atom(element(ligand.atoms[name]))
        # Hydrogens are explicit in the .rtf; stop RDKit adding a second set.
        atom.SetNoImplicit(True)
        idx = rw.AddAtom(atom)
        rw.GetAtomWithIdx(idx).SetProp("charmm_name", name)
        index[name] = idx

    for u, v in ligand.bonds:
        if u in index and v in index:
            rw.AddBond(index[u], index[v], Chem.BondType.SINGLE)

    mol = rw.GetMol()
    if ligand.xyz:
        conf = Chem.Conformer(mol.GetNumAtoms())
        for name, idx in index.items():
            x, y, z = ligand.xyz[name]
            conf.SetAtomPosition(idx, (float(x), float(y), float(z)))
        mol.AddConformer(conf)

    mol.UpdatePropertyCache(strict=False)
    charge = round(ligand.total_charge)
    try:
        rdDetermineBonds.DetermineBondOrders(mol, charge=charge)
        Chem.SanitizeMol(mol)
    except Exception as exc:                      # noqa: BLE001 - rdkit raises widely
        raise Unavailable(
            f"could not infer bond orders for {ligand.stem} "
            f"(charge {charge:+d}): {exc}") from None

    if not keep_h:
        mol = Chem.RemoveHs(mol)
    return mol


def tanimoto(ref: Ligand, mut: Ligand) -> dict[str, float]:
    """Tanimoto scores for the pair.  Raises `Unavailable` if it cannot score.

    `ecfp4` is the one that counts; the rest are reported so a surprising verdict
    can be diagnosed rather than argued with.
    """
    from rdkit import DataStructs
    from rdkit.Chem import MACCSkeys

    try:
        from rdkit.Chem import rdFingerprintGenerator as G
    except ImportError:
        # The legacy AllChem.GetMorganFingerprintAsBitVect would still work, but
        # it warns on every call and the warnings would bury this report.
        raise Unavailable(
            "this RDKit predates rdFingerprintGenerator, which this check needs"
        ) from None

    a, b = build_mol(ref), build_mol(mut)

    def morgan(radius: int, features: bool = False):
        kwargs = {"radius": radius, "fpSize": NBITS}
        if features:
            kwargs["atomInvariantsGenerator"] = G.GetMorganFeatureAtomInvGen()
        gen = G.GetMorganGenerator(**kwargs)
        return gen.GetFingerprint(a), gen.GetFingerprint(b)

    score = DataStructs.TanimotoSimilarity
    e4a, e4b = morgan(2)
    e6a, e6b = morgan(3)
    f4a, f4b = morgan(2, features=True)
    return {
        "ecfp4": score(e4a, e4b),
        "ecfp6": score(e6a, e6b),
        "fcfp4": score(f4a, f4b),
        "maccs": score(MACCSkeys.GenMACCSKeys(a), MACCSkeys.GenMACCSKeys(b)),
    }


def _unavailable(reason: str) -> None:
    """Say plainly, on stderr, that the check did not happen."""
    print(f"  similarity NOT checked: {reason}", file=sys.stderr)
    print("           this run has NOT been screened for RBFE suitability.",
          file=sys.stderr)


def check(ref: Ligand, mut: Ligand, *, threshold: float = DEFAULT_THRESHOLD,
          printer=print) -> dict | None:
    """Score the pair and report it.  Returns the result, or None if unchecked.

    Warns and returns; it never raises on a low score and never changes an exit
    code.  A pair below the threshold is reported loudly and the build continues,
    because the threshold is a heuristic and the operator may know better -- the
    point is that they are told.
    """
    try:
        scores = tanimoto(ref, mut)
    except Unavailable as exc:
        _unavailable(str(exc))
        return None
    except ImportError:
        _unavailable("RDKit is not installed")
        return None

    value = scores["ecfp4"]
    passes = value >= threshold
    verdict = "ok" if passes else "BELOW THRESHOLD"
    printer(f"ligand similarity    : ECFP4 {value:.3f} "
            f"(threshold {threshold:.3f}) -- {verdict}")
    printer(f"                       ECFP6 {scores['ecfp6']:.3f}   "
            f"FCFP4 {scores['fcfp4']:.3f}   MACCS {scores['maccs']:.3f}")

    if not passes:
        printer("")
        printer(f"  WARNING: {ref.stem} and {mut.stem} are too different for a "
                f"relative free energy calculation.")
        printer(f"           ECFP4 Tanimoto {value:.3f} is below {threshold:.3f}. "
                f"RBFE mutates one ligand into the")
        printer("           other and relies on the two sharing enough "
                "configuration space for that")
        printer("           to be meaningful; past this point the two end states "
                "are effectively")
        printer("           different molecules and the resulting ddG is not a "
                "relative binding free")
        printer("           energy. ABFE (absolute binding free energy) is the "
                "appropriate method")
        printer("           for a change this large -- it is out of scope for "
                "this repo.")
        printer("")
        printer("           Continuing anyway: this is a warning, not a refusal. "
                "Raise")
        printer("           [similarity] threshold to silence it if the pair is "
                "deliberate.")
        printer("")

    return {**scores, "threshold": threshold, "passes": passes}
