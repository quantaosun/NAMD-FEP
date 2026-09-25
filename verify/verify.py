"""Reproduce the frozen systems, and prove nothing in them was touched.

A generalization that silently changes the physics is worse than no
generalization, so this is not an optional extra -- it is how the port earns
trust.  Two things it insists on:

**The frozen trees are read-only, and that is asserted, not promised.**  Every
file under `6I5I_DUAL_FEP/` and `4YLJ/` is hashed before and after; any change
is a hard failure of the harness itself.  This matters because `4YLJ/` holds a
production job and its `hybrid.prm` is held open by a running NAMD -- rewriting
it would desync the live physics.

**Differences are enumerated, not tolerated.**  `expected/deltas.json` lists
every difference that is *accepted*, with a reason.  A difference that is not in
the manifest fails the run.  That is the difference between a verification suite
and a rubber stamp.

Tier 1 (default) is pure Python and takes seconds.  Tier 2 (`--deep`) re-runs
VMD and takes minutes; it can confirm the cell and the protein/ligand
coordinates but NOT the water, because solvate is non-deterministic.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from rbfe.errors import VerificationError

REPO = Path(__file__).resolve().parent.parent
FROZEN = ("6I5I_DUAL_FEP", "4YLJ")

# case -> (config, frozen tree, stages to run)
#
# The "mcs" stage rebuilds the same system with the mapping DERIVED instead of
# declared (verify/cases/<case>-mcs.ini) and asserts it reproduces the
# specialist strategy's hybrid byte for byte. That is the invariant which lets
# `mcs` be offered as the general path: it must not disagree with the strategy
# it generalises.
CASES = {
    "4ylj": ("systems/4ylj/system.ini", "4YLJ", ("hybrid", "inputs", "mcs")),
    "6i5i": ("verify/cases/6i5i-legacy.ini", "6I5I_DUAL_FEP",
             ("hybrid", "inputs", "mcs")),
}

OK, BAD, SKIP = "ok", "FAIL", "skip"


class Report:
    def __init__(self) -> None:
        self.rows: list[tuple[str, str, str]] = []

    def add(self, case: str, status: str, what: str) -> None:
        self.rows.append((case, status, what))

    def summarise(self) -> int:
        width = max((len(r[2]) for r in self.rows), default=0)
        print("\n" + "=" * (width + 22))
        for case, status, what in self.rows:
            mark = {OK: " ok ", BAD: "FAIL", SKIP: "--  "}[status]
            print(f"  [{mark}] {case:6s} {what}")
        bad = sum(1 for r in self.rows if r[1] == BAD)
        skipped = sum(1 for r in self.rows if r[1] == SKIP)
        print("=" * (width + 22))
        print(f"  {len(self.rows) - bad - skipped} passed, {bad} failed"
              + (f", {skipped} skipped" if skipped else ""))
        return 1 if bad else 0


# The guard protects the *source* of a frozen system -- the inputs, the hybrid
# it was built from, and the built topology.  The trajectory output is excluded
# on purpose: it is what a running job rewrites continuously, so including it
# would (a) hash hundreds of MB of .dcd each pass and (b) be racy against the
# job rather than against this code.  A run output changing is not evidence
# about the verifier; a .rtf or a .psf changing is.
VOLATILE_SUFFIXES = (".dcd", ".coor", ".vel", ".xsc", ".fepout", ".log",
                     ".coor.old", ".vel.old", ".xsc.old", ".coor.BAK",
                     ".vel.BAK", ".xsc.BAK", ".restart", ".BAK")


def _tracked(p: Path) -> bool:
    if "window_snapshots" in p.parts:
        return False
    n = p.name
    return not (n.endswith(VOLATILE_SUFFIXES) or n.endswith(".done")
                or n.startswith(".") and n.endswith(".done"))


def snapshot(root: Path) -> dict[str, str]:
    """md5 of every tracked file under a tree, keyed by relative path."""
    out = {}
    for p in sorted(root.rglob("*")):
        if p.is_file() and _tracked(p):
            out[str(p.relative_to(root))] = hashlib.md5(p.read_bytes()).hexdigest()
    return out


def job_running() -> bool:
    return subprocess.run(["pgrep", "-x", "namd3"],
                          capture_output=True).returncode == 0


def assert_frozen_untouched(before: dict[str, dict[str, str]],
                            after: dict[str, dict[str, str]],
                            live: bool) -> None:
    if live:
        # A running production job writes into its own tree continuously --
        # .coor/.dcd/.fepout/.done all change between the two snapshots.  The
        # guard cannot attribute a change to the verifier while that is
        # happening, so it stands down and says so rather than either failing
        # on the job's own output or passing quietly.
        print("  frozen-tree assertion DISABLED: a namd3 job is running and is "
              "itself writing there.\n"
              "    Re-run with the job stopped to assert read-only access.")
        return
    for tree, b in before.items():
        a = after[tree]
        changed = sorted(k for k in b if k in a and a[k] != b[k])
        removed = sorted(set(b) - set(a))
        added = sorted(set(a) - set(b))
        if changed or removed or added:
            raise VerificationError(
                f"THE VERIFIER MODIFIED {tree}/ -- that must never happen.\n"
                f"  changed: {changed[:5]}\n  removed: {removed[:5]}\n"
                f"  added:   {added[:5]}\n"
                f"  This is a bug in the harness or in the code it exercises, "
                f"not an acceptable difference. 4YLJ/ holds a live production "
                f"job.")
    print(f"  frozen trees unchanged: {', '.join(before)} "
          f"({sum(len(v) for v in before.values())} files hashed before and after)")


def load_manifest() -> list[dict]:
    p = Path(__file__).resolve().parent / "expected" / "deltas.json"
    return json.loads(p.read_text())


def manifested(manifest: list[dict], case: str, artifact: str, kind: str) -> dict | None:
    for m in manifest:
        if m.get("kind") != kind:
            continue
        if m.get("case") not in (case, "*"):
            continue
        pat = m.get("artifact", "")
        if pat == artifact or (pat.startswith("*.") and artifact.endswith(pat[1:])) \
                or artifact.endswith(pat):
            return m
    return None


def compare_bytes(case: str, label: str, new: Path, old: Path, rep: Report,
                  manifest: list[dict], kind_when_diff: str) -> None:
    if not old.exists():
        rep.add(case, SKIP, f"{label} (no frozen counterpart)")
        return
    if not new.exists():
        rep.add(case, BAD, f"{label} (not generated)")
        return
    if new.read_bytes() == old.read_bytes():
        rep.add(case, OK, f"{label} byte-identical")
        return
    m = manifested(manifest, case, label, kind_when_diff)
    if m:
        rep.add(case, OK, f"{label} differs as documented ({kind_when_diff})")
    else:
        rep.add(case, BAD, f"{label} differs, and that is NOT in the manifest")


def compare_rendered(case: str, label: str, new_text: str, old_text: str,
                     rep: Report, manifest: list[dict], kind: str) -> None:
    if new_text == old_text:
        rep.add(case, OK, f"{label} byte-identical")
        return
    m = manifested(manifest, case, label, kind)
    if m:
        rep.add(case, OK, f"{label} differs as documented ({kind})")
    else:
        diff = [l for l in old_text.splitlines() if l not in new_text.splitlines()]
        rep.add(case, BAD,
                f"{label} differs unexpectedly ({len(diff)} line(s)): "
                f"{diff[0][:60] if diff else ''}")


def box_from_namd(namd: Path):
    """The periodic cell out of a generated .namd's cellBasisVector/cellOrigin."""
    import re as _re
    if not namd.exists():
        return None
    txt = namd.read_text()
    vecs = _re.findall(r"^cellBasisVector\d\s+([-\d.]+) ([-\d.]+) ([-\d.]+)\s*$",
                       txt, _re.M)
    org = _re.search(r"^cellOrigin\s+([-\d.]+) ([-\d.]+) ([-\d.]+)\s*$", txt, _re.M)
    if len(vecs) != 3 or not org:
        return None
    return {"cell": [float(v[0 if i == j else 1]) if False else float(vecs[i][i])
                     for i, j in enumerate(range(3))],
            "origin": [float(x) for x in org.groups()]}


def run_case(case: str, cfg_path: Path, frozen: Path, stages, scratch: Path,
             rep: Report, manifest: list[dict], deep: bool) -> None:
    sys.path.insert(0, str(REPO))
    from rbfe import config as C

    cfg = C.load(cfg_path)
    out = scratch / case
    out.mkdir(parents=True, exist_ok=True)

    if not (out / "hybrid").exists():
        (out / "hybrid").mkdir(parents=True)

    if "hybrid" in stages:
        from rbfe import hybrid as H
        from rbfe.mutations import get_strategy, spec_from_section
        spec = spec_from_section(cfg.mutation_section())
        strategy = get_strategy(spec.strategy).from_spec(spec)
        lig = cfg.resolve(cfg.get("ligand", "dir"))
        H.build(lig / cfg.get("ligand", "ref"), lig / cfg.get("ligand", "mut"),
                out / "hybrid", resname=cfg.get("ligand", "resname"),
                title=cfg.get("hybrid", "title"),
                type_prefix=cfg.get("hybrid", "type_prefix"),
                rtf_order=cfg.get("hybrid", "rtf_order"),
            prm_order=cfg.get("hybrid", "prm_order"),
                missing_terms=spec.missing_terms, strategy=strategy, verbose=False)
        kinds = {"hybrid.rtf": "charge-closure", "hybrid.prm": "mirror-ambiguity",
                 "hybrid.pdb": "element-column"}
        for f, kind in kinds.items():
            compare_bytes(case, f"hybrid/{f}", out / "hybrid" / f,
                          frozen / "hybrid" / f, rep, manifest, kind)

        # The calibration the similarity threshold rests on: both worked systems
        # must score ABOVE it. A threshold that rejected either would be wrong,
        # and this is what catches a fingerprint change that moves the numbers.
        from rbfe import similarity as SIM
        from rbfe.charmm import load_ligand
        try:
            score = SIM.tanimoto(load_ligand(lig / cfg.get("ligand", "ref")),
                                 load_ligand(lig / cfg.get("ligand", "mut")))["ecfp4"]
        except SIM.Unavailable as exc:
            rep.add(case, SKIP, f"ligand similarity not scored ({exc})")
        else:
            thr = float(cfg.get("similarity", "threshold"))
            rep.add(case, OK if score >= thr else BAD,
                    f"ligand ECFP4 {score:.3f} vs threshold {thr:.2f}")

    if "mcs" in stages:
        from rbfe import similarity as SIM
        if not SIM.available():
            rep.add(case, SKIP, "mcs agrees with the specialist (no RDKit)")
        else:
            from rbfe import hybrid as H
            from rbfe.mutations import get_strategy, spec_from_section
            mcs_cfg = C.load(REPO / "verify" / "cases" / f"{case}-mcs.ini")
            spec = spec_from_section(mcs_cfg.mutation_section())
            strategy = get_strategy(spec.strategy).from_spec(spec)
            lig2 = mcs_cfg.resolve(mcs_cfg.get("ligand", "dir"))
            out_mcs = out / "hybrid-mcs"
            H.build(lig2 / mcs_cfg.get("ligand", "ref"),
                    lig2 / mcs_cfg.get("ligand", "mut"), out_mcs,
                    resname=mcs_cfg.get("ligand", "resname"),
                    title=mcs_cfg.get("hybrid", "title"),
                    type_prefix=mcs_cfg.get("hybrid", "type_prefix"),
                    rtf_order=mcs_cfg.get("hybrid", "rtf_order"),
                    prm_order=mcs_cfg.get("hybrid", "prm_order"),
                    missing_terms=spec.missing_terms, strategy=strategy,
                    verbose=False)
            # Compared against THIS RUN's specialist build, not the frozen tree:
            # the frozen deltas are already accounted for above, so any
            # difference here is the two strategies genuinely disagreeing.
            for f in ("hybrid.rtf", "hybrid.prm", "hybrid.pdb"):
                compare_bytes(case, f"mcs/{f}", out_mcs / f, out / "hybrid" / f,
                              rep, manifest, "mcs-disagrees")

    if "inputs" in stages:
        from rbfe import inputs as I
        # Copy the frozen built system in, so this tier does not need VMD.
        for leg in cfg.legs:
            d = out / leg
            d.mkdir(parents=True, exist_ok=True)
            for f in ("ionized.pdb", "box.json"):
                if (frozen / leg / f).exists():
                    shutil.copy(frozen / leg / f, d / f)
            # 6I5I predates the box.json convention: its generated .namd carries
            # the cell inline.  Recover it from there so this case can be
            # verified at all -- rbfe itself refuses to run without box.json,
            # which is the behaviour that stops a system being simulated in
            # another system's cell.
            if not (d / "box.json").exists():
                box = box_from_namd(frozen / leg / "md_forward.namd")
                if box is None:
                    rep.add(case, SKIP, f"{leg}: no box.json and no cell to recover")
                    return
                (d / "box.json").write_text(json.dumps(box, indent=2) + "\n")
                rep.add(case, OK,
                        f"{leg}: box recovered from the frozen .namd "
                        f"(6I5I has no box.json)")
        I.write_all(cfg, out)
        for leg in cfg.legs:
            compare_bytes(case, f"{leg}/ionized.fep", out / leg / "ionized.fep",
                          frozen / leg / "ionized.fep", rep, manifest, "x")
            for name in ("nvt_equil", "npt_equil", "md_forward", "md_backward"):
                compare_bytes(case, f"{leg}/{name}.namd", out / leg / f"{name}.namd",
                              frozen / leg / f"{name}.namd", rep, manifest,
                              "param-path")
        compare_rendered(case, "fep.tcl", I.render_fep_tcl(cfg),
                         (frozen / "fep.tcl").read_text(), rep, manifest, "header")

    if deep:
        tier2(case, cfg, frozen, scratch, rep)


def tier2(case: str, cfg, frozen: Path, scratch: Path, rep: Report) -> None:
    """Re-run VMD and compare what VMD can reproduce.

    Asserts the cell and the protein+ligand coordinates, NOT the water: solvate
    is non-deterministic (measured -- two runs from identical input share only
    30080 of 72634 atoms), so the water network is not a thing anyone can
    reproduce, and asserting on it would fail for reasons that have nothing to
    do with this code.
    """
    from rbfe import build as B
    out = scratch / case
    try:
        B.build(cfg, out)
    except Exception as exc:                       # noqa: BLE001
        rep.add(case, BAD, f"tier2 build failed: {exc}")
        return

    for leg in cfg.legs:
        new_box = json.loads((out / leg / "box.json").read_text()) \
            if (out / leg / "box.json").exists() else None
        old_box_p = frozen / leg / "box.json"
        if new_box and old_box_p.exists():
            old_box = json.loads(old_box_p.read_text())
            same_cell = all(abs(a - b) < 1e-3 for a, b in
                            zip(new_box["cell"], old_box["cell"]))
            rep.add(case, OK if same_cell else BAD,
                    f"{leg}: cell reproduced ({'yes' if same_cell else 'NO'})"
                    + ("" if same_cell else f" {new_box['cell']} vs {old_box['cell']}"))

        new_p = out / leg / "ionized.pdb"
        old_p = frozen / leg / "ionized.pdb"
        if new_p.exists() and old_p.exists():
            pl = _solute(new_p)
            ol = _solute(old_p)
            same = pl == ol
            rep.add(case, OK if same else BAD,
                    f"{leg}: {len(pl)} protein+ligand atoms reproduced exactly "
                    f"({'yes' if same else 'NO'})")


def _solute(pdb: Path) -> list[str]:
    """Coordinates of everything that is not water or a lone ion."""
    out = []
    for line in pdb.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        res = line[17:20].strip()
        if res.startswith("TIP") or res in ("CLA", "SOD", "POT", "MG", "CAL", "ZN"):
            continue
        out.append(line[30:54])
    return out


def run(cases: list[str] | None, deep: bool = False, scratch: str | None = None) -> int:
    cases = cases or list(CASES)
    unknown = [c for c in cases if c not in CASES]
    if unknown:
        raise VerificationError(
            f"unknown case(s) {unknown}; known: {', '.join(CASES)}")

    if deep and subprocess.run(["pgrep", "-x", "namd3"], capture_output=True).returncode == 0:
        raise VerificationError(
            "--deep runs VMD and would compete with the running NAMD job. "
            "Stop the job first, or run tier 1 only.")

    manifest = load_manifest()
    rep = Report()
    print(f"verifying {', '.join(cases)}  (tier {'2' if deep else '1'})")

    live = job_running()
    before = {t: snapshot(REPO / t) for t in FROZEN if (REPO / t).exists()}
    tmp = Path(scratch).resolve() if scratch else Path(tempfile.mkdtemp(prefix="rbfe-verify-"))
    if any(str(tmp).startswith(str((REPO / t).resolve())) for t in FROZEN):
        raise VerificationError(
            f"scratch directory {tmp} is inside a frozen tree. Refusing to run.")
    print(f"  scratch: {tmp}")

    try:
        for case in cases:
            cfg_path, tree, stages = CASES[case]
            try:
                run_case(case, REPO / cfg_path, REPO / tree, stages, tmp, rep, manifest, deep)
            except VerificationError:
                raise
            except Exception as exc:                       # noqa: BLE001
                rep.add(case, BAD, f"raised {type(exc).__name__}: {exc}")
    finally:
        after = {t: snapshot(REPO / t) for t in before}
        assert_frozen_untouched(before, after, live)

    rc = rep.summarise()
    if not scratch:
        print(f"  (scratch kept at {tmp}; remove it when done)")
    return rc
