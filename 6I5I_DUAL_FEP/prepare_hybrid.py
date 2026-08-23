#!/usr/bin/env python3
"""
Dual-topology hybrid ligand builder for the 6I5I N-CH3 -> N-H FEP.

Reference ligand (H3E "12H") carries an N-methyl; the desmethyl mutant carries
an N-H in its place.  In NAMD's default *dual topology* scheme
(singleTopology off) the two ligands share a common core (present once) while
the atoms that differ are BOTH present and are switched on/off through the
alchFile B-factor column:

    B = -1  ->  reference-only atoms (the methyl group)  : vanish at lam=1
    B =  0  ->  common atoms (identical in both states)  : always on
    B = +1  ->  mutant-only atoms (the N-H hydrogen)     : appear at lam=1

Only VMD/NAMD-adjacent plain files (CHARMM rtf/prm/pdb from acpype) + numpy
are used.  No FEPrepare, no RDKit.

Outputs (./hybrid/):  hybrid.pdb (B-factors), hybrid.rtf, hybrid.prm.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
INP = HERE / "inputs"
OUT = HERE / "hybrid"


# ---------------------------------------------------------------------------
# CHARMM rtf / prm / pdb parsers
# ---------------------------------------------------------------------------

def parse_rtf(path: Path) -> dict:
    txt = Path(path).read_text().splitlines()
    mass: dict[str, float] = {}
    atoms: list[tuple[str, str, float]] = []          # (name, type, charge)
    bonds: list[tuple[str, str]] = []
    angles: list[tuple[str, str, str]] = []
    dihe: list[tuple[str, str, str, str]] = []
    impr: list[tuple[str, str, str, str]] = []
    for line in txt:
        s = line.strip()
        if not s or s.startswith("*"):
            continue
        p = s.split()
        kw = p[0].upper()
        if kw == "MASS" and len(p) >= 4:
            mass[p[2].lower()] = float(p[3])
        elif kw == "ATOM" and len(p) >= 4:
            atoms.append((p[1], p[2].lower(), float(p[3])))
        elif kw == "BOND" and len(p) >= 3:
            bonds.append((p[1], p[2]))
        elif kw == "ANGL" and len(p) >= 4:
            angles.append((p[1], p[2], p[3]))
        elif kw == "DIHE" and len(p) >= 5:
            dihe.append((p[1], p[2], p[3], p[4]))
        elif kw in ("IMPR", "IMPH") and len(p) >= 5:
            impr.append((p[1], p[2], p[3], p[4]))
    return {"mass": mass, "atoms": atoms, "bonds": bonds,
            "angles": angles, "dihedrals": dihe, "impropers": impr}


def _nt(t: str) -> str:
    """Normalise a CHARMM atom type: lowercase, but keep the X wildcard as 'X'."""
    return "X" if t.upper() == "X" else t.lower()


def parse_prm(path: Path) -> dict:
    txt = Path(path).read_text().splitlines()
    bonds: dict[tuple, tuple] = {}
    angles: dict[tuple, tuple] = {}
    dihe: dict[tuple, tuple] = {}
    impr: dict[tuple, tuple] = {}
    nonb: dict[str, str] = {}          # type -> raw line (verbatim)
    mode = None
    for line in txt:
        s = line.strip()
        if not s or s.startswith(("*", "!")):
            continue
        u = s.upper()
        if u.startswith("BOND"):
            mode = "bond"; continue
        if u.startswith("ANGL"):
            mode = "angl"; continue
        if u.startswith(("DIHE", "DIHEDRAL")):
            mode = "dihe"; continue
        if u.startswith(("IMPR", "IMPHI")):
            mode = "impr"; continue
        if u.startswith("NONBONDED"):
            mode = "nonb"; continue
        if u.startswith(("CUTNB", "HBOND", "END")):
            continue                     # settings / footer, keep mode
        p = s.split()
        if mode == "bond" and len(p) >= 4:
            bonds[(_nt(p[0]), _nt(p[1]))] = (float(p[2]), float(p[3]))
        elif mode == "angl" and len(p) >= 5:
            angles[(_nt(p[0]), _nt(p[1]), _nt(p[2]))] = (float(p[3]), float(p[4]))
        elif mode == "dihe" and len(p) >= 7:
            dihe[tuple(_nt(x) for x in p[:4])] = (float(p[4]), int(p[5]), float(p[6]))
        elif mode == "impr" and len(p) >= 7:
            impr[tuple(_nt(x) for x in p[:4])] = (float(p[4]), int(p[5]), float(p[6]))
        elif mode == "nonb" and len(p) >= 7:
            nonb[_nt(p[0])] = s
    return {"bonds": bonds, "angles": angles, "dihedrals": dihe,
            "impropers": impr, "nonbonded": nonb}


def parse_pdb(path: Path) -> dict[str, tuple[float, float, float]]:
    xyz: dict[str, tuple[float, float, float]] = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            name = line[12:16].strip()
            xyz[name] = (float(line[30:38]), float(line[38:46]),
                         float(line[46:54]))
    return xyz


def element(typ: str) -> str:
    c = typ[0]
    return c.upper() if c.isalpha() else typ[1].upper()


def adjacency(atoms, bonds):
    adj = {a[0]: [] for a in atoms}
    for a, b in bonds:
        adj[a].append(b)
        adj[b].append(a)
    return adj


# ---------------------------------------------------------------------------
# Unique-atom identification (chemical, not graph-isomorphism)
# ---------------------------------------------------------------------------

def identify_methyl(atoms, bonds):
    """Return (methyl_C_name, [H_names]) for the N-methyl in the reference.
    A methyl carbon is an sp3 carbon (type c3) bonded to exactly 3 hydrogens."""
    adj = adjacency(atoms, bonds)
    typ = {a[0]: a[1] for a in atoms}
    for name, t, _ in atoms:
        if t != "c3":
            continue
        h = [nb for nb in adj[name] if element(typ[nb]) == "H"]
        if len(h) == 3:
            return name, h
    raise RuntimeError("methyl carbon not found")


def identify_nh(atoms, bonds):
    """Return (N_name, H_name) for the N-H in the mutant (H type hn)."""
    adj = adjacency(atoms, bonds)
    typ = {a[0]: a[1] for a in atoms}
    for name, t, _ in atoms:
        if t == "hn":
            n_nb = [nb for nb in adj[name] if element(typ[nb]) == "N"]
            assert len(n_nb) == 1, "hn hydrogen must have one N neighbour"
            return n_nb[0], name
    raise RuntimeError("N-H hydrogen (hn) not found")


# ---------------------------------------------------------------------------
# Kabsch + ICP rigid alignment (numpy)
# ---------------------------------------------------------------------------

def kabsch(P: np.ndarray, Q: np.ndarray):
    P = np.asarray(P, float); Q = np.asarray(Q, float)
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = qc - R @ pc
    return R, t


def icp_align(P, Q, elems_P, elems_Q, iters=20):
    """Align points P onto Q (element-constrained nearest-neighbour ICP)."""
    cur = np.asarray(P, float).copy()
    R = np.eye(3); t = np.zeros(3)
    for _ in range(iters):
        pairs = []
        used = set()
        for i in range(len(cur)):
            best = None
            for j in range(len(Q)):
                if elems_P[i] != elems_Q[j] or j in used:
                    continue
                d = np.sum((cur[i] - Q[j]) ** 2)
                if best is None or d < best[0]:
                    best = (d, j)
            if best is not None:
                used.add(best[1])
                pairs.append((i, best[1]))
        if len(pairs) < 3:
            break
        pidx = [i for i, _ in pairs]
        qidx = [j for _, j in pairs]
        Ri, ti = kabsch(P[pidx], Q[qidx])
        cur = P @ Ri.T + ti
        R = Ri @ R
        t = Ri @ t + ti
    return R, t


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ref = parse_rtf(INP / "ref.rtf")
    mut = parse_rtf(INP / "mut.rtf")
    ref_prm = parse_prm(INP / "ref.prm")
    mut_prm = parse_prm(INP / "mut.prm")
    ref_xyz = parse_pdb(INP / "ref.pdb")
    mut_xyz = parse_pdb(INP / "mut.pdb")

    ref_typ = {a[0]: a[1] for a in ref["atoms"]}
    mut_typ = {a[0]: a[1] for a in mut["atoms"]}
    ref_adj = adjacency(ref["atoms"], ref["bonds"])

    # --- unique atoms ---
    methyl_C, methyl_H = identify_methyl(ref["atoms"], ref["bonds"])
    ref_only = [methyl_C] + methyl_H
    nh_N, nh_H = identify_nh(mut["atoms"], mut["bonds"])
    mut_only = [nh_H]
    # the shared (common) nitrogen = the atom the methyl / N-H hang off
    anchor_N = [nb for nb in ref_adj[methyl_C] if element(ref_typ[nb]) == "N"][0]
    print(f"ref-only methyl   : {ref_only}")
    print(f"mut-only N-H      : {mut_only}   (bonded to mut N '{nh_N}')")
    print(f"anchor N (ref)    : {anchor_N}")

    # --- place the N-H via a local N-frame alignment ---
    # align mut anchor N + its two ring neighbours onto ref anchor N + its two
    # ring neighbours (by matching the CHARMM atom types cc<->cc, nc<->nc).
    mut_adj = adjacency(mut["atoms"], mut["bonds"])
    mut_nb = {mut_typ[n]: n for n in mut_adj[nh_N]}
    ref_nb = {ref_typ[n]: n for n in ref_adj[anchor_N]}
    P = np.array([mut_xyz[nh_N], mut_xyz[mut_nb["cc"]], mut_xyz[mut_nb["nc"]]])
    Q = np.array([ref_xyz[anchor_N], ref_xyz[ref_nb["cc"]], ref_xyz[ref_nb["nc"]]])
    R, t = kabsch(P, Q)
    nh_xyz = tuple(R @ np.array(mut_xyz[nh_H]) + t)
    d = float(np.linalg.norm(np.array(nh_xyz) - np.array(ref_xyz[anchor_N])))
    print(f"N-H placed at : ({nh_xyz[0]:.3f}, {nh_xyz[1]:.3f}, {nh_xyz[2]:.3f})"
          f"  (N-H bond {d:.3f} A)")

    # --- Build hybrid rtf ---
    nh_ref_name = "H17"
    nh_type = mut_typ[nh_H]
    mut_atom = {a[0]: a for a in mut["atoms"]}
    nh_chg = mut_atom[nh_H][2]

    masses = dict(ref["mass"])
    if nh_type not in masses:
        masses[nh_type] = mut["mass"].get(nh_type, 1.008)

    lines = ["* Dual-topology hybrid ligand (6I5I: N-CH3 -> N-H)", "*", "   99   1"]
    for i, (t, m) in enumerate(sorted(masses.items()), 1):
        lines.append(f"MASS {i:4d} {t:<5s} {m:10.6f}")
    lines.append("")
    lines.append("RESI UNL  0.000")
    lines.append("GROUP")
    for name, typ, chg in ref["atoms"]:
        lines.append(f"ATOM {name:<5s} {typ:<5s} {chg: .6f}")
    lines.append(f"ATOM {nh_ref_name:<5s} {nh_type:<5s} {nh_chg: .6f}")
    lines.append("")
    # BOND section: ref bonds + the new N-H bond (one record per line)
    for a, b in ref["bonds"]:
        lines.append(f"BOND {a:<5s} {b:<5s}")
    lines.append(f"BOND {anchor_N:<5s} {nh_ref_name:<5s}")

    # ANGL: ref angles + new H17 angles (X-anchor-H17 for each X bonded to anchor)
    new_angles = [(nb, anchor_N, nh_ref_name) for nb in ref_adj[anchor_N]]
    for a, b, c in ref["angles"] + new_angles:
        lines.append(f"ANGL {a:<5s} {b:<5s} {c:<5s}")

    # DIHE: ref dihedrals + new proper dihedrals with H17 as a terminal atom
    new_dihe = []
    for B in ref_adj[anchor_N]:
        for A in ref_adj[B]:
            if A == anchor_N:
                continue
            new_dihe.append((A, B, anchor_N, nh_ref_name))
            new_dihe.append((nh_ref_name, anchor_N, B, A))
    for a, b, c, d in ref["dihedrals"] + new_dihe:
        lines.append(f"DIHE {a:<5s} {b:<5s} {c:<5s} {d:<5s}")

    # IMPH: ref impropers + one N-H planarity improper
    #   mirror the mutant's "C H N N1" improper (types cc hn na nc) onto the
    #   reference naming: (cc-neighbour, H17, anchor_N, nc-neighbour).
    new_impr = []
    if "cc" in ref_nb and "nc" in ref_nb:
        new_impr.append((ref_nb["cc"], nh_ref_name, anchor_N, ref_nb["nc"]))
    for a, b, c, d in ref["impropers"] + new_impr:
        lines.append(f"IMPH {a:<5s} {b:<5s} {c:<5s} {d:<5s}")
    lines.append("END")
    lines.append("")
    (OUT / "hybrid.rtf").write_text("\n".join(lines))
    print("wrote hybrid.rtf")

    # --- Build hybrid prm (ref params + N-H terms from mut) ---
    pl = ["* Dual-topology hybrid parameters", "*"]
    pl.append("BOND")
    for k, v in ref_prm["bonds"].items():
        pl.append(f"{k[0]:<5s} {k[1]:<5s} {v[0]:9.2f} {v[1]:8.3f}")
    for k, v in mut_prm["bonds"].items():
        if nh_type in k and k not in ref_prm["bonds"]:
            pl.append(f"{k[0]:<5s} {k[1]:<5s} {v[0]:9.2f} {v[1]:8.3f}")
    pl.append("")
    pl.append("ANGLE")
    for k, v in ref_prm["angles"].items():
        pl.append(f"{k[0]:<5s} {k[1]:<5s} {k[2]:<5s} {v[0]:9.2f} {v[1]:8.3f}")
    for k, v in mut_prm["angles"].items():
        if nh_type in k and k not in ref_prm["angles"]:
            pl.append(f"{k[0]:<5s} {k[1]:<5s} {k[2]:<5s} {v[0]:9.2f} {v[1]:8.3f}")
    # cross angle c3-na-hn (methyl C - N - H): both endpoints off, small in between
    pl.append(f"{'c3':<5s} {'na':<5s} {nh_type:<5s} {46.8:9.2f} {125.5:8.3f}   ! approx cross term")
    pl.append("")
    pl.append("DIHEDRAL")
    for k, v in ref_prm["dihedrals"].items():
        pl.append(f"{k[0]:<4s} {k[1]:<4s} {k[2]:<4s} {k[3]:<4s} {v[0]:8.3f} {v[1]:3d} {v[2]:8.1f}")
    for k, v in mut_prm["dihedrals"].items():
        if nh_type in k and k not in ref_prm["dihedrals"]:
            pl.append(f"{k[0]:<4s} {k[1]:<4s} {k[2]:<4s} {k[3]:<4s} {v[0]:8.3f} {v[1]:3d} {v[2]:8.1f}")
    pl.append("")
    pl.append("IMPROPER")
    for k, v in ref_prm["impropers"].items():
        pl.append(f"{k[0]:<4s} {k[1]:<4s} {k[2]:<4s} {k[3]:<4s} {v[0]:8.3f} {v[1]:3d} {v[2]:8.1f}")
    for k, v in mut_prm["impropers"].items():
        if nh_type in k and k not in ref_prm["impropers"]:
            pl.append(f"{k[0]:<4s} {k[1]:<4s} {k[2]:<4s} {k[3]:<4s} {v[0]:8.3f} {v[1]:3d} {v[2]:8.1f}")
    pl.append("")
    pl.append("NONBONDED  NBXMOD 5  GROUP SWITCH CDIEL -")
    pl.append("CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4")
    pl.append("!                Emin     Rmin/2              Emin/2     Rmin  (for 1-4's)")
    pl.append("!             (kcal/mol)    (A)")
    for t, raw in sorted(ref_prm["nonbonded"].items()):
        pl.append(raw)
    if nh_type not in ref_prm["nonbonded"] and nh_type in mut_prm["nonbonded"]:
        pl.append(mut_prm["nonbonded"][nh_type])
    pl.append("")
    (OUT / "hybrid.prm").write_text("\n".join(pl))
    print("wrote hybrid.prm")

    # --- Build hybrid pdb with B-factors ---
    bfac = {n: -1.0 for n in ref_only}
    bfac[nh_ref_name] = 1.0
    pdbl = []
    serial = 0
    for name, typ, chg in ref["atoms"]:
        serial += 1
        x, y, z = ref_xyz[name]
        b = bfac.get(name, 0.0)
        pdbl.append(f"HETATM{serial:5d} {name:<4s} UNL Z   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{b:6.2f}          {element(typ):>2s}")
    serial += 1
    x, y, z = nh_xyz
    pdbl.append(f"HETATM{serial:5d} {nh_ref_name:<4s} UNL Z   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{1.0:6.2f}          {element(nh_type):>2s}")
    pdbl.append("END")
    (OUT / "hybrid.pdb").write_text("\n".join(pdbl) + "\n")
    print("wrote hybrid.pdb")
    print(f"\nHybrid ligand: {len(ref['atoms']) + 1} atoms "
          f"({len(ref_only)} vanish, {len(mut_only)} appear, "
          f"{len(ref['atoms']) - len(ref_only)} common)")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    main()
