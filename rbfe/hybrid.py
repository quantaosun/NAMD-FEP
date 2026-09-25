"""Build the dual-topology hybrid ligand: hybrid.{pdb,rtf,prm}.

The chemistry lives in a `MutationStrategy` (rbfe/mutations/); this module is
the plumbing every mutation needs, and contains no chemistry of its own:

  1. map the two ligands to one hybrid            <- strategy
  2. close the charges so both endpoints are right
  3. enumerate the hybrid's connectivity           <- rbfe.topology
  4. source a parameter for every term             <- rbfe.params
  5. emit the .rtf / .prm / .pdb

Step 2 is the one worth stating plainly.  The two ligands' charge sets disagree
over the shared core (they were fitted separately), so the appearing atoms
absorb the difference.  Without it the hybrid would not sum to the net charge of
either physical ligand and lambda=0 and lambda=1 would describe something that
does not exist.
"""

from __future__ import annotations

from pathlib import Path

from rbfe import params as P
from rbfe import similarity as S
from rbfe import topology as T
from rbfe.charmm import element, load_ligand
from rbfe.errors import ChemistryError
from rbfe.model import Atom, Ligand, Mapping

# Namespaces.  NAMD matches atom types with strcasecmp and, on a duplicate, the
# LAST file read wins (src/Parameters.C: add_to_vdw_tree / add_bond).  GAFF
# types are lowercase and CHARMM36's are uppercase, but that does NOT keep them
# apart: a ligand's `ca c o ha hn na nb` collide with the protein's CA C O HA HN
# and with the nucleic NA NB.  hybrid.prm is listed last in the .namd files, so
# without a prefix the LIGAND silently overwrites the PROTEIN's parameters and
# NAMD reports only "Warning: DUPLICATE vdW ENTRY FOR CA".
#
# `_prefix` returns the string to prepend: "" for the declared compatibility
# mode, otherwise the configured prefix.
RAW = "raw"

FILLER = {"mirror": "mirrored", "synth_from_geometry": "synthesised",
          "mirror_then_synth": "mirrored or synthesised"}


class Builder:
    def __init__(self, ref: Ligand, mut: Ligand, mapping: Mapping,
                 resname: str, title: str, type_prefix: str,
                 rtf_order: str, prm_order: str,
                 missing_terms: str) -> None:
        self.ref, self.mut, self.m = ref, mut, mapping
        self.resname, self.title = resname, title
        self.rtf_order = rtf_order
        self.prm_order = prm_order
        self.missing_terms = missing_terms
        self.prefix = "" if type_prefix == RAW else type_prefix

        # MASS lines from both ligands; the reference wins on a clash.
        self.masses: dict[str, float] = {t: a.mass for t, a in
                                         ((a.type, a) for a in ref.atoms.values())}
        for a in mut.atoms.values():
            self.masses.setdefault(a.type, a.mass)

    # -- helpers -------------------------------------------------------
    def px(self, t: str) -> str:
        """Prefix a type, leaving the X wildcard alone (prefixing it would stop
        it matching anything)."""
        return t if t.upper() == "X" else self.prefix + t

    def _mass(self, typ: str) -> float:
        """A mass for the type, falling back to the other ligand's."""
        if typ in self.masses:
            return self.masses[typ]
        raise ChemistryError(
            f"no MASS line for atom type '{typ}' in either ligand's .rtf; "
            f"the hybrid would be written with mass 0 and NAMD would treat the "
            f"atoms as massless.")

    def _atom(self, name: str) -> Atom:
        return Atom(name=name, type=self.m.types[name],
                    charge=self.m.charges[name], mass=self._mass(self.m.types[name]))

    # -- 2. charge closure ---------------------------------------------
    def close_charges(self) -> dict:
        m, ref, mut = self.m, self.ref, self.mut
        ref_total = ref.total_charge
        mut_total = mut.total_charge
        core = sum(m.charges[n] for n in m.common)
        vanish = sum(m.charges[n] for n in m.vanish)
        # The strategy has already put each appearing atom's *native* charge in
        # the mapping: for element_swap the mutant's own value under the same
        # name, for atom_addition the mutant's value under the newly invented
        # name -- which by construction is not a key in mut.atoms.  Read it from
        # the mapping, not from the mutant.
        native = sum(m.charges[n] for n in m.appear)
        if not m.appear:
            raise ChemistryError(
                "nothing is appearing in this mutation, so there is no atom to "
                "absorb the charge difference between the two ligands. A pure "
                "deletion is not supported; declare the mutation the other way "
                "round (reference and mutant swapped) or add an appearing atom.")
        corr = (mut_total - core - native) / len(m.appear)
        for n in m.appear:
            m.charges[n] += corr

        appear = sum(m.charges[n] for n in m.appear)
        lam0, lam1 = core + vanish, core + appear

        for label, total in (("reference", ref_total), ("mutant", mut_total)):
            if abs(total - round(total)) > 0.02:
                raise ChemistryError(
                    f"the {label} ligand's charges sum to {total:+.6f}, which is "
                    f"not within 0.02 of an integer. Charges from acpype/AM1-BCC "
                    f"are rounded, so this usually means the .rtf is truncated "
                    f"or was edited by hand.")

        if abs(lam0 - ref_total) > 1e-3 or abs(lam1 - mut_total) > 1e-3:
            raise ChemistryError(
                f"charge closure failed: lambda=0 {lam0:+.6f} (want {ref_total:+.6f}), "
                f"lambda=1 {lam1:+.6f} (want {mut_total:+.6f})")

        return {"ref_total": ref_total, "mut_total": mut_total, "core": core,
                "vanish": vanish, "native": native, "appear": appear,
                "corr": corr, "lam0": lam0, "lam1": lam1}

    # -- 5. emission ---------------------------------------------------
    @staticmethod
    def _key(t) -> tuple:
        """Comparison form of a term.

        A BOND is an unordered pair; a DIHE reads the same backwards; an ANGL
        and an IMPH are written vertex-second and only normalise as themselves.
        """
        t = tuple(t)
        if len(t) == 2:
            return tuple(sorted(t))
        if len(t) in (3, 4):
            # ANGL and DIHE both read identically backwards; IMPH does not
            # (its first three atoms define a plane and the fourth is the
            # out-of-plane one).
            return min(t, tuple(reversed(t))) if len(t) == 4 else min(t, tuple(reversed(t)))
        return t

    def _order(self, table: str, terms) -> list:
        """Order one term set for the .rtf, per `[hybrid] rtf_order`.

        `sorted`       full enumeration of the merged connectivity, sorted
                       (new systems, and 4YLJ).
        `source_then_new`
                       the reference's own .rtf entries in file order, then the
                       terms the perturbation adds (6I5I's builder).

        Subtracting the reference's *listed* topology would be wrong: a CHARMM
        `.rtf` does not spell out every angle its bonds imply, so about twenty
        pre-existing angles would be misreported as new and written twice.  The
        comparison is against what the reference already implies.
        """
        if self.rtf_order == "sorted":
            return sorted(terms)

        m = self.m
        universe = set(m.types)
        listed, seen = [], set()
        for t in getattr(self.ref, table):
            k = self._key(t)
            if all(x in universe for x in k) and k not in seen:
                seen.add(k)
                listed.append(tuple(t))

        # The reference ALONE -- deliberately without mapping.extra_bonds: a
        # bond the perturbation invents is not something the reference already
        # implied, and folding it in here would silently drop the very terms
        # this function exists to append.
        ref_bonds = {tuple(sorted(b)) for b in self.ref.bonds
                     if all(x in universe for x in b)}
        base = {self._key(x) for x in
                T.implied(ref_bonds, set(m.vanish), universe)[table]}
        already = {self._key(t) for t in listed}
        return listed + sorted(t for t in terms
                               if self._key(t) not in base
                               and self._key(t) not in already)

    def render_rtf(self, topo: T.HybridTopology) -> str:
        m = self.m
        used = sorted({m.types[n] for n in m.types})
        lines = [f"* {self.title}", "*", "   99   1"]
        lines += [f"MASS {i:4d} {self.px(t):<5s} {self._mass(t):10.6f}"
                  for i, t in enumerate(used, 1)]
        lines += ["", f"RESI {self.resname} {sum(m.charges.values()): .6f}", "GROUP"]
        lines += [f"ATOM {n:<5s} {self.px(m.types[n]):<5s} {m.charges[n]: .6f}"
                  for n in m.order]
        lines += [""]
        self.emitted = {}
        lines += [f"BOND {a:<5s} {b:<5s}" for a, b in self._order("bonds", topo.bonds)]
        lines += [f"ANGL {a:<5s} {b:<5s} {c:<5s}"
                  for a, b, c in self._order("angles", topo.angles)]
        lines += [f"DIHE {a:<5s} {b:<5s} {c:<5s} {d:<5s}"
                  for a, b, c, d in self._order("dihedrals", topo.dihedrals)]
        lines += [f"IMPH {a:<5s} {b:<5s} {c:<5s} {d:<5s}"
                  for a, b, c, d in self._order("impropers", topo.impropers)]
        lines += ["END", ""]
        return "\n".join(lines)

    @staticmethod
    def _f3(key, v) -> str:
        return " ".join(f"{x:<5s}" for x in key) + f" {v[0]:9.2f} {v[1]:8.3f}"

    @staticmethod
    def _f4(key, v) -> str:
        return " ".join(f"{x:<4s}" for x in key) + f" {v[0]:8.3f} {v[1]:3d} {v[2]:8.1f}"

    def render_prm(self, tables: dict, missing: list[P.Missing], emitted: dict) -> str:
        m = self.m
        ref_prm, mut_prm = tables["ref"], tables["mut"]
        p3 = lambda k, v: self._f3([self.px(x) for x in k], v)
        p4 = lambda k, v: self._f4([self.px(x) for x in k], v)

        used = {m.types[n] for n in m.types}
        # The type-tuples the hybrid's own topology actually spells out.  The
        # mutant side is filtered to exactly these: not "every type the hybrid
        # contains" (still too many -- a mutant parameter file describes a
        # molecule the hybrid is not) and not "terms mentioning the new atom"
        # (too few in general, and it only worked for 6I5I by accident).
        # Exactly the terms the .rtf lists -- not the full enumeration, which
        # for `source_then_new` is a superset.  A parameter is needed iff the
        # topology references it, so the topology is the thing to compare with.
        spelled = {tag: {tuple(m.types[x] for x in t) for t in terms}
                   for tag, terms in emitted.items()}

        def keep(rows, table):
            """Only mutant parameters the hybrid can actually reference.

            The REFERENCE side is copied verbatim -- it is the ligand the
            hybrid is built around, and its file is reproduced as-is.  The
            mutant side is filtered, because a mutant parameter file describes
            a molecule the hybrid is not: carrying its entries for atom types
            the hybrid does not contain is not merely untidy, in `raw` mode
            those extra GAFF types collide case-insensitively with the
            protein's and NAMD takes the last duplicate, so padding the file
            with them *adds* collisions rather than being harmless.
            """
            # Normalised comparison: a parameter file may key an angle as
            # `nc na hn` where the topology spells it `hn na nc`, and a bond in
            # either order.  Comparing raw tuples would silently drop the term
            # the new atom needs.
            want = {self._key(x) for x in spelled[table]}
            return {k: v for k, v in rows.items() if self._key(k) in want}

        def block(heading: str, table: str, fmt):
            r, mu = ref_prm[table], keep(mut_prm[table], table)
            if self.prm_order == "sorted":
                merged = dict(mu)
                merged.update(r)          # the reference wins on a clash
                rows = [fmt(k, v) for k, v in sorted(merged.items())]
            else:
                # "already in the reference" must be a normalised test: the
                # reference may key the same angle the other way round, and a
                # raw dict lookup would then emit the mutant's copy as well.
                have = {self._key(k) for k in r}
                rows = [fmt(k, v) for k, v in r.items()]
                rows += [fmt(k, v) for k, v in mu.items() if self._key(k) not in have]
            return ["", heading] + rows

        pl = [f"* {self.title} -- parameters", "*", "BOND"]
        if self.prm_order == "sorted":
            merged_b = keep(mut_prm["bonds"], "bonds"); merged_b.update(ref_prm["bonds"])
            pl += [p3(k, v) for k, v in sorted(merged_b.items())]
        else:
            rb, mb = ref_prm["bonds"], keep(mut_prm["bonds"], "bonds")
            have = {self._key(k) for k in rb}
            pl += [p3(k, v) for k, v in rb.items()]
            pl += [p3(k, v) for k, v in mb.items() if self._key(k) not in have]
        pl += block("ANGLE", "angles", p3)

        # Every synthesised term is written, in all four categories.  An earlier
        # version emitted only the ANGL ones, so a missing bond or dihedral was
        # *reported* as approximated and then silently absent from the file.
        synth_lines = {"BOND": [], "ANGL": [], "DIHE": [], "IMPH": []}
        for miss in missing:
            key = tuple(self.m.types[x] for x in miss.term)
            synth_lines[miss.tag].append(
                f"{(p3 if miss.tag in ('BOND', 'ANGL') else p4)(key, miss.value)}"
                f"   ! approx cross term")
        if synth_lines["BOND"]:
            pl += ["", "! synthesised bonds"] + synth_lines["BOND"]
        pl += synth_lines["ANGL"]

        pl += block("DIHEDRAL", "dihedrals", p4)
        pl += synth_lines["DIHE"]
        pl += block("IMPROPER", "impropers", p4)
        pl += synth_lines["IMPH"]

        pl += ["", "NONBONDED  NBXMOD 5  GROUP SWITCH CDIEL -",
               "CUTNB 14.0  CTOFNB 12.0  CTONNB 10.0  EPS 1.0  E14FAC 0.83333333  WMIN 1.4",
               "!                Emin     Rmin/2              Emin/2     Rmin  (for 1-4's)",
               "!             (kcal/mol)    (A)"]
        seen = set()
        for t in sorted({self.m.types[n] for n in self.m.types}):
            raw = ref_prm["nonbonded"].get(t) or mut_prm["nonbonded"].get(t)
            if raw and t not in seen:
                head, rest = raw.split(None, 1)
                pl.append(f"{self.px(head):<5s} {rest}")
                seen.add(t)
        pl.append("")
        return "\n".join(pl)

    def render_pdb(self) -> str:
        m = self.m
        bfac = {n: 0.0 for n in m.types}
        for n in m.vanish:
            bfac[n] = -1.0
        for n in m.appear:
            bfac[n] = 1.0
        self._bfac = bfac
        rows = []
        for i, n in enumerate(m.order, 1):
            x, y, z = m.coords[n]
            rows.append(
                f"HETATM{i:5d} {n:<4s} {self.resname:<3s} Z   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.0:6.2f}{bfac[n]:6.2f}          {element(self._atom(n)):>2s}")
        rows.append("END")
        return "\n".join(rows) + "\n"


def build(ref_stem: Path, mut_stem: Path, out: Path, *, resname: str, title: str,
          type_prefix: str, rtf_order: str, prm_order: str, missing_terms: str,
          strategy, verbose: bool = True,
          similarity_threshold: float | None = S.DEFAULT_THRESHOLD) -> dict:
    """Run the whole builder.  Returns a report dict; writes three files.

    `similarity_threshold` screens the ligand pair before anything is emitted --
    RBFE is only meaningful between similar ligands, and nothing downstream would
    notice if it were not.  It warns and continues; see `rbfe.similarity`.  Pass
    None to skip the screen entirely (the verifier does, via `verbose=False`).
    """
    out = Path(out)
    ref, mut = load_ligand(ref_stem), load_ligand(mut_stem)
    p = print if verbose else (lambda *a, **k: None)

    if verbose and similarity_threshold is not None:
        S.check(ref, mut, threshold=similarity_threshold, printer=p)

    mapping = strategy.map(ref, mut)
    p(f"atom mapping         : {len(mapping.common)} common, "
      f"{len(mapping.vanish)} vanish {mapping.vanish}, "
      f"{len(mapping.appear)} appear {mapping.appear}")
    for note in mapping.provenance:
        p(f"                       {note}")

    b = Builder(ref, mut, mapping, resname, title, type_prefix, rtf_order,
                prm_order, missing_terms)
    closure = b.close_charges()
    p(f"charge closure       : ref {closure['ref_total']:+.6f}  "
      f"mut {closure['mut_total']:+.6f}")
    p(f"                       core {closure['core']:+.6f}  "
      f"vanish {closure['vanish']:+.6f}")
    p(f"                       appearing {closure['native']:+.6f} -> "
      f"{closure['appear']:+.6f} (correction {closure['corr']:+.6f} "
      f"spread over {mapping.appear})")
    p(f"                       lambda=0 {closure['lam0']:+.6f}   "
      f"lambda=1 {closure['lam1']:+.6f}")

    topo = T.build(mapping, ref, mut)
    alch = mapping.alchemical
    p(f"hybrid               : {len(mapping.types)} atoms, {len(topo.bonds)} bonds, "
      f"alchemical {sorted(alch)}")

    if type_prefix == RAW:
        colliding = _collisions(topo, mapping, ref)
        if colliding:
            p(f"  WARNING: {len(colliding)} atom type(s) collide case-insensitively "
              f"with the protein's: {', '.join(sorted(colliding))}")
            p(f"           NAMD takes the last duplicate, and hybrid.prm is listed "
              f"last -- the ligand will override the PROTEIN here. Set "
              f"[hybrid] type_prefix to a namespace (e.g. L) to avoid it.")

    ref_prm = _load_prm(ref_stem)
    mut_prm = _load_prm(mut_stem)
    tables = {"ref": ref_prm, "mut": mut_prm}
    sourced, missing = {}, []
    for tag, terms, key in (("BOND", topo.bonds, "bonds"),
                            ("ANGL", topo.angles, "angles"),
                            ("DIHE", topo.dihedrals, "dihedrals"),
                            ("IMPH", topo.impropers, "impropers")):
        tbl, miss = P.resolve(tag, sorted(terms), mapping.types, alch,
                              ref_prm, mut_prm, mapping.coords, missing_terms)
        sourced[tag] = tbl
        missing.extend(miss)

    p(f"parameters           : {len(topo.bonds)} bonds, {len(topo.angles)} angles, "
      f"{len(topo.dihedrals)} dihedrals, {len(topo.impropers)} impropers")
    cross = [x for x in missing if x.cross]
    p(f"                       {len(missing)} term(s) with no source parameter "
      f"({len(cross)} alchemical)")
    for miss in missing:
        p(f"                         {miss.describe()}")

    out.mkdir(parents=True, exist_ok=True)
    for label, val in (("rtf_order", rtf_order), ("prm_order", prm_order)):
        if val not in ("sorted", "source_then_new"):
            raise ChemistryError(f"unknown [hybrid] {label} '{val}'")
    (out / "hybrid.rtf").write_text(b.render_rtf(topo))
    emitted = {"bonds": b._order("bonds", topo.bonds),
               "angles": b._order("angles", topo.angles),
               "dihedrals": b._order("dihedrals", topo.dihedrals),
               "impropers": b._order("impropers", topo.impropers)}
    (out / "hybrid.prm").write_text(b.render_prm(tables, missing, emitted))
    (out / "hybrid.pdb").write_text(b.render_pdb())
    p(f"wrote {out}/hybrid.rtf, hybrid.prm, hybrid.pdb")

    if verbose:
        _geometry_check(b, topo, mapping)
    return {"mapping": mapping, "topology": topo, "missing": missing,
            "closure": closure, "builder": b}


def _load_prm(stem: Path) -> dict:
    from rbfe.charmm import parse_prm
    return parse_prm(Path(stem).with_suffix(".prm"))


def _collisions(topo, mapping, ref: Ligand) -> set[str]:
    """Types that would collide with the protein's, for the raw-prefix warning."""
    protein = {"CA", "C", "O", "HA", "HN", "N", "CB", "CG", "CD", "CE", "CZ",
               "NA", "NB", "OS", "OH", "CT", "HP", "H"}
    return {t for t in {mapping.types[n] for n in mapping.types}
            if t.upper() in protein and t != t.upper()}


def _geometry_check(b: Builder, topo, mapping) -> None:
    """Print the distances and angles psfgen is about to receive."""
    from rbfe.charmm import dist
    m = mapping
    print("geometry check:")
    for n in sorted(m.alchemical):
        nb = sorted(topo.adjacency[n])
        if not nb:
            continue
        d = nb[0]
        print(f"  {n:<4s}({m.types[n]:<3s}) - {d}({m.types[d]}) = "
              f"{dist(m.coords[n], m.coords[d]):.3f} A   B-factor {b._bfac[n]:+.0f}")
    for n in m.types:
        hs = sorted(set(topo.adjacency[n]) & m.alchemical)
        if len(hs) == 2:
            from rbfe.charmm import angle_deg
            ds = ", ".join(f"{h}={dist(m.coords[n], m.coords[h]):.3f} A" for h in hs)
            print(f"  both alchemical on {n}({m.types[n]}): {ds}; "
                  f"X-{n}-X = {angle_deg(m.coords[hs[0]], m.coords[n], m.coords[hs[1]]):.1f} deg")
    print(f"\nHybrid ligand: {len(m.types)} atoms "
          f"({len(m.vanish)} vanish, {len(m.appear)} appear, {len(m.common)} common)")
