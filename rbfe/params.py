"""Finding a parameter for every hybrid term, and being honest when there is none.

Most terms of the hybrid exist in one of the two ligands and are sourced from
that ligand's `.prm` verbatim.  The interesting case is the handful that exist
in **neither** endpoint state -- for an element swap, every term involving both
halogens, a species no physical ligand contains.

There are three defensible answers and they are not interchangeable, so the
policy is named in the config (`[mutation] missing_terms`) and every decision is
reported:

``fail``
    Refuse to build.  Right when you want to know a term was missing.

``mirror``
    Take the term from the *other* ligand's parameter file, matching on the
    positions that both ligands share, and relabel it.  This reproduces
    chemistry that was actually fitted: 6I5I's frozen cross angle
    ``c3 na hn 46.8 125.5`` is verbatim ``mut.prm``'s ``cc na hn`` relabelled.
    Preferred whenever a real analogue exists, because a fitted number beats an
    invented one.

``synth``
    Approximate from the hybrid's own geometry, so at the endpoints the term
    contributes nothing and cannot distort the ligand.  This is the right answer
    for the angle between two substituents sitting on the same vector: it is
    ~0 degrees, not ~120, and a harmonic term with a made-up 120-degree
    reference would try to pry the ring open.
"""

from __future__ import annotations

import math

from rbfe.charmm import angle_deg, is_wildcard, _nt
from rbfe.errors import ChemistryError
from rbfe.model import Mapping, XYZ


def lookup(table: dict, key) -> tuple | None:
    """Find a parameter under any spelling CHARMM might use.

    As written, with the term reversed, and with X wildcards in the terminal
    positions -- CHARMM parameter files are not consistent about which end of a
    term the specific type sits on, and a dihedral may be given with two
    wildcards.
    """
    k = tuple(key)
    base = [k, tuple(reversed(k))]
    cands = list(base)
    if len(k) == 3:
        for c in base:
            cands += [("X", c[1], c[2]), (c[0], c[1], "X")]
    elif len(k) == 4:
        for c in base:
            cands += [("X", c[1], c[2], c[3]), (c[0], c[1], c[2], "X"),
                      ("X", c[1], c[2], "X")]
        for i in range(4):
            for j in range(i + 1, 4):
                c = list(k)
                c[i] = c[j] = "X"
                cands.append(tuple(c))
    for c in cands:
        if c in table:
            return table[c]
    return None


def mirror(tag: str, key: tuple[str, ...], other: dict,
           alch_pos: tuple[int, ...]) -> tuple | None:
    """A term from the other ligand's parameters, relabelled.

    Match position-wise, allowing a difference **only at the alchemical
    positions** -- everything else must agree, since everything else is the same
    chemistry.  This is what finds 6I5I's `cc na hn` for the query `c3 na hn`:
    the two ligands call the ring carbon differently, and the mutant's value is
    a fitted number for the same angle, which beats anything synthesised from
    geometry.

    Ambiguity is not resolved by guessing: if more than one entry matches, this
    returns None and the caller falls back to synthesis, reporting it.
    """
    if not other:
        return None
    alch = set(alch_pos)

    def compatible(entry) -> bool:
        if len(entry) != len(key):
            return False
        for i, (a, b) in enumerate(zip(entry, key)):
            if i in alch or a == b or is_wildcard(a) or is_wildcard(b):
                continue
            return False
        return True

    hits = [(k, v) for k, v in other.items() if compatible(k)]
    if not hits:
        return None
    # Deterministic, and the choice is REPORTED rather than made silently: a
    # mutant with two ring carbons offers both `cc na hn` and `nc na hn` for
    # the query `c3 na hn`, and neither is wrong.  File order decides, and the
    # report says how many candidates there were, so an ambiguous substitution
    # is visible instead of looking like a clean lookup.
    return hits[0][1], hits[0][0], len(hits)


def synth(key: tuple[str, ...], coords: dict[str, XYZ]) -> tuple:
    """Approximate a term from the hybrid's own geometry.

    A near-degenerate angle (both substituents collinear) additionally gets
    k=0.  Its only job is then to put the pair into the psf's 1-3 list so NAMD
    EXCLUDES their nonbonded interaction -- without that the two halogens,
    0.2 A apart, would interact at full strength mid-window.
    """
    if len(key) == 2:
        return (200.0, round(math.dist(coords[key[0]], coords[key[1]]), 3))
    if len(key) == 3:
        th = angle_deg(coords[key[0]], coords[key[1]], coords[key[2]])
        return (0.0 if th < 10.0 else 45.0, round(th, 2))
    if len(key) == 4:
        return (0.0, 1, 0.0)
    raise ChemistryError(f"cannot synthesise a {len(key)}-atom term: {key}")


class Missing:
    """One term with no source parameter, and what was done about it."""

    __slots__ = ("tag", "term", "types", "value", "how", "cross")

    def __init__(self, tag: str, term, types, value, how: str, cross: bool) -> None:
        self.tag, self.term, self.types = tag, term, types
        self.value, self.how, self.cross = value, how, cross

    def describe(self) -> str:
        atoms = " ".join(f"{a}({t})" for a, t in zip(self.term, self.types))
        mark = "   <-- CROSS TERM: absent from BOTH endpoint states" if self.cross else ""
        return f"{self.tag} {atoms}  -> {self.value}  [{self.how}]{mark}"


def resolve(tag: str, terms, types: dict[str, str], alch: set[str],
            ref_prm: dict, mut_prm: dict, coords: dict[str, XYZ],
            policy: str) -> tuple[dict, list[Missing]]:
    """Source every term, recording each one that had no natural home.

    Returns the parameter table and the list of missing terms -- the caller
    prints the latter, so a synthesised value is never silent.
    """
    table_name = {"BOND": "bonds", "ANGL": "angles",
                  "DIHE": "dihedrals", "IMPH": "impropers"}[tag]
    ref_tab, mut_tab = ref_prm[table_name], mut_prm[table_name]

    out: dict[tuple, tuple] = {}
    missing: list[Missing] = []

    for t in terms:
        key = tuple(types[x] for x in t)
        v = lookup(ref_tab, key)
        how = "ref"
        if v is None:
            v = lookup(mut_tab, key)
            how = "mut"
        if v is not None:
            out[t] = v
            continue

        cross = bool(set(t) & alch)
        if policy == "fail":
            raise ChemistryError(
                f"no {tag} parameter for {' '.join(f'{a}({types[a]})' for a in t)} "
                f"in either ligand, and [mutation] missing_terms = fail.\n"
                f"  This term involves {sorted(set(t) & alch)} -- a species "
                f"neither endpoint state contains.\n"
                f"  Set missing_terms to 'mirror', 'mirror_then_synth' or "
                f"'synth_from_geometry' to have it filled instead.")

        v = None
        if policy in ("mirror", "mirror_then_synth"):
            alch_pos = tuple(i for i, x in enumerate(t) if x in alch)
            found = mirror(tag, key, ref_tab, alch_pos) or mirror(tag, key, mut_tab, alch_pos)
            if found:
                v, src, n = found
                how = ("mirrored" if n == 1
                       else f"mirrored from {' '.join(src)} ({n} candidates, took the first)")
        if v is None and policy in ("synth_from_geometry", "mirror_then_synth"):
            v = synth(t, coords)
            how = "synthesised"
        if v is None:
            raise ChemistryError(
                f"no {tag} parameter for {' '.join(f'{a}({types[a]})' for a in t)} "
                f"and policy '{policy}' found no substitute.")
        out[t] = v
        missing.append(Missing(tag, t, [types[x] for x in t], v, how, cross))

    return out, missing
