"""Connectivity of the hybrid ligand: which bonds, angles, dihedrals, impropers.

Ported from `4YLJ/prepare_hybrid.py` step 4, which is the *generic* half of the
builder -- it derives every term implied by the merged bond list instead of
copying the reference topology and appending hand-written terms.

That this reproduces 6I5I's hand-written new angles and dihedrals exactly is
what makes `atom_addition` a small strategy rather than a second builder: the
only term it needs beyond the generic enumeration is the planarity improper at
the new atom, which it supplies explicitly.

The rules, in one place so they can be argued with:

* bonds      union of both ligands', remapped to hybrid names
* angles     every pair of neighbours of every atom
* dihedrals  those in the inputs, plus any spanning a bond where at least one
             end is alchemical -- a dihedral that does not touch the
             perturbation has no reason to be invented
* impropers  those in the inputs, plus a planarity term at any genuine
             3-coordinate sp2 centre that involves an alchemical atom
"""

from __future__ import annotations

import itertools

from rbfe.model import Ligand, Mapping, HybridTopology


def remap(terms, universe: set[str]) -> list[tuple[str, ...]]:
    """Keep only terms whose every atom survives into the hybrid."""
    return [tuple(t) for t in terms if all(x in universe for x in t)]


def implied(bonds, alch: set[str], universe: set[str]) -> dict[str, set]:
    """The terms a bond set implies, using the same rules as `build`.

    Exists so `source_then_new` can ask "what does this term set add?" -- a
    reference `.rtf` lists most of its implied angles but not all of them, so
    subtracting the *listed* terms would misreport ~20 pre-existing angles as
    new.  The comparison has to be against what the reference already implies.
    """
    adj: dict[str, set[str]] = {n: set() for n in universe}
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)

    angles: set[tuple[str, str, str]] = set()
    for centre, neighbours in adj.items():
        for a, c in itertools.combinations(sorted(neighbours), 2):
            angles.add((a, centre, c))

    dihedrals: set[tuple[str, str, str, str]] = set()
    for b, c in bonds:
        for a in adj[b] - {c}:
            for d in adj[c] - {b}:
                if a != d and (a in alch or d in alch):
                    dihedrals.add((a, b, c, d))

    impropers: set[tuple[str, str, str, str]] = set()
    for centre, neighbours in adj.items():
        if len(neighbours) == 3 and (centre in alch or any(x in alch for x in neighbours)):
            for a, c, d in itertools.combinations(sorted(neighbours), 3):
                impropers.add((a, centre, c, d))

    return {"bonds": set(bonds), "angles": angles,
            "dihedrals": dihedrals, "impropers": impropers}


def build(mapping: Mapping, ref: Ligand, mut: Ligand) -> HybridTopology:
    """The hybrid's connectivity, from the mapping and the two input ligands."""
    universe = set(mapping.types)
    alch = mapping.alchemical

    other = remap(mut.bonds, universe) if mapping.share_names else []
    bonds: set[tuple[str, str]] = set()
    for t in remap(ref.bonds, universe) + other:
        bonds.add(tuple(sorted(t)))
    bonds |= {tuple(sorted(b)) for b in mapping.extra_bonds}

    adjacency: dict[str, set[str]] = {n: set() for n in universe}
    for a, b in bonds:
        adjacency[a].add(b)
        adjacency[b].add(a)

    angles: set[tuple[str, str, str]] = set()
    for centre, neighbours in adjacency.items():
        for a, c in itertools.combinations(sorted(neighbours), 2):
            angles.add((a, centre, c))

    dihedrals = {tuple(d) for d in remap(ref.dihedrals, universe)}
    if mapping.share_names:
        dihedrals |= {tuple(d) for d in remap(mut.dihedrals, universe)}
    for b, c in bonds:
        for a in adjacency[b] - {c}:
            for d in adjacency[c] - {b}:
                if a != d and (a in alch or d in alch):
                    dihedrals.add((a, b, c, d))

    impropers = {tuple(d) for d in remap(ref.impropers, universe)}
    if mapping.share_names:
        impropers |= {tuple(d) for d in remap(mut.impropers, universe)}
    for centre, neighbours in adjacency.items():
        # Only a genuine 3-coordinate sp2 centre gets a synthesised planarity
        # term.  The carbon carrying BOTH halogens has four neighbours and is
        # not a normal sp2 centre -- the input impropers already cover it, and
        # inventing more would constrain a species that never exists.
        if len(neighbours) == 3 and (centre in alch or any(x in alch for x in neighbours)):
            for a, c, d in itertools.combinations(sorted(neighbours), 3):
                impropers.add((a, centre, c, d))

    # Terms a strategy adds beyond the generic enumeration -- for `atom_addition`
    # this is the planarity improper at the newly placed atom, which the mutant's
    # own topology implies but the merged bond list cannot see.
    impropers |= set(mapping.extra_impropers)

    return HybridTopology(bonds=bonds, angles=angles, dihedrals=dihedrals,
                          impropers=impropers, adjacency=adjacency)
