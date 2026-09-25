"""General mutation: map the two ligands by their maximum common substructure.

`element_swap` matches by atom *name*, which only works when both ligands came
from the same parameterisation run.  `atom_addition` places exactly one new atom,
and the reference has to be the larger ligand.  Neither covers a substitution that
adds a group, removes a group, and renames things at once -- a methyl that becomes
an ethyl, a chlorine that becomes a hydroxyl, an edit to a ring.

This strategy derives the mapping instead of being told it.  RDKit finds the
maximum common substructure (MCS) of the two ligands; the atoms in it are the
hybrid's common core, and everything else vanishes or appears.  `vanish` and
`appear` are then *not* declared in `[mutation]`, and declaring them is an error:
the whole point is that they come from the structures.

Checked against both worked systems (no per-system tuning) it recovers exactly
what the hand-written strategies declare:

    4YLJ   vanish {I1}                 appear {Br1}
    6I5I   vanish {C12, H7, H8, H9}    appear {H -> H17}

Placement is the part that had to be got right, and the obvious answer is wrong.
Aligning the whole mutant onto the reference over the common core -- a single
global Kabsch -- puts 6I5I's appearing hydrogen **5.29 A** from the nitrogen it
is meant to bond to, because the two 6I5I ligands are only a 3.03 A RMSD match
and a least-squares fit over a non-rigid core smears the new atom.  So placement
is **local**, exactly as in `atom_addition`: each appearing group is aligned
through a frame at its own attachment point, so the new atom inherits the
mutant's own bond length and angle.  The same rule gives 0.96 A for 6I5I -- a
correct N-H -- and for 4YLJ, whose two ligands are co-registered to 0.0000 A, it
is the identity transform, which is why this strategy reproduces that frozen
hybrid byte for byte.

What it does NOT do
-------------------
MCS is not a chemistry oracle.  On a molecule with a symmetric core it can pick a
different common substructure than the one intended, and the result will still be
a structurally valid hybrid that looks entirely normal.  The derived mapping is
therefore printed loudly (`provenance`), and the match count is reported when the
MCS maps in more than one way.  If it picks wrong, the remedy is a strategy that
states the answer -- `element_swap` or `atom_addition` with an explicit
declaration -- not a knob here.  This is the general fallback, not a replacement.
"""

from __future__ import annotations

import numpy as np

from rbfe.charmm import angle_deg, element
from rbfe.errors import ChemistryError
from rbfe.model import Ligand, Mapping
from rbfe.mutations.addition import kabsch, next_h_name
from rbfe.mutations.base import MutationSpec, register

# Seconds RDKit may spend on the MCS before giving up.  A timeout is a hard
# failure: a truncated search returns a *smaller* common core, which would
# silently turn a modest mutation into a much larger alchemical one.
DEFAULT_TIMEOUT = 60

# How many symmetry-equivalent correspondences to enumerate per ligand before
# scoring them. 8 is typical; this is a runaway guard, not a tuning knob.
MAX_MATCHES = 64

# A frame triple must bend by at least this much at the anchor. Three collinear
# points do not fix a rotation, and least-squares through them is undefined in
# the direction that matters.
MIN_FRAME_ANGLE = 15.0

# van der Waals radii (A).  Used only to *score* a candidate correspondence,
# never to build a topology -- the force field's own radii come from the .prm.
VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
             "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98}


def worst_overlap(ref: Ligand, mut: Ligand, cand: list[tuple[str, str]]) -> float:
    """Largest fractional vdW overlap between the placed appearing atoms and
    every reference atom, for one candidate correspondence.

    Type agreement cannot separate two correspondences that differ only in
    *which* of several equivalent hydrogens they leave out -- every pair still
    matches on type.  Geometry can.

    3HTB is the case that needs it.  Its core ends at a reference methyl (the
    end of the ethyl) that maps to the mutant's *middle* CH2, so only two of the
    reference's three methyl hydrogens have a partner and one must vanish.  All
    three choices score 18/18 on type, and the one RDKit returns first keeps the
    hydrogen that points straight down the propyl chain -- a retained H **0.43 A**
    inside the appearing methyl carbon.  That builds, passes every existing
    check, and then dies at timestep 1 of the backward leg with a 614394 kcal/mol
    vdW energy.

    The appearing atoms are placed here with a global Kabsch over the common
    heavy atoms.  That is *not* how `_place` positions them -- it uses a local
    frame -- but it is the right tool for ranking, because it needs only to be
    consistently wrong across candidates, and for the co-registered pair it is
    very nearly the identity.
    """
    heavy = [(mut.xyz[m], ref.xyz[r]) for r, m in cand
             if element(mut.atoms[m]) != "H"]
    if len(heavy) < 3:
        return 0.0
    P = np.array([p for p, _ in heavy], float)
    Q = np.array([q for _, q in heavy], float)
    R, t = kabsch(P, Q)

    common = {m for _, m in cand}
    appearing = [n for n in mut.names if n not in common]
    # Only the atoms that are *retained* are checked against.  The vanishing
    # ones are switched off at lambda=1, exactly like the appearing ones, so a
    # vanishing atom sitting near the new group is not a clash -- and counting
    # it would score every candidate the same and defeat the whole exercise.
    fixed = [(r, np.array(ref.xyz[r], float)) for r, _ in cand]

    worst = 0.0
    for n in appearing:
        p = R @ np.array(mut.xyz[n], float) + t
        rn = VDW_RADII.get(element(mut.atoms[n]), 1.70)
        for m, q in fixed:
            rs = rn + VDW_RADII.get(element(ref.atoms[m]), 1.70)
            d = float(np.linalg.norm(p - q))
            if d < rs:
                worst = max(worst, (rs - d) / rs)
    return worst


@register
class MolecularMCS:
    name = "mcs"
    # vanish/appear are derived, so they are deliberately not declared here.
    # `check_declared` is not called.
    required_keys = ()
    optional_keys = ("new_atom_name", "mcs_timeout")
    key_types = {"new_atom_name": "str", "mcs_timeout": "int"}

    def __init__(self, spec: MutationSpec) -> None:
        self.spec = spec

    @classmethod
    def from_spec(cls, spec: MutationSpec) -> "MolecularMCS":
        declared = [n for n in (*spec.vanish, *spec.appear)]
        if declared:
            raise ChemistryError(
                f"[mutation] strategy 'mcs' derives the mapping from the two "
                f"structures, so 'vanish'/'appear' must not be declared (got "
                f"{declared}).\n"
                f"  Remove them, or use a strategy that takes them explicitly: "
                f"element_swap, or atom_addition.")
        if spec.placement not in ("", "native", "local_frame"):
            raise ChemistryError(
                f"[mutation] placement = '{spec.placement}' is not valid for mcs, "
                f"which places every appearing atom by a local frame alignment.")
        return cls(spec)

    # -- the MCS itself ----------------------------------------------------

    def _pairs(self, ref: Ligand, mut: Ligand) -> tuple[list[tuple[str, str]], dict]:
        """[(ref atom name, mutant atom name)] for the common core, plus a report.

        The MCS fixes *which* atoms are common but not *which maps to which*: a
        symmetric core admits several equally large correspondences, and taking
        the first is not safe.  On 6I5I the first correspondence swaps the
        pyrazole's two nitrogens, which builds a structurally valid hybrid of the
        wrong mutation -- it moves the methyl rather than removing it.  So every
        correspondence of the maximal core is scored by how many atom pairs share
        a CHARMM type, and the best-scoring one wins.  On 6I5I that recovers the
        same pairing `atom_addition` derives by hand, at 34 of 40 atoms agreeing.
        """
        from rbfe import similarity as S

        try:
            ma = S.build_mol(ref, keep_h=True)
            mb = S.build_mol(mut, keep_h=True)
        except S.Unavailable as exc:
            raise ChemistryError(
                f"[mutation] strategy 'mcs' needs RDKit: {exc}") from None

        from rdkit import Chem
        from rdkit.Chem import rdFMCS

        res = rdFMCS.FindMCS(
            [ma, mb],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            timeout=self.spec.opt("mcs_timeout", DEFAULT_TIMEOUT),
            completeRingsOnly=True,
        )
        if res.canceled:
            raise ChemistryError(
                f"the maximum common substructure search timed out after "
                f"{self.spec.opt('mcs_timeout', DEFAULT_TIMEOUT)}s.\n"
                f"  A truncated search would return a SMALLER common core and "
                f"silently enlarge the alchemical transformation. Raise "
                f"[mutation] mcs_timeout, or use a strategy that declares the "
                f"mapping explicitly.")
        if res.numAtoms == 0:
            raise ChemistryError(
                f"{ref.stem} and {mut.stem} share no common substructure at all. "
                f"There is no hybrid to build and no meaningful relative free "
                f"energy here; ABFE is the method for a change this large.")

        query = Chem.MolFromSmarts(res.smartsString)
        # uniquify=False: the symmetry-equivalent correspondences are exactly
        # what has to be discriminated between, so they must not be collapsed.
        # maxMatches bounds the search on a highly symmetric molecule.
        hits_a = ma.GetSubstructMatches(query, uniquify=False,
                                        maxMatches=MAX_MATCHES)
        hits_b = mb.GetSubstructMatches(query, uniquify=False,
                                        maxMatches=MAX_MATCHES)
        if not hits_a or not hits_b:
            raise ChemistryError(
                f"the common substructure RDKit found could not be located back "
                f"in {ref.stem} or {mut.stem}; refusing to guess a mapping.")

        name = lambda m, i: m.GetAtomWithIdx(i).GetProp("charmm_name")  # noqa: E731
        best: list[tuple[str, str]] = []
        best_score, best_overlap, considered, ties = -1, 0.0, 0, 0
        for a in hits_a:
            for b in hits_b:
                considered += 1
                cand = [(name(ma, i), name(mb, j)) for i, j in zip(a, b)]
                score = sum(1 for r, m in cand
                            if ref.atoms[r].type == mut.atoms[m].type)
                if score < best_score:
                    continue
                # Type agreement decides; geometry only separates an exact tie.
                # A correspondence that buries a retained atom inside an
                # appearing one is never the one that was meant, however well
                # its types line up.
                ov = worst_overlap(ref, mut, cand)
                if score > best_score:
                    ties = 1
                    best, best_score, best_overlap = cand, score, ov
                else:
                    ties += 1
                    if ov < best_overlap:
                        best, best_overlap = cand, ov

        truncated = (len(hits_a) >= MAX_MATCHES or len(hits_b) >= MAX_MATCHES)
        return best, {
            "core": res.numAtoms,
            "candidates": considered,
            "type_agreement": best_score,
            "tied_on_type": ties,
            "worst_overlap": best_overlap,
            "truncated": truncated,
        }

    # -- placement ---------------------------------------------------------

    @staticmethod
    def _anchor(adj, start: str, common_mut: set[str]) -> str:
        """The common atom the appearing group containing `start` hangs off.

        Walks outward through appearing atoms, so a whole appearing group (a
        methyl's carbon *and* its three hydrogens) resolves to the one atom the
        group is attached to.
        """
        seen, frontier = set(), [start]
        while frontier:
            nxt = []
            for node in frontier:
                if node in seen:
                    continue
                seen.add(node)
                # sorted(): `adj` values are sets, and set iteration order for
                # strings varies with PYTHONHASHSEED. Without this the anchor --
                # and therefore the frame and the placed coordinate -- would
                # differ between runs of the same input.
                for nb in sorted(adj[node]):
                    if nb in common_mut:
                        return nb
                    if nb not in seen:
                        nxt.append(nb)
            frontier = nxt
        raise ChemistryError(
            f"appearing atom '{start}' is not bonded to the common core, so "
            f"there is nothing to place it against. The MCS may have chosen a "
            f"disconnected core; use a strategy that declares the mapping.")

    def _frame(self, mut: Ligand, adj_mut, anchor_m: str, common_mut: set[str],
               group: list[str]) -> list[str]:
        """Two more common heavy atoms that, with the anchor, fix a rigid frame.

        `atom_addition` takes the anchor's two heavy neighbours, which works
        because its anchor is a ring atom with two.  A **chain** anchor has only
        one once the appearing group is excluded -- growing 2-ethylphenol into
        2-propylphenol anchors on the chain's CH2, whose only remaining heavy
        neighbour is the CH2 before it.  So the search walks outward through the
        common core until it has a non-collinear triple.  For the two worked
        systems the anchor already has two neighbours, the first level is taken,
        and the frame is unchanged.
        """
        levels, seen, frontier = [], {anchor_m} | set(group), [anchor_m]
        while frontier and sum(len(l) for l in levels) < 8:
            level, nxt = [], []
            for n in frontier:
                for nb in adj_mut[n]:
                    if nb in seen or nb not in common_mut:
                        continue
                    if element(mut.atoms[nb]) == "H":
                        continue
                    seen.add(nb)
                    level.append(nb)
                    nxt.append(nb)
            if level:
                level.sort(key=lambda x: (mut.atoms[x].type, x))
                levels.append(level)
            frontier = nxt
        order = [a for lvl in levels for a in lvl]

        for i in range(len(order)):
            for j in range(i + 1, len(order)):
                ang = angle_deg(mut.xyz[anchor_m], mut.xyz[order[i]], mut.xyz[order[j]])
                if MIN_FRAME_ANGLE < ang < 180.0 - MIN_FRAME_ANGLE:
                    return [order[i], order[j]]
        raise ChemistryError(
            f"cannot define a local frame at anchor '{anchor_m}': no two of its "
            f"{len(order)} nearby heavy neighbour(s) {order} form a non-collinear "
            f"triple with it.\n"
            f"  This is the symmetric-centre case, where the frame is genuinely "
            f"ambiguous. Use a strategy that declares the mapping explicitly.")

    def _place(self, ref: Ligand, mut: Ligand, adj_mut, group: list[str],
               anchor_m: str, common_mut: set[str],
               ref_of: dict[str, str]) -> np.ndarray:
        """Aligned coordinates for one appearing group, in the reference frame."""
        anchor_r = ref_of[anchor_m]
        frame_m = self._frame(mut, adj_mut, anchor_m, common_mut, group)
        frame_r = [ref_of[n] for n in frame_m]
        rt = [ref.atoms[n].type for n in frame_r]
        mt = [mut.atoms[n].type for n in frame_m]
        if rt != mt:
            raise ChemistryError(
                f"the frames at anchor '{anchor_m}' do not correspond: reference "
                f"types {rt} vs mutant types {mt}. Use a strategy that declares "
                f"the mapping.")

        P = np.array([mut.xyz[anchor_m], *(mut.xyz[n] for n in frame_m)], float)
        Q = np.array([ref.xyz[anchor_r], *(ref.xyz[n] for n in frame_r)], float)
        R, t = kabsch(P, Q)
        return np.array([R @ np.array(mut.xyz[n], float) + t for n in group]), frame_r

    def _planarity_improper(self, ref: Ligand, adj_mut, anchor_m: str,
                            anchor_r: str, group: list[str], vanish: list[str],
                            hybrid_name: dict[str, str]):
        """The planarity improper for a new atom on an sp2 anchor, or None.

        `topology.implied` enumerates impropers at three-coordinate centres, but
        the centre here has an *alchemical* neighbour, so the term it needs is
        the one the generic pass cannot see.  Only a genuinely three-coordinate
        anchor wants one: 6I5I's nitrogen has two ring neighbours plus the new
        H.  A four-coordinate anchor -- 3HTB's chain CH2, with two carbons and
        two hydrogens -- is tetrahedral and must NOT get an improper, and the
        two frame atoms the search may have walked out to are not necessarily
        bonded to the anchor, so they cannot be used to build one.
        """
        substituents = [n for n in ref.adjacency()[anchor_r]
                        if n not in set(vanish)]
        bonded = [n for n in group if anchor_m in adj_mut[n]]
        if len(substituents) + len(bonded) != 3:
            return None
        heavy = sorted((n for n in substituents if element(ref.atoms[n]) != "H"),
                       key=lambda n: (ref.atoms[n].type, n))
        if len(heavy) != 2 or len(bonded) != 1:
            return None
        return (heavy[0], hybrid_name[bonded[0]], anchor_r, heavy[1])

    # -- the mapping -------------------------------------------------------

    def map(self, ref: Ligand, mut: Ligand) -> Mapping:
        pairs, info = self._pairs(ref, mut)

        ref_of = {m: r for r, m in pairs}         # mutant name -> reference name
        common = [n for n in ref.names if n in {r for r, _ in pairs}]
        common_mut = set(ref_of)
        vanish = [n for n in ref.names if n not in set(common)]
        appear_mut = [n for n in mut.names if n not in common_mut]

        if not common:
            raise ChemistryError(
                f"the common substructure of {ref.stem} and {mut.stem} is empty.")
        if not appear_mut:
            raise ChemistryError(
                f"nothing appears in this mutation -- {mut.stem} has no atom "
                f"outside the common core, so there is no atom to absorb the "
                f"charge difference between the two ligands. A pure deletion is "
                f"not supported: swap reference and mutant.")

        # Names for the appearing atoms.  A mutant name is reusable only if the
        # hybrid does not already use it; hydrogens continue the H<n> sequence,
        # which is where 6I5I's H17 comes from.
        used = set(ref.names)
        force = self.spec.opt("new_atom_name")
        if force and len(appear_mut) != 1:
            raise ChemistryError(
                f"[mutation] new_atom_name names one atom, but {len(appear_mut)} "
                f"appear: {appear_mut}")
        hybrid_name: dict[str, str] = {}
        for n in appear_mut:
            if force:
                new = force
            elif element(mut.atoms[n]) == "H":
                new = next_h_name(used)
            elif n not in used:
                new = n
            else:
                i = 2
                while f"{n}_{i}" in used:
                    i += 1
                new = f"{n}_{i}"
            if new in used:
                raise ChemistryError(
                    f"the name '{new}' for appearing atom '{n}' is already in the "
                    f"hybrid. Set [mutation] new_atom_name.")
            used.add(new)
            hybrid_name[n] = new

        # Local frame per attachment point.  Groups are placed independently so
        # a multi-site substitution does not drag one group to fit another.
        adj_mut = mut.adjacency()
        groups: dict[str, list[str]] = {}
        for n in appear_mut:
            anchor = self._anchor(adj_mut, n, common_mut)
            groups.setdefault(anchor, []).append(n)

        coords = {n: ref.xyz[n] for n in ref.names}
        provenance = [
            f"core is the maximum common substructure: {len(pairs)} atoms "
            f"matched by element and bond order",
            f"chosen from {info['candidates']} correspondence(s) of that core "
            f"by CHARMM type agreement ({info['type_agreement']} of "
            f"{info['core']} atom pairs share a type)",
            f"derived vanish {sorted(vanish)}, appear "
            f"{[hybrid_name[n] for n in appear_mut]}",
        ]
        if info["truncated"]:
            provenance.append(
                f"note: the correspondence search hit its {MAX_MATCHES}-match cap, "
                f"so a better-scoring mapping may exist. If this mapping looks "
                f"wrong, declare the mutation with element_swap or atom_addition.")

        bad_impropers: set[tuple[str, str, str, str]] = set()
        for anchor_m, group in groups.items():
            placed, frame_r = self._place(ref, mut, adj_mut, group, anchor_m,
                                          common_mut, ref_of)
            anchor_r = ref_of[anchor_m]
            for n, xyz in zip(group, placed):
                coords[hybrid_name[n]] = tuple(float(v) for v in xyz)
                provenance.append(
                    f"'{hybrid_name[n]}' placed at ({xyz[0]:.3f} {xyz[1]:.3f} "
                    f"{xyz[2]:.3f}), {np.linalg.norm(xyz - np.array(ref.xyz[anchor_r])):.3f} A "
                    f"from anchor '{anchor_r}' (local frame {frame_r})")
            term = self._planarity_improper(ref, adj_mut, anchor_m, anchor_r,
                                            group, vanish, hybrid_name)
            if term:
                bad_impropers.add(term)

        # The two ligands name their atoms consistently only if every common pair
        # shares a name AND every appearing atom kept its mutant name.  When they
        # do, the mutant's own bond/improper lists describe the hybrid and merge
        # cleanly; when they do not, merging them invents connectivity, and every
        # bond touching an appearing atom has to be supplied by hand.
        share_names = (all(r == m for r, m in pairs)
                       and all(hybrid_name[n] == n for n in appear_mut))
        hname = {**ref_of, **hybrid_name}

        extra_bonds: set[tuple[str, str]] = set()
        if not share_names:
            for u, v in mut.bonds:
                if u in appear_mut or v in appear_mut:
                    extra_bonds.add(tuple(sorted((hname[u], hname[v]))))  # type: ignore[arg-type]

        types = {n: ref.atoms[n].type for n in ref.names}
        charges = {n: ref.atoms[n].charge for n in ref.names}
        for n in appear_mut:
            types[hybrid_name[n]] = mut.atoms[n].type
            charges[hybrid_name[n]] = mut.atoms[n].charge   # native, before closure

        return Mapping(
            common=common,
            vanish=vanish,
            appear=[hybrid_name[n] for n in appear_mut],
            types=types,
            charges=charges,
            coords=coords,
            ref_order=list(ref.names),
            share_names=share_names,
            provenance=provenance,
            extra_bonds=extra_bonds,
            # Only when the mutant's terms are discarded -- otherwise its own
            # improper for the appearing atom is merged in already.
            extra_impropers=set() if share_names else bad_impropers,
        )
