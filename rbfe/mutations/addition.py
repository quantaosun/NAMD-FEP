"""Atom addition: the mutant has an atom the reference does not.

The 6I5I case -- a pyrazole N-CH3 that becomes N-H.  A methyl is removed and a
hydrogen appears in its place, so unlike an element swap there is no mutant atom
standing where the new one goes: its position has to be *computed*.

Where an element swap matches by atom name, this strategy cannot.  acpype
renumbers independently parameterised ligands, so 6I5I's ref and mut share
almost no atom names.  Instead:

* the hybrid's core is the **reference's** atoms, minus the vanishing group,
  keeping the reference's names, types and charges;
* the appearing atom is **newly named** (the next free `H<n>`, which is how
  6I5I's `H17` arises) and typed from the mutant;
* its position comes from a rigid alignment of the mutant's local frame onto the
  reference's, so the N-H sits at the mutant's own bond length and angle but in
  the reference's frame of coordinates.

The frame is derived, not configured.  The historical builder hardcoded the
types `cc`/`nc` as "the two ring neighbours"; here they are the anchor's
non-vanishing, non-hydrogen neighbours, ordered by type name -- which for 6I5I
yields exactly `cc` then `nc`, so the alignment is reproduced bit-for-bit.  If
that set is not exactly two atoms with matching types, the strategy refuses and
asks for `frame_atoms` rather than guessing: an alignment built on the wrong
three points produces a plausible hybrid with a wrong N-H position, which is
far worse than stopping.

`anchor`, `frame_atoms` and `new_atom_name` are the explicit overrides for the
cases where the derivation is legitimately ambiguous (a symmetric centre).
"""

from __future__ import annotations

import re

import numpy as np

from rbfe.charmm import element
from rbfe.errors import ChemistryError
from rbfe.model import Ligand, Mapping
from rbfe.mutations.base import MutationSpec, check_declared, register


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Rotation + translation taking P onto Q (least squares, proper rotation)."""
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def next_h_name(names) -> str:
    """The next unused `H<n>`, which is where `H17` comes from."""
    used = {int(m.group(1)) for n in names
            if (m := re.fullmatch(r"H(\d+)", n.strip(), re.I))}
    return f"H{max(used) + 1}" if used else "H1"


@register
class AtomAddition:
    name = "atom_addition"
    required_keys = ("vanish", "appear")
    optional_keys = ("placement", "missing_terms", "frame_atoms",
                     "new_atom_name", "anchor", "improper")

    key_types = {"vanish": "list", "appear": "list", "placement": "str",
                 "missing_terms": "str", "frame_atoms": "list",
                 "new_atom_name": "str", "anchor": "str", "improper": "list"}

    def __init__(self, spec: MutationSpec) -> None:
        self.spec = spec

    @classmethod
    def from_spec(cls, spec: MutationSpec) -> "AtomAddition":
        if spec.placement not in ("local_frame", ""):
            raise ChemistryError(
                f"[mutation] placement = '{spec.placement}' is not valid for "
                f"atom_addition; the new atom is placed by a local-frame "
                f"alignment. Use 'local_frame' (or omit it).")
        if len(spec.appear) != 1:
            raise ChemistryError(
                f"atom_addition places exactly one new atom per run; "
                f"{len(spec.appear)} were declared ({list(spec.appear)}).")
        return cls(spec)

    # -- derivation ----------------------------------------------------
    @staticmethod
    def _heavy(lig: Ligand, adj, name: str, exclude: set[str]) -> list[str]:
        return sorted((n for n in adj[name]
                       if n not in exclude and element(lig.atoms[n]) != "H"),
                      key=lambda n: (lig.atoms[n].type, n))

    def _anchor(self, ref: Ligand, mut: Ligand) -> tuple[str, str]:
        """(reference anchor, mutant anchor) -- the atom the group hangs off."""
        if self.spec.opt("anchor"):
            return self.spec.need("anchor"), self.spec.need("anchor")
        ra, ma = ref.adjacency(), mut.adjacency()
        vanish = set(self.spec.vanish)
        appear = self.spec.appear[0]

        # The anchor is the one atom *outside* the vanishing group that the group
        # hangs off -- not an atom bonded to every vanishing atom, which would be
        # wrong for a group like a methyl: its hydrogens bond only to its carbon.
        # (The historical builder expressed this as "the N among the methyl
        # carbon's neighbours".)
        outside = {n for v in vanish for n in ra[v]} - vanish
        if len(outside) != 1:
            raise ChemistryError(
                f"cannot identify the anchor: the vanishing group "
                f"{sorted(vanish)} hangs off {len(outside)} atom(s) "
                f"{sorted(outside)} rather than one. Set [mutation] anchor "
                f"explicitly to the atom in the reference the group is attached to.")
        cand = sorted(outside)
        mcand = [n for n in mut.names if appear in ma[n]]
        if len(mcand) != 1:
            raise ChemistryError(
                f"cannot identify the mutant anchor: {len(mcand)} atom(s) "
                f"{mcand} are bonded to the appearing atom '{appear}'. "
                f"Set [mutation] anchor explicitly.")
        return cand[0], mcand[0]

    def _frames(self, ref: Ligand, mut: Ligand, anchor_r: str, anchor_m: str):
        """The two atoms defining the local frame, in matching order."""
        forced = self.spec.opt("frame_atoms")
        if forced:
            if len(forced) != 2:
                raise ChemistryError(
                    f"[mutation] frame_atoms needs exactly 2 reference atom "
                    f"names, got {len(forced)}")
            rf = list(forced)
            ra = ref.adjacency()
            mf = [n for n in mut.adjacency()[anchor_m]
                  if mut.atoms[n].type in {ref.atoms[x].type for x in rf}
                  and element(mut.atoms[n]) != "H"]
            if len(mf) != 2:
                raise ChemistryError(
                    f"[mutation] frame_atoms {rf} could not be matched to two "
                    f"mutant atoms by type (found {mf})")
            return rf, mf

        rf = self._heavy(ref, ref.adjacency(), anchor_r, set(self.spec.vanish))
        mf = self._heavy(mut, mut.adjacency(), anchor_m, set(self.spec.appear))
        if len(rf) != 2 or len(mf) != 2:
            raise ChemistryError(
                f"cannot derive the local frame: the reference anchor "
                f"'{anchor_r}' has {len(rf)} non-vanishing heavy neighbour(s) "
                f"{rf} and the mutant anchor '{anchor_m}' has {len(mf)} {mf}; "
                f"two each are needed to define a rigid frame.\n"
                f"  Set [mutation] frame_atoms to the two reference atom names "
                f"explicitly -- guessing here would place the new atom wrongly "
                f"and look entirely normal doing it.")
        rt = [ref.atoms[n].type for n in rf]
        mt = [mut.atoms[n].type for n in mf]
        if rt != mt:
            raise ChemistryError(
                f"the derived frames do not correspond: reference types {rt} "
                f"vs mutant types {mt}. Set [mutation] frame_atoms explicitly.")
        return rf, mf

    # -- mapping -------------------------------------------------------
    def map(self, ref: Ligand, mut: Ligand) -> Mapping:
        check_declared(self.spec, ref, mut)
        vanish = list(self.spec.vanish)
        appear_name = self.spec.appear[0]
        common = [n for n in ref.names if n not in set(vanish)]

        anchor_r, anchor_m = self._anchor(ref, mut)
        frame_r, frame_m = self._frames(ref, mut, anchor_r, anchor_m)

        P = np.array([mut.xyz[anchor_m], *(mut.xyz[n] for n in frame_m)])
        Q = np.array([ref.xyz[anchor_r], *(ref.xyz[n] for n in frame_r)])
        R, t = kabsch(P, Q)
        new_xyz = tuple(R @ np.array(mut.xyz[appear_name]) + t)
        d = float(np.linalg.norm(np.array(new_xyz) - np.array(ref.xyz[anchor_r])))

        new_name = self.spec.opt("new_atom_name") or next_h_name(ref.names)
        if new_name in ref.atoms:
            raise ChemistryError(
                f"the new atom name '{new_name}' already exists in the "
                f"reference; set [mutation] new_atom_name.")

        # Every reference atom keeps its own type, charge and coordinates --
        # including the vanishing group, which stays in the hybrid at B-factor
        # -1 and must therefore be in these tables even though it is not in
        # `common`.
        types = {n: ref.atoms[n].type for n in ref.names}
        charges = {n: ref.atoms[n].charge for n in ref.names}
        coords = {n: ref.xyz[n] for n in ref.names}
        types[new_name] = mut.atoms[appear_name].type
        charges[new_name] = mut.atoms[appear_name].charge     # native, before closure
        coords[new_name] = new_xyz

        provenance = [
            f"core is the reference's {len(common)} atoms (ref and mut share few "
            f"names, so matching is by structure, not by name)",
            f"anchor '{anchor_r}' (ref) <-> '{anchor_m}' (mut)",
            f"frame {frame_r} (ref types {[ref.atoms[n].type for n in frame_r]}) "
            f"<-> {frame_m} (mut types {[mut.atoms[n].type for n in frame_m]})",
            f"placed '{new_name}' at ({new_xyz[0]:.3f}, {new_xyz[1]:.3f}, "
            f"{new_xyz[2]:.3f}), {d:.3f} A from the anchor",
        ]
        for n, a in zip(frame_r, frame_m):
            if a != n:
                provenance.append(f"frame '{n}' (ref) <-> '{a}' (mut)")

        improper = self._improper(ref, frame_r, new_name, anchor_r)
        return Mapping(common=common, vanish=vanish, appear=[new_name],
                       types=types, charges=charges, coords=coords,
                       ref_order=list(ref.names), share_names=False,
                       extra_bonds={tuple(sorted((anchor_r, new_name)))},
                       provenance=provenance,
                       extra_impropers={improper} if improper else set())

    def _improper(self, ref: Ligand, frame_r: list[str], new_name: str,
                  anchor_r: str):
        """The N-H planarity improper, mirrored from the mutant.

        The merged bond list cannot imply this term: H17 is terminal and bonded
        to a 3-coordinate nitrogen, and the generic enumeration only synthesises
        impropers at 3-coordinate centres *involving* an alchemical atom -- this
        centre's alchemical neighbour is H17 itself, so it would be missed.
        6I5I's builder writes `C H N N1` using the two ring neighbours; that is
        the same term.
        """
        forced = self.spec.opt("improper")
        if forced:
            if len(forced) != 4:
                raise ChemistryError("[mutation] improper needs exactly 4 atom names")
            return tuple(forced)
        if len(frame_r) == 2:
            return (frame_r[0], new_name, anchor_r, frame_r[1])
        return None
