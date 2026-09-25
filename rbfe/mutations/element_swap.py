"""Element swap: one atom becomes a different element at the same ring position.

The 4YLJ case -- ligand 4E1, 10-iodo -> 10-bromo.  Nothing is added and nothing
is removed: both ligands carry the same atoms under the same names (acpype
guarantees this), exactly one differs, and the mutant's `.pdb` already holds it
at the right place.  So the mapping is a name match and no geometry is invented.

Contrast with `atom_addition`, which has to identify a chemical group and place a
new atom by rigid alignment because the reference has no such atom at all.
"""

from __future__ import annotations

from rbfe.errors import ChemistryError
from rbfe.model import Ligand, Mapping
from rbfe.mutations.base import MutationSpec, check_declared, register


@register
class ElementSwap:
    name = "element_swap"
    required_keys = ("vanish", "appear")
    optional_keys = ("placement", "missing_terms")

    # vanish/appear are atom names, so a whitespace list; the rest are scalars.
    key_types = {"vanish": "list", "appear": "list",
                 "placement": "str", "missing_terms": "str"}

    def __init__(self, spec: MutationSpec) -> None:
        self.spec = spec

    @classmethod
    def from_spec(cls, spec: MutationSpec) -> "ElementSwap":
        if spec.placement not in ("native", ""):
            raise ChemistryError(
                f"[mutation] placement = '{spec.placement}' is not valid for "
                f"element_swap, which always keeps each atom's own input geometry. "
                f"Use 'native' (or omit it).")
        return cls(spec)

    def map(self, ref: Ligand, mut: Ligand) -> Mapping:
        check_declared(self.spec, ref, mut)

        common = [n for n in ref.names if n in mut.atoms]
        ref_only = [n for n in ref.names if n not in mut.atoms]
        mut_only = [n for n in mut.names if n not in ref.atoms]

        if not common:
            raise ChemistryError(
                f"no atoms matched by name between {ref.stem} and {mut.stem}.\n"
                f"  element_swap matches by atom NAME, which requires both ligands "
                f"to come from the same parameterisation run (acpype output "
                f"satisfies this). If they were parameterised independently the "
                f"names will differ throughout -- use a strategy that matches by "
                f"chemical environment instead.")

        if set(ref_only) != set(self.spec.vanish) or set(mut_only) != set(self.spec.appear):
            raise ChemistryError(
                f"the declared mutation does not describe these ligands.\n"
                f"  declared vanish {sorted(self.spec.vanish)}  appear {sorted(self.spec.appear)}\n"
                f"  actual   vanish {sorted(ref_only)}  appear {sorted(mut_only)}\n"
                f"  ([mutation] must account for every atom that differs between "
                f"the two ligands -- an unlisted one would be silently dropped "
                f"from the hybrid.)")

        if not mut_only:
            raise ChemistryError(
                "the two ligands have identical atom names, so nothing appears "
                "or vanishes -- there is no alchemical perturbation to run.")

        provenance = [f"matched {len(common)} atom(s) by name"]
        for n in common:
            rt, mt = ref.atoms[n].type, mut.atoms[n].type
            if rt != mt:
                provenance.append(f"'{n}' changes type {rt} -> {mt}")

        types: dict[str, str] = {}
        charges: dict[str, float] = {}
        coords: dict[str, tuple[float, float, float]] = {}

        for n in common + ref_only:
            types[n] = ref.atoms[n].type
            charges[n] = ref.atoms[n].charge
            coords[n] = ref.xyz[n]

        for n in mut_only:
            if n not in mut.xyz:
                raise ChemistryError(
                    f"appearing atom '{n}' has no coordinates in {mut.stem}.pdb. "
                    f"element_swap keeps each atom's own input geometry and "
                    f"invents nothing; a strategy that places atoms by alignment "
                    f"is the one for a mutation that adds a group.")
            types[n] = mut.atoms[n].type
            charges[n] = mut.atoms[n].charge
            coords[n] = mut.xyz[n]
            provenance.append(
                f"'{n}' appears at its own mutant geometry "
                f"({coords[n][0]:.3f} {coords[n][1]:.3f} {coords[n][2]:.3f})")

        return Mapping(common=common, vanish=ref_only, appear=mut_only,
                       ref_order=list(ref.names),
                       types=types, charges=charges, coords=coords,
                       provenance=provenance)
