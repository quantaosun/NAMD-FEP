"""Plain data types shared by the parsers, the strategies and the builder.

Deliberately dependency-free (no other `rbfe` imports) so that anything may
import it -- in particular the mutation strategies, which must not import the
config layer or we would get a cycle.
"""

from __future__ import annotations

from dataclasses import dataclass, field

XYZ = tuple[float, float, float]


@dataclass(frozen=True)
class Atom:
    name: str
    type: str
    charge: float
    mass: float


@dataclass
class Ligand:
    """One input ligand: topology, parameters and coordinates, as parsed.

    `atoms` preserves file order, which matters because the hybrid must list its
    atoms in the same order as the reference for the emitted `.rtf` to be
    byte-comparable with the one this pipeline has historically produced.
    """

    stem: str
    atoms: dict[str, Atom] = field(default_factory=dict)
    bonds: tuple[tuple[str, str], ...] = ()
    angles: tuple[tuple[str, str, str], ...] = ()
    dihedrals: tuple[tuple[str, str, str, str], ...] = ()
    impropers: tuple[tuple[str, str, str, str], ...] = ()
    xyz: dict[str, XYZ] = field(default_factory=dict)

    @property
    def total_charge(self) -> float:
        return sum(a.charge for a in self.atoms.values())

    @property
    def names(self) -> list[str]:
        return list(self.atoms)

    def adjacency(self) -> dict[str, set[str]]:
        adj: dict[str, set[str]] = {n: set() for n in self.atoms}
        for a, b in self.bonds:
            if a in adj and b in adj:
                adj[a].add(b)
                adj[b].add(a)
        return adj


@dataclass
class Mapping:
    """How the reference and mutant ligands become one hybrid.

    Produced by a `MutationStrategy`, consumed by the builder.  Everything the
    builder needs to emit the hybrid is here, so the builder itself contains no
    chemistry.
    """

    common: list[str]
    vanish: list[str]          # ref-only  -> B-factor -1
    appear: list[str]          # mut-only  -> B-factor +1
    types: dict[str, str]      # hybrid atom name -> source CHARMM type
    charges: dict[str, float]  # AFTER charge closure
    coords: dict[str, XYZ]
    # The order the reference ligand lists its atoms in.  The hybrid must keep
    # it -- the vanishing atoms belong *where the reference had them*, not
    # appended at the end, or the emitted .rtf and .pdb list every atom after
    # the perturbation site in a different order from the reference.
    ref_order: list[str] = field(default_factory=list)
    # Whether the two ligands name their atoms consistently.  element_swap: yes,
    # acpype gave both the same names, so the mutant's bond list describes the
    # same molecule and must be merged.  atom_addition: no -- the ligands were
    # parameterised independently and share almost no names, so a "bond" C-H1
    # read out of the mutant is a coincidence of labels, not a real bond, and
    # merging it invents connectivity the hybrid does not have.
    share_names: bool = True
    provenance: list[str] = field(default_factory=list)
    # Impropers the strategy wants added beyond the generic enumeration, as
    # ordered 4-tuples of hybrid names.  Strategies that need none leave it empty.
    extra_impropers: set[tuple[str, str, str, str]] = field(default_factory=set)
    # Bonds the strategy invents.  A merged bond list cannot imply them: for
    # atom_addition the new atom's bond exists in the *mutant* under a different
    # name, and the mutant's bonds are not merged precisely because its names do
    # not correspond.  Without this the hybrid has an atom with no bond to
    # anything, and psfgen builds a disconnected fragment.
    extra_bonds: set[tuple[str, str]] = field(default_factory=set)

    @property
    def order(self) -> list[str]:
        """Hybrid atom order: the reference's order, then the appearing atoms.

        The appearing atoms go last because the reference has no position for
        them; everything the reference did have keeps its original slot.
        """
        base = self.ref_order or (self.common + self.vanish)
        return base + self.appear

    @property
    def alchemical(self) -> set[str]:
        return set(self.vanish) | set(self.appear)


@dataclass
class HybridTopology:
    """Connectivity implied by the merged bond list."""

    bonds: set[tuple[str, str]] = field(default_factory=set)
    angles: set[tuple[str, str, str]] = field(default_factory=set)
    dihedrals: set[tuple[str, str, str, str]] = field(default_factory=set)
    impropers: set[tuple[str, str, str, str]] = field(default_factory=set)
    adjacency: dict[str, set[str]] = field(default_factory=dict)

    @property
    def counts(self) -> tuple[int, int, int, int]:
        return (len(self.bonds), len(self.angles),
                len(self.dihedrals), len(self.impropers))
