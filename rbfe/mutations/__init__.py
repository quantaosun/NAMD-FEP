"""The mutation-strategy registry.

Importing this package registers every strategy.  `rbfe.config` imports
`mutation_key_schema` from here to validate the `[mutation]` section, so this
module must not import `rbfe.config` (it does not).
"""

from __future__ import annotations

from rbfe.mutations.base import (
    BASE_KEYS,
    MISSING_TERM_POLICIES,
    REGISTRY,
    MutationSpec,
    MutationStrategy,
    check_declared,
    mutation_key_schema,
    register,
)
from rbfe.mutations.addition import AtomAddition
from rbfe.mutations.element_swap import ElementSwap
from rbfe.mutations.mcs import MolecularMCS

__all__ = [
    "BASE_KEYS", "MISSING_TERM_POLICIES", "REGISTRY", "MutationSpec",
    "MutationStrategy", "check_declared", "mutation_key_schema", "register",
    "AtomAddition", "ElementSwap", "MolecularMCS", "get_strategy",
    "spec_from_section",
]

# The implemented strategies are `element_swap`, `atom_addition` and `mcs`. A
# system naming anything else gets the registry's loud error rather than a
# silent fallback, which is the intended failure mode for an unrecognised
# mutation.
#
# `element_swap` and `atom_addition` are the specialists: each states its
# mapping, and each reproduces one frozen system byte for byte. `mcs` derives
# the mapping from the structures and is the general path. Note that importing
# this module must not import RDKit -- `rbfe.config` imports it for the key
# schema, so a top-level RDKit import here would make every subcommand fail on a
# machine that does not have it. `mcs.py` imports RDKit inside its methods.


def get_strategy(name: str) -> type:
    """The strategy class for `name`, or a loud failure listing what exists."""
    if name not in REGISTRY:
        # Reuse the schema call purely for its error message.
        mutation_key_schema(name)
    return REGISTRY[name]


def spec_from_section(section: dict) -> MutationSpec:
    """Build a MutationSpec from a validated `[mutation]` section."""
    strategy = section["strategy"]
    declared = mutation_key_schema(strategy)          # validates the name
    known = set(declared)
    options = {k: v for k, v in section.items()
               if k not in ("strategy", "vanish", "appear") and k in known}
    return MutationSpec(
        strategy=strategy,
        vanish=tuple(section.get("vanish", [])),
        appear=tuple(section.get("appear", [])),
        placement=section.get("placement", "native"),
        missing_terms=section.get("missing_terms", "fail"),
        options=options,
    )
