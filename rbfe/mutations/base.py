"""The mutation-strategy interface.

The problem this solves
-----------------------
The two systems in this repo need structurally different hybrid-builders: 6I5I
(N-CH3 -> N-H) *adds* an atom the reference does not have, so that builder hunts
the methyl and the N-H by local chemical environment and invents the new
position with a rigid alignment; 4YLJ (I -> Br) adds and removes nothing, so its
builder just matches the two ligands by atom name.  Historically that meant two
copies of `prepare_hybrid.py` that diverged by 632 lines.

Here the *chemistry* lives in a strategy and the *plumbing* -- parsing,
connectivity, charge closure, parameter sourcing, emission -- lives once in
`rbfe.hybrid`.  A strategy is small: it says which atoms vanish, which appear,
what they are called, and where the appearing ones go.

Failing loudly
--------------
A strategy declares `REQUIRED_KEYS` and `OPTIONAL_KEYS`.  The config loader
rejects any `[mutation]` key not claimed by the named strategy, and any strategy
name not in the registry.  So a mutation the engine does not understand cannot
be silently approximated -- which is the whole reason this interface exists.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar, Protocol, runtime_checkable

from rbfe.errors import ChemistryError
from rbfe.model import Ligand, Mapping

# Keys every strategy may read.  `strategy` selects the strategy; the rest
# describe the mutation itself.
BASE_KEYS: dict[str, str] = {
    "strategy": "str",
    "vanish": "list",
    "appear": "list",
    "placement": "str",
    "missing_terms": "str",
}

MISSING_TERM_POLICIES = ("fail", "mirror", "mirror_then_synth", "synth_from_geometry")


@dataclass(frozen=True)
class MutationSpec:
    """The `[mutation]` section, already coerced, handed to a strategy.

    Strategies read `options` by key; anything they do not declare is rejected
    earlier, by the config loader, so a strategy can trust every key it reads.
    """

    strategy: str
    vanish: tuple[str, ...]
    appear: tuple[str, ...]
    placement: str
    missing_terms: str = "mirror_then_synth"
    options: dict[str, Any] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.options is None:
            object.__setattr__(self, "options", {})

    def opt(self, key: str, default: Any = None) -> Any:
        return self.options.get(key, default)

    def need(self, key: str) -> Any:
        if key not in self.options:
            raise ChemistryError(
                f"[mutation] strategy '{self.strategy}' requires '{key}'")
        return self.options[key]


def check_declared(spec: MutationSpec, ref: Ligand, mut: Ligand) -> None:
    """The declared vanish/appear names must exist, once each, in the right ligand.

    Run by every strategy before it maps anything, so the error names the
    offending atom rather than surfacing later as a mysterious charge-closure
    failure or an unmarked atom in the built system.
    """
    for label, names, lig in (("vanish", spec.vanish, ref),
                              ("appear", spec.appear, mut)):
        unknown = [n for n in names if n not in lig.atoms]
        if unknown:
            raise ChemistryError(
                f"[mutation] {label} names atom(s) not in "
                f"{Path(lig.stem).name}: {', '.join(unknown)}\n"
                f"  available: {', '.join(lig.names)}")
    both = set(spec.vanish) & set(spec.appear)
    if both:
        raise ChemistryError(
            f"[mutation] {', '.join(sorted(both))} is listed as both vanishing "
            f"and appearing. In a dual topology those are opposite roles; the "
            f"same name cannot be both.")
    dup = [n for n in spec.vanish if spec.vanish.count(n) > 1]
    dup += [n for n in spec.appear if spec.appear.count(n) > 1]
    if dup:
        raise ChemistryError(f"[mutation] duplicate name(s): {', '.join(sorted(set(dup)))}")


@runtime_checkable
class MutationStrategy(Protocol):
    """What `rbfe.hybrid` calls.  Implementations live beside this file."""

    name: ClassVar[str]
    required_keys: ClassVar[tuple[str, ...]]
    optional_keys: ClassVar[tuple[str, ...]]

    @classmethod
    def from_spec(cls, spec: MutationSpec) -> "MutationStrategy": ...

    def map(self, ref: Ligand, mut: Ligand) -> Mapping: ...


REGISTRY: dict[str, type] = {}


def register(cls: type) -> type:
    """Register a strategy under its `name`.  Called at import time."""
    name = getattr(cls, "name", None)
    if not name:
        raise AssertionError(f"{cls.__name__} has no `name`")
    if name in REGISTRY:
        raise AssertionError(f"duplicate mutation strategy {name!r}")
    REGISTRY[name] = cls
    return cls


def mutation_key_schema(strategy: str | None = None) -> dict[str, str]:
    """Legal `[mutation]` keys -> their config types.

    With no argument, the union across every registered strategy -- used only to
    bound the universe.  With a strategy named, exactly that strategy's keys,
    which is what makes an unrecognised key fatal rather than ignored.
    """
    if strategy is None:
        keys = dict(BASE_KEYS)
        for cls in REGISTRY.values():
            for k in cls.required_keys + cls.optional_keys:
                keys.setdefault(k, "list")
        return keys

    try:
        cls = REGISTRY[strategy]
    except KeyError:
        known = ", ".join(sorted(REGISTRY)) or "(none implemented)"
        raise ChemistryError(
            f"unknown [mutation] strategy '{strategy}'.\n"
            f"  implemented: {known}\n"
            f"  A mutation the engine does not understand is not approximated -- "
            f"add a strategy in rbfe/mutations/ rather than editing the builder."
        ) from None

    keys = dict(BASE_KEYS)
    for k in cls.required_keys + cls.optional_keys:
        keys[k] = getattr(cls, "key_types", {}).get(k, "list")
    return keys
