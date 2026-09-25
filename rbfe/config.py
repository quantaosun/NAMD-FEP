"""Per-system configuration: one `system.ini`, validated against a schema.

Why INI and not TOML/YAML/JSON
------------------------------
`tomllib` is Python 3.11+ and this repo supports 3.8 (it is absent on the box
this was written on, 3.10.10).  YAML would mean adding a dependency to a project
whose documented contract is "Python 3.8+ with numpy" -- and the whole point of
this repo is that it runs from a bare Python with no stack.  JSON is stdlib but
has no comments, which is fatal for a file whose entire purpose is to be
hand-edited and explained.  `configparser` is stdlib, takes comments, has
sections, and has been in the language since forever.

The one thing this module exists to guarantee
---------------------------------------------
**An unknown key is fatal.**  `vansih = C12 H7 H8 H9` must not quietly produce a
hybrid with no vanishing set: that hybrid would have no -1 atoms, the run would
produce a one-directional alchFile, and the resulting ddG would be meaningless
while looking perfectly normal.  So the schema below is a closed world: every
legal (section, key) pair is listed, and anything else stops the job.  For the
same reason unknown *sections* are fatal, and `[DEFAULT]` is refused outright
(configparser leaks DEFAULT keys into every section, which would defeat the
whole check).

The `[mutation]` section is the one dynamic part: which keys are legal depends on
which strategy is named.  Those keys are supplied by the strategy registry in
`rbfe.mutations`, so a strategy that does not understand a key cannot silently
ignore it.
"""

from __future__ import annotations

import configparser
import os
import re
from pathlib import Path
from typing import Any

from rbfe.errors import ConfigError
from rbfe.mutations import mutation_key_schema

CONFIG_NAME = "system.ini"

# section -> key -> type
#
# Types:
#   str    -- verbatim, whitespace-trimmed
#   int/float/bool
#   list   -- whitespace-separated (newlines from continuation lines included)
#   lines  -- one item per line, internal spaces preserved
#   floats -- whitespace-separated, parsed as float
SCHEMA: dict[str, dict[str, str]] = {
    "system": {
        "name": "str",
        "description": "str",
    },
    "ligand": {
        "dir": "str",          # where ref/mut live, relative to the system root
        "ref": "str",          # stem: <dir>/<ref>.{rtf,prm,pdb}
        "mut": "str",
        "resname": "str",      # residue name the hybrid is written under
    },
    "protein": {
        "pdb": "str",
        "chain": "str",        # "" means every chain
        "segment_prefix": "str",
        "first": "str",        # psfgen terminus, e.g. NTER
        "last": "str",         # e.g. CTER
    },
    "mutation": {},            # filled in from the strategy registry
    "hybrid": {
        "title": "str",
        "type_prefix": "str",  # "raw" or a literal prefix such as "L"
        # The .rtf and the .prm are ordered independently -- 4YLJ sorts its
        # topology terms but emits its parameters in source order.  One key
        # cannot say that, and forcing one to serve both silently reorders a
        # file that was correct.
        "rtf_order": "str",    # "sorted" | "source_then_new"
        "prm_order": "str",    # "sorted" | "source_then_new"
    },
    "similarity": {
        # Screens the ligand pair for RBFE suitability.  A property of the pair,
        # not of a mutation strategy, so it lives in its own section rather than
        # in the strategy-driven [mutation] schema.
        "threshold": "float",  # ECFP4 Tanimoto below which the pair is queried
        "enabled": "bool",
    },
    "build": {
        "padding": "float",    # solvate -t
        "neutralize": "bool",
        "salt": "float",       # autoionize -c ; 0.0 omits the flag
        "pdbalias": "lines",
    },
    "topology": {
        "prot_rtf": "str",     # path to top_all36_prot.rtf, relative to system root
        "params": "lines",     # emitted VERBATIM into every .namd (see inputs.py)
    },
    "run": {
        "legs": "list",
        "lambdas": "floats",
        "steps_per_window": "int",
        "equil_steps": "int",
        "alch_equil_steps": "int",
        "outfreq": "int",
        "min_steps": "int",
        "timestep": "float",
        "temperature": "float",
        "cutoff": "float",
        "switchdist": "float",
        "pairlistdist": "float",
        "elec_lambda_start": "float",
        "vdw_lambda_end": "float",
        "vdw_shift_coeff": "float",
        "decouple": "str",
        "max_hours": "float",
    },
    "binaries": {
        "namd": "str",
        "vmd": "str",
        "flags": "list",       # e.g. +p1 +devices 0  (separate tokens)
        "gpu_resident": "str",
    },
}

# Keys with no sensible universal default.  Everything else has one below.
REQUIRED: set[tuple[str, str]] = {
    ("system", "name"),
    ("ligand", "ref"),
    ("ligand", "mut"),
    ("ligand", "resname"),
    ("protein", "pdb"),
    ("protein", "chain"),
    ("mutation", "strategy"),
    ("topology", "prot_rtf"),
    ("topology", "params"),
    ("run", "lambdas"),
    ("run", "steps_per_window"),
    ("binaries", "namd"),
}

DEFAULTS: dict[tuple[str, str], str] = {
    ("system", "description"): "",
    ("ligand", "dir"): "inputs",
    ("protein", "segment_prefix"): "P",
    ("protein", "first"): "NTER",
    ("protein", "last"): "CTER",
    ("mutation", "vanish"): "",
    ("mutation", "appear"): "",
    ("mutation", "placement"): "native",
    ("mutation", "missing_terms"): "fail",
    ("hybrid", "title"): "dual-topology hybrid ligand",
    ("hybrid", "type_prefix"): "L",
    ("hybrid", "rtf_order"): "sorted",
    ("hybrid", "prm_order"): "source_then_new",
    ("similarity", "threshold"): "0.60",
    ("similarity", "enabled"): "yes",
    ("build", "padding"): "15.0",
    ("build", "neutralize"): "yes",
    ("build", "salt"): "0.0",
    ("build", "pdbalias"): "atom ILE CD1 CD\natom SER HG HG1\nresidue HIS HSD",
    ("run", "legs"): "complex solvent",
    ("run", "equil_steps"): "50000",
    ("run", "alch_equil_steps"): "50000",
    ("run", "outfreq"): "500",
    ("run", "min_steps"): "5000",
    ("run", "timestep"): "2.0",
    ("run", "temperature"): "300.0",
    ("run", "cutoff"): "12.0",
    ("run", "switchdist"): "10.0",
    ("run", "pairlistdist"): "14.0",
    ("run", "elec_lambda_start"): "0.5",
    ("run", "vdw_lambda_end"): "1.0",
    ("run", "vdw_shift_coeff"): "5.0",
    ("run", "decouple"): "off",
    ("run", "max_hours"): "0.0",      # 0 disables the wall-clock check
    ("binaries", "vmd"): "vmd",
    ("binaries", "flags"): "+p1 +devices 0",
    ("binaries", "gpu_resident"): "on",
}


_ENV_DEFAULT = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*):-([^}]*)\}")


def _expand_env(s: str) -> str:
    """`$VAR`, `${VAR}` and `${VAR:-default}`.

    `os.path.expandvars` handles the first two but not the `:-` default, which
    is the form a config file wants: the machine-specific binary paths get a
    working default while staying overridable by environment, matching the
    `os.environ.get("NAMD", <path>)` pattern the shell scripts already use.
    """
    s = _ENV_DEFAULT.sub(lambda m: os.environ.get(m.group(1), m.group(2)), s)
    return os.path.expandvars(s)


def _coerce(value: str, kind: str, where: str) -> Any:
    v = value.strip()
    try:
        if kind == "str":
            return v
        if kind == "int":
            return int(v)
        if kind == "float":
            return float(v)
        if kind == "floats":
            return [float(x) for x in v.split()]
        if kind == "list":
            return v.split()
        if kind == "lines":
            return [ln.strip() for ln in v.splitlines() if ln.strip()]
        if kind == "bool":
            low = v.lower()
            if low in ("yes", "true", "on", "1"):
                return True
            if low in ("no", "false", "off", "0"):
                return False
            raise ConfigError(f"{where}: expected a yes/no value, got {value!r}")
    except ValueError as exc:
        raise ConfigError(f"{where}: cannot read {value!r} as {kind} ({exc})") from None
    raise AssertionError(f"unknown schema type {kind!r}")


class Config:
    """A validated system.ini, plus the directory it lives in."""

    def __init__(self, root: Path, values: dict[tuple[str, str], Any],
                 raw: dict[tuple[str, str], str] | None = None) -> None:
        self.root = root
        self._v = values
        self._raw = raw or {}

    def __repr__(self) -> str:
        return f"<Config {self.get('system', 'name')} at {self.root}>"

    def has(self, section: str, key: str) -> bool:
        return (section, key) in self._v

    def get(self, section: str, key: str) -> Any:
        if (section, key) in self._v:
            return self._v[(section, key)]
        raise ConfigError(f"[{section}] {key} is not set and has no default")

    # -- derived values ----------------------------------------------------
    @property
    def name(self) -> str:
        return self.get("system", "name")

    @property
    def legs(self) -> list[str]:
        return self.get("run", "legs")

    @property
    def lambdas(self) -> list[float]:
        return self.get("run", "lambdas")

    @property
    def lambda_tokens(self) -> list[str]:
        """The lambda values exactly as the config wrote them.

        `fep.tcl` is emitted from these rather than from the parsed floats, so
        the generated Tcl is byte-identical to what was written for the frozen
        systems -- `repr()` of a float is normally the shortest round-trip form,
        but that is a property of the interpreter and not something to bet a
        reproduction check on.
        """
        raw = self._raw.get(("run", "lambdas"))
        return raw.split() if raw else [repr(x) for x in self.lambdas]

    @property
    def nwin(self) -> int:
        """Number of lambda windows = n_lambdas - 1."""
        return len(self.lambdas) - 1

    @property
    def steps_per_window(self) -> int:
        return self.get("run", "steps_per_window")

    @property
    def equil_rows(self) -> int:
        """Rows per window that are equilibration, not production.

        NAMD writes `alchEquilSteps` steps of equilibration at `alchOutFreq`,
        and the first row lands one interval in, so the count is one less than
        the naive division.  For 50000/500 that is 99 -- the number that used to
        be copy-pasted into three separate files.
        """
        return self.get("run", "alch_equil_steps") // self.get("run", "outfreq") - 1

    @property
    def expect_rows(self) -> int:
        """Rows a complete window should contain."""
        return self.steps_per_window // self.get("run", "outfreq")

    def resolve(self, rel: str) -> Path:
        """A path from the config, resolved against the system root.

        `{root}` expands to the system directory, so a config can name a
        sibling tree (the shared `toppar/`) without depending on how deep the
        system directory happens to sit -- which is exactly the bug that made
        6I5I's `../../toppar/...` strings wrong under `systems/<name>/`.
        `{hybrid}` is left alone here: it depends on `--out`, which the config
        does not know.
        """
        s = _expand_env(rel.replace("{root}", str(self.root)))
        p = Path(s)
        return p if p.is_absolute() else (self.root / p)

    def leg_dir(self, leg: str) -> Path:
        return self.root / leg

    def mutation_section(self) -> dict[str, Any]:
        """The [mutation] section as a plain dict, for the strategy to consume."""
        return {k: v for (s, k), v in self._v.items() if s == "mutation"}


def _mutation_schema_from(parser: configparser.ConfigParser, path: Path) -> dict[str, str]:
    """The legal `[mutation]` keys, given the strategy this file names.

    Resolved before any other key checking, so a key that belongs to a
    *different* strategy is rejected too -- not just an outright typo.
    """
    if not parser.has_section("mutation"):
        raise ConfigError(
            f"{path}: no [mutation] section. Every system must declare what the "
            f"perturbation is; there is no default.")
    raw = parser["mutation"].get("strategy")
    if raw is None or not raw.strip():
        raise ConfigError(
            f"{path}: [mutation] strategy is required. "
            f"(An empty strategy would mean the engine guessed at the mutation, "
            f"which is exactly what this pipeline refuses to do.)")
    return mutation_key_schema(raw.strip())


def _validate_keys(parser: configparser.ConfigParser, path: Path) -> dict[str, str]:
    """Reject unknown sections and keys before anything reads a value."""
    if parser.defaults():
        raise ConfigError(
            f"{path}: [DEFAULT] is not allowed. configparser copies DEFAULT keys "
            f"into every section, which would defeat the unknown-key check that "
            f"stops a typo'd key from silently disabling part of the mutation.")

    mut_schema = _mutation_schema_from(parser, path)
    schema = {s: dict(k) for s, k in SCHEMA.items()}
    schema["mutation"] = mut_schema

    known = set(schema)
    for section in parser.sections():
        if section not in known:
            raise ConfigError(
                f"{path}: unknown section [{section}]. "
                f"Known sections: {', '.join(sorted(known))}.")
        allowed = set(schema[section])
        for key in parser[section]:
            if key not in allowed:
                near = _closest(key, allowed)
                hint = f"  did you mean '{near}'?" if near else ""
                extra = ""
                if section == "mutation":
                    extra = (f"\n  (keys are checked against strategy "
                             f"'{parser['mutation'].get('strategy', '').strip()}'; "
                             f"a key from a different strategy is not accepted)")
                raise ConfigError(
                    f"{path}: unknown key '{key}' in [{section}].{hint}{extra}\n"
                    f"  allowed here: {', '.join(sorted(allowed)) or '(none)'}")
    return mut_schema


def _closest(word: str, options: set[str]) -> str | None:
    """Cheap typo suggestion -- enough for a transposition or a missing letter."""
    def dist(a: str, b: str) -> int:
        prev = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            cur = [i]
            for j, cb in enumerate(b, 1):
                cur.append(min(prev[j] + 1, cur[j - 1] + 1,
                               prev[j - 1] + (ca != cb)))
            prev = cur
        return prev[-1]

    best = min(options, key=lambda o: dist(word, o), default=None)
    if best is None:
        return None
    return best if dist(word, best) <= max(2, len(word) // 3) else None


def load(path: Path) -> Config:
    """Read and fully validate one system.ini."""
    path = Path(path).resolve()
    if not path.exists():
        raise ConfigError(f"no such config: {path}")

    parser = configparser.ConfigParser(
        interpolation=None,   # configs mention % in prose comments
        strict=True,          # a duplicate key is a mistake, not a merge
        delimiters=("=",),    # never split on ':'
    )
    try:
        with path.open() as fh:
            parser.read_file(fh, source=str(path))
    except configparser.Error as exc:
        raise ConfigError(f"{path}: {exc}") from None

    mut_schema = _validate_keys(parser, path)
    root = path.parent

    def kind_of(section: str, key: str) -> str:
        return SCHEMA.get(section, {}).get(key) or mut_schema[key]

    values: dict[tuple[str, str], Any] = {}
    raw_values: dict[tuple[str, str], str] = {}
    for section in parser.sections():
        for key, raw in parser[section].items():
            values[(section, key)] = _coerce(
                raw, kind_of(section, key), f"{path}: [{section}] {key}")
            raw_values[(section, key)] = raw

    for (section, key), raw in DEFAULTS.items():
        if (section, key) in values:
            continue
        if section == "mutation" and key not in mut_schema:
            continue          # this strategy does not declare that key
        values[(section, key)] = _coerce(
            raw, kind_of(section, key), f"default for [{section}] {key}")

    missing = [f"[{s}] {k}" for (s, k) in sorted(REQUIRED) if (s, k) not in values]
    if missing:
        raise ConfigError(f"{path}: required keys absent: {', '.join(missing)}")

    cfg = Config(root, values, raw_values)
    _validate_semantics(cfg, path)
    return cfg


def _validate_semantics(cfg: Config, path: Path) -> None:
    """Checks that need more than one value at a time."""
    lam = cfg.lambdas
    if len(lam) < 2:
        raise ConfigError(f"{path}: [run] lambdas needs at least 2 values, got {len(lam)}")
    if lam != sorted(lam):
        raise ConfigError(f"{path}: [run] lambdas must be ascending, got {lam}")
    if abs(lam[0]) > 0.0 or abs(lam[-1] - 1.0) > 0.0:
        raise ConfigError(
            f"{path}: [run] lambdas must run 0.0 -> 1.0 (got {lam[0]} -> {lam[-1]}); "
            f"an incomplete alchemical path gives a dG that is missing its endpoints.")
    if len(set(lam)) != len(lam):
        raise ConfigError(f"{path}: [run] lambdas has duplicates: {lam}")
    if cfg.get("run", "outfreq") <= 0:
        raise ConfigError(f"{path}: [run] outfreq must be positive")
    if cfg.steps_per_window % cfg.get("run", "outfreq"):
        raise ConfigError(
            f"{path}: [run] steps_per_window ({cfg.steps_per_window}) is not a multiple "
            f"of outfreq ({cfg.get('run', 'outfreq')}); the window would end mid-interval "
            f"and the driver's per-window row count would be wrong.")
    if cfg.equil_rows < 0:
        raise ConfigError(
            f"{path}: [run] alch_equil_steps ({cfg.get('run', 'alch_equil_steps')}) is "
            f"smaller than outfreq ({cfg.get('run', 'outfreq')}), so no equilibration row "
            f"is ever written.")
    if not cfg.legs:
        raise ConfigError(f"{path}: [run] legs is empty")
    if cfg.get("hybrid", "type_prefix") == "":
        raise ConfigError(
            f"{path}: [hybrid] type_prefix is empty. Use 'raw' to reproduce a "
            f"legacy build explicitly; an empty prefix is almost certainly a mistake.")

    # the input files the pipeline actually opens
    lig = cfg.resolve(cfg.get("ligand", "dir"))
    for stem_key in ("ref", "mut"):
        stem = lig / cfg.get("ligand", stem_key)
        for ext in (".rtf", ".prm", ".pdb"):
            if not stem.with_suffix(ext).exists():
                raise ConfigError(
                    f"{path}: [ligand] {stem_key} = {cfg.get('ligand', stem_key)!r} "
                    f"-> missing {stem.with_suffix(ext)}")
    prot = cfg.resolve(cfg.get("protein", "pdb"))
    if not prot.exists():
        raise ConfigError(f"{path}: [protein] pdb -> missing {prot}")
    prot_rtf = cfg.resolve(cfg.get("topology", "prot_rtf"))
    if not prot_rtf.exists():
        raise ConfigError(f"{path}: [topology] prot_rtf -> missing {prot_rtf}")


def find(start: Path | None = None) -> Config:
    """Load the system.ini governing `start` (default: cwd), walking upward.

    Commands therefore work from anywhere inside a system directory, so the
    step-by-step commands never need a path argument.
    """
    here = Path(start or Path.cwd()).resolve()
    for d in [here, *here.parents]:
        candidate = d / CONFIG_NAME
        if candidate.exists():
            return load(candidate)
    raise ConfigError(
        f"no {CONFIG_NAME} found in {here} or any parent directory.\n"
        f"  Run this inside a system directory, or create one with "
        f"`rbfe init systems/<name>`.")
