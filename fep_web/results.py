"""Read analysis results for a system.

Shells out to the system's own `analyze_fep.py` with `--json`, so nothing here
parses human-readable output. The analyzer lives next to the data it describes
(each system carries its own copy), which is why this is a subprocess rather
than an import.
"""
from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

ANALYZER = "analyze_fep.py"

# The four legs BAR needs, in the order analyze_fep.py expects them.
COMBINED = (
    "complex/md_forward_combined.fepout",
    "complex/md_backward_combined.fepout",
    "solvent/md_forward_combined.fepout",
    "solvent/md_backward_combined.fepout",
)


class AnalysisError(RuntimeError):
    """The analysis could not be run or did not produce usable output."""


@dataclass
class Readiness:
    ready: bool
    missing: list[str] = field(default_factory=list)
    reason: str = ""


def readiness(system: Path) -> Readiness:
    """Are the four combined fepouts present, and is the analyzer there?"""
    system = Path(system)
    analyzer = system / ANALYZER
    if not analyzer.exists():
        return Readiness(False, [], f"{ANALYZER} not found in {system.name}")
    missing = [f for f in COMBINED if not (system / f).exists()]
    if missing:
        return Readiness(False, missing,
                         f"{len(missing)} of {len(COMBINED)} combined fepouts "
                         f"not assembled yet")
    return Readiness(True)


def bar_result(system: Path, runner=subprocess.run, timeout: float = 600.0,
               python: str | None = None) -> dict:
    """Run the two-sided BAR analysis and return its JSON.

    Raises AnalysisError with a message that is safe to show a user.
    """
    system = Path(system)
    r = readiness(system)
    if not r.ready:
        raise AnalysisError(r.reason)

    cmd = [python or sys.executable, ANALYZER, "bar", *COMBINED, "--json"]
    try:
        proc = runner(cmd, cwd=str(system), capture_output=True, text=True,
                      timeout=timeout)
    except subprocess.TimeoutExpired:
        raise AnalysisError(f"analysis timed out after {timeout:.0f}s")
    except OSError as exc:
        raise AnalysisError(f"could not run the analyzer: {exc}")

    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        raise AnalysisError(
            f"analyzer exited {proc.returncode}: "
            f"{tail[-1] if tail else 'no output'}")

    try:
        data = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise AnalysisError(f"analyzer produced no usable JSON ({exc})")

    for key in ("dG_complex", "dG_solvent", "ddG"):
        if key not in data:
            raise AnalysisError(f"analyzer JSON is missing '{key}'")
    return data
