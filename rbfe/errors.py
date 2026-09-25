"""Failures that should stop the pipeline with a message, not a traceback.

Every class here is a SystemExit, so `raise ConfigError(...)` prints the message
and exits non-zero.  The point of this module is the house rule that a mistake in
the inputs must never be silently absorbed: a typo'd config key, an atom name
that is not in the ligand, or a mutation that does not add up must all stop the
job *before* it burns GPU time on a plausible-looking wrong answer.
"""

from __future__ import annotations


class RbfeError(SystemExit):
    """Base: print the message to stderr and exit 1."""

    def __init__(self, msg: str) -> None:
        super().__init__(f"rbfe: {msg}")


class ConfigError(RbfeError):
    """system.ini is malformed, incomplete, or names something that does not exist."""


class ChemistryError(RbfeError):
    """The ligands or the mutation declaration do not describe a valid hybrid."""


class BuildError(RbfeError):
    """psfgen / solvate / autoionize failed, or produced the wrong thing."""


class VerificationError(RbfeError):
    """A reproduction check did not match what it was supposed to reproduce."""
