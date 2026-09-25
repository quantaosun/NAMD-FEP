"""A general, config-driven relative binding free energy workflow for NAMD.

One system is one directory containing a `system.ini`; every step of the
pipeline is the same command for every system.  See RUNBOOK.md.

The engine is deliberately dependency-free beyond numpy: charmm parsing,
topology merging, charge closure and NAMD config rendering are all stdlib.
"""

from __future__ import annotations

__version__ = "0.1.0"
