"""
fep_web — host-agnostic job control for the NAMD RBFE workflow.

This package deliberately knows nothing about Flask, HTTP, or where it is
running. It wraps the proven on-disk conventions (`fep_run.py` per-window
drivers, `.done` markers, `run_checkpointed.sh` / `respawn_controller.sh`) behind
a small interface:

    status(system)              -> structured state of every leg and stage
    submit(system)              -> start or resume the controller (detached)
    cancel(system)              -> stop the controller and any namd3
    logs(system, name, tail)    -> recent lines of a leg's log
    preflight(...)              -> is there a usable GPU right now?

A "system" is a directory laid out like `6I5I_DUAL_FEP/`: it contains the leg
sub-directories (`complex/`, `solvent/`), `fep_run.py`, and the controller
scripts. Everything the UI needs comes from those files, so the UI can be
restarted, moved to another host, or replaced without touching the job state.

The GPU check matters here because the box is a single **shared** V100: another
process may already hold memory, and starting NAMD on a contended card either
fails or crawls. `preflight()` therefore checks *free* memory, not mere presence.
"""
from .gpu import GpuDevice, PreflightResult, preflight, query_devices
from .jobs import JobError, cancel, logs, status, submit
from .state import SystemState, read_state

__version__ = "0.1.0"

__all__ = [
    "GpuDevice", "PreflightResult", "preflight", "query_devices",
    "JobError", "cancel", "logs", "status", "submit",
    "SystemState", "read_state",
]
