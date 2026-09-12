"""Start / stop / inspect a system's FEP job.

This wraps the scripts that are already proven on this box rather than
reimplementing them:

    respawn_controller.sh  -> relaunches run_checkpointed.sh whenever the box
                              kills it (~4.5 h), until ALL LEGS DONE
      run_checkpointed.sh  -> marker-guarded, resumable, strictly sequential legs
        fep_run.py         -> one NAMD process per lambda window

So `submit()` means "start or resume", and it is safe to call repeatedly. The
one thing added here is the **GPU preflight**, which the shell scripts do not do:
they assume a free card, and on this shared V100 that assumption is often wrong.
"""
from __future__ import annotations

import os
import signal
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .gpu import DEFAULT_REQUIRED_MIB, PreflightResult, preflight
from .state import SystemState, read_state

CONTROLLER = "respawn_controller.sh"
CONTROLLER_LOG = "controller_run.log"


class JobError(RuntimeError):
    """A refusal or a failure to launch. The message is safe to show a user."""


@dataclass
class SubmitResult:
    started: bool
    pid: int | None
    reason: str
    gpu: PreflightResult


def status(system: Path, runner=subprocess.run) -> SystemState:
    """Current state of every leg and stage. Read-only, never raises."""
    return read_state(Path(system), runner=runner)


def logs(system: Path, name: str = CONTROLLER_LOG, tail: int = 200) -> list[str]:
    """Last `tail` lines of a log file inside the system directory."""
    path = Path(system) / name
    if not path.exists():
        return []
    try:
        with open(path, errors="replace") as fh:
            return fh.read().splitlines()[-tail:]
    except OSError:
        return []


def submit(system: Path, *, required_mib: float = DEFAULT_REQUIRED_MIB,
           allow_no_gpu: bool = False, controller: str = CONTROLLER,
           smi: str = "nvidia-smi", runner=subprocess.run,
           popen=subprocess.Popen) -> SubmitResult:
    """Start or resume the job. Refuses on a bad GPU or if one is already live.

    `allow_no_gpu=True` skips the preflight, for exercising the plumbing on a
    machine without a card. It is NOT a way to run GPU-resident NAMD on CPU --
    the configs set `GPUresident on` and NAMD will refuse to start.
    """
    system = Path(system)
    if not system.is_dir():
        raise JobError(f"{system} is not a directory")

    # Check "already running" BEFORE the script-exists check: if a job is live,
    # that is the message the user needs, and reporting "script not found" for
    # a directory that is obviously working would be actively misleading.
    current = status(system, runner=runner)
    if current.running:
        who = (f"controller pid {current.controller_pid}" if current.controller_pid
               else f"namd3 pid {current.namd_pids[0]}")
        raise JobError(f"already running ({who}). Use cancel() first.")

    script = system / controller
    if not script.exists():
        raise JobError(f"{script} not found — is this a prepared system directory?")

    gpu = preflight(required_mib, smi)
    if not gpu.ok and not allow_no_gpu:
        raise JobError(f"GPU preflight failed: {gpu.reason}")

    log = open(system / CONTROLLER_LOG, "a")           # controller appends too
    log.write(f"\n--- launched by fep_web at pid {os.getpid()} ---\n")
    log.flush()
    try:
        proc = popen(["bash", script.name],
                     cwd=str(system), stdout=log, stderr=subprocess.STDOUT,
                     stdin=subprocess.DEVNULL,
                     start_new_session=True)           # survive the UI/session
    except OSError as exc:
        log.close()
        raise JobError(f"could not start {script.name}: {exc}")
    log.close()

    note = "GPU ok" if gpu.ok else "GPU check SKIPPED (allow_no_gpu)"
    return SubmitResult(True, proc.pid, f"started {script.name} (pid {proc.pid}); {note}",
                        gpu)


def cancel(system: Path, *, timeout: float = 10.0, runner=subprocess.run) -> str:
    """Stop the controller and any namd3 it launched. Safe when nothing runs."""
    system = Path(system)
    current = status(system, runner=runner)
    pids = ([current.controller_pid] if current.controller_pid else []) + \
           list(current.namd_pids)
    if not pids:
        return "nothing was running"

    for pid in pids:
        for sig in (signal.SIGTERM,):
            try:
                os.kill(pid, sig)
            except ProcessLookupError:
                pass
            except PermissionError:
                raise JobError(f"not permitted to signal pid {pid}")
    return f"sent SIGTERM to {', '.join(str(p) for p in pids)}"
