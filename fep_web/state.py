"""Read job state off the filesystem.

The single source of truth is the set of `.done` markers that `fep_run.py`
writes -- one per completed lambda window, only on a clean exit. That is
deliberately *not* the same as "the .fepout file exists": a killed window leaves
a partial .fepout behind, and trusting it is what would let a truncated window be
analysed as if it were complete.

Everything here is read-only and cheap enough to call on every UI request.
"""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

WINDOWS = 15                     # lambdas = 16, windows = 15
LEGS = ("complex", "solvent")
ROWS_PER_WINDOW = 250000 // 500  # alchOutFreq 500

# stage -> how its completion is recorded
SINGLE_STAGES = ("nvt_equil", "npt_equil")
WINDOW_STAGES = ("md_forward", "md_backward")

_FEP_ROW = re.compile(r"^FepEnergy:")


@dataclass
class StageState:
    name: str
    done: int
    total: int
    running: bool = False
    rows_in_progress: int = 0
    note: str = ""          # e.g. "finished without a marker" — shown in the UI

    @property
    def complete(self) -> bool:
        return self.done >= self.total

    @property
    def fraction(self) -> float:
        if self.complete:
            return 1.0
        if self.running and self.rows_in_progress:
            # sub-window progress, so the UI moves during a 500 ps window
            part = min(self.rows_in_progress / ROWS_PER_WINDOW, 1.0)
            return (self.done + part) / self.total
        return self.done / self.total if self.total else 0.0


@dataclass
class LegState:
    name: str
    stages: dict[str, StageState] = field(default_factory=dict)

    @property
    def complete(self) -> bool:
        return all(s.complete for s in self.stages.values())


@dataclass
class SystemState:
    path: Path
    legs: dict[str, LegState] = field(default_factory=dict)
    controller_pid: int | None = None
    namd_pids: list[int] = field(default_factory=list)
    active_window: str | None = None
    error: str | None = None

    @property
    def running(self) -> bool:
        return self.controller_pid is not None or bool(self.namd_pids)

    @property
    def complete(self) -> bool:
        return all(leg.complete for leg in self.legs.values()) and bool(self.legs)


def _namd_processes(runner=subprocess.run) -> list[tuple[int, str]]:
    """[(pid, cmdline)] for running namd3 processes, or [] if none/unavailable."""
    try:
        proc = runner(["pgrep", "-x", "namd3"], capture_output=True, text=True,
                      timeout=5)
    except (OSError, subprocess.TimeoutExpired):
        return []
    if proc.returncode != 0:
        return []
    out = []
    for token in proc.stdout.split():
        try:
            pid = int(token)
        except ValueError:
            continue
        try:
            cmd = Path(f"/proc/{pid}/cmdline").read_bytes().replace(b"\0", b" ").decode(
                errors="replace")
        except OSError:
            cmd = ""
        out.append((pid, cmd))
    return out


def _controller_pid(system: Path, runner=subprocess.run) -> int | None:
    """PID of a live controller/wrapper for this system, if any.

    `respawn_controller.sh` writes `controller_run.pid`; it is treated as a hint
    only, and confirmed against /proc so a stale file does not report a job as
    running forever.
    """
    for pidfile in (system / "controller_run.pid", system / "run_all.pid"):
        try:
            pid = int(pidfile.read_text().split()[0])
        except (OSError, ValueError, IndexError):
            continue
        if Path(f"/proc/{pid}").exists():
            return pid
    return None


def _count_rows(path: Path) -> int:
    try:
        with open(path, errors="replace") as fh:
            return sum(1 for ln in fh if _FEP_ROW.match(ln))
    except OSError:
        return 0


def single_run_windows(path: Path) -> set[int]:
    """Indices of lambda windows completed by a *single-run* (pre-controller) fepout.

    NAMD writes one '#NEW FEP WINDOW' line per process, so a multi-window file
    has to be split by step range: window i spans (i*250000, (i+1)*250000], and
    it is complete only if the row at exactly (i+1)*250000 is present. This is
    the same rule `fep_run.py single_window_data` uses.

    It matters because 6I5I's complex/forward windows 0-5 were salvaged from a
    crashed single run and therefore have no `.done` markers at all -- counting
    markers alone would report a finished leg as 9/15.
    """
    have: set[int] = set()
    try:
        with open(path, errors="replace") as fh:
            for ln in fh:
                if not _FEP_ROW.match(ln):
                    continue
                parts = ln.split()
                if len(parts) < 2 or not parts[1].isdigit():
                    continue
                step = int(parts[1])
                if step % 250000 == 0:
                    have.add(step // 250000 - 1)
    except OSError:
        return set()
    return {i for i in have if 0 <= i < WINDOWS}


def _stage_finished_without_marker(legdir: Path, stage: str) -> bool:
    """True if a stage clearly finished before markers existed.

    The original `run_all.sh` wrote no markers, so the complex leg's nvt/npt
    equilibration looks 'not started' to a marker-only check. Require BOTH the
    restart file (which downstream stages read as `bincoordinates`) and a log
    that reached NAMD's clean end -- output alone is not proof of completion.
    """
    if not (legdir / f"{stage}.coor").exists():
        return False
    log = legdir / f"{stage}.log"
    try:
        tail = log.read_text(errors="replace")[-4000:]
    except OSError:
        return False
    return "End of program" in tail


def read_state(system: Path, runner=subprocess.run) -> SystemState:
    """Snapshot of every leg and stage. Never raises for missing files."""
    system = Path(system)
    st = SystemState(path=system)

    if not system.is_dir():
        st.error = f"{system} is not a directory"
        return st

    namd = _namd_processes(runner)
    st.namd_pids = [pid for pid, _ in namd]
    st.controller_pid = _controller_pid(system, runner)

    # which window is the running namd3 working on?
    for _pid, cmd in namd:
        m = re.search(r"(md_(?:forward|backward)|nvt_equil|npt_equil)_w?\d*\.namd", cmd)
        if m:
            st.active_window = m.group(0).removesuffix(".namd")
            break

    for leg in LEGS:
        legdir = system / leg
        if not legdir.is_dir():
            continue
        ls = LegState(name=leg)
        for stage in SINGLE_STAGES:
            marked = (legdir / f".{stage}.done").exists()
            done = marked or _stage_finished_without_marker(legdir, stage)
            ls.stages[stage] = StageState(
                name=stage, done=int(done), total=1,
                running=st.active_window == stage,
                note="" if marked or not done else "finished without a marker")
        for stage in WINDOW_STAGES:
            marked = {i for i in range(WINDOWS)
                      if (legdir / f".{stage}_w{i:02d}.done").exists()}
            # a crashed pre-controller single run can still have completed
            # windows, which carry no markers
            salvaged = single_run_windows(legdir / f"{stage}.fepout") - marked
            running = st.active_window is not None and \
                st.active_window.startswith(stage)
            rows = 0
            if running:
                idx = st.active_window.rsplit("_w", 1)[-1]
                if idx.isdigit():
                    rows = _count_rows(legdir / f"{stage}_w{int(idx):02d}.fepout")
            ls.stages[stage] = StageState(
                name=stage, done=len(marked) + len(salvaged), total=WINDOWS,
                running=running, rows_in_progress=rows,
                note=(f"{len(salvaged)} salvaged single-run window(s)"
                      if salvaged else ""))
        st.legs[leg] = ls

    return st
