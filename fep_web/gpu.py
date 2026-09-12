"""GPU detection and a launch preflight.

The NAMD build used here is `--with-single-node-cuda` (GPU-resident), so a job
started without a usable GPU does not merely run slowly -- it fails after
consuming the setup time. This module answers one question before a launch:
*is there a GPU, and does it have room for us right now?*

Both `smi` (the executable) and the command runner are injectable so the tests
can simulate "no nvidia-smi", "no devices", and "not enough free memory"
without a GPU present.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field

# NAMD peaks at ~590 MiB on the 6I5I complex (see CLAUDE.md). The margin is
# generous on purpose: this card is SHARED, and a launch that fits at t=0 can
# be squeezed out by another process a few minutes in.
DEFAULT_REQUIRED_MIB = 1500.0

_QUERY = "index,name,memory.total,memory.used,compute_cap"


@dataclass(frozen=True)
class GpuDevice:
    index: int
    name: str
    total_mib: float
    used_mib: float
    compute_cap: str

    @property
    def free_mib(self) -> float:
        return self.total_mib - self.used_mib


@dataclass
class PreflightResult:
    ok: bool
    reason: str
    devices: list[GpuDevice] = field(default_factory=list)
    required_mib: float = DEFAULT_REQUIRED_MIB

    @property
    def best_free_mib(self) -> float:
        return max((d.free_mib for d in self.devices), default=0.0)


class GpuUnavailable(RuntimeError):
    """Raised by require_gpu() when no usable device is present."""


def query_devices(smi: str = "nvidia-smi", timeout: float = 10.0) -> list[GpuDevice]:
    """Return the visible CUDA devices. Raises RuntimeError with a human reason."""
    exe = shutil.which(smi) if not smi.startswith("/") else smi
    if not exe or not shutil.os.path.exists(exe):
        raise RuntimeError(
            f"'{smi}' not found. There is no NVIDIA driver/SMI on this host, so "
            f"GPU-resident NAMD cannot run here.")

    try:
        proc = subprocess.run(
            [exe, f"--query-gpu={_QUERY}", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"'{smi}' timed out after {timeout}s — the driver "
                           f"may be wedged.")
    except OSError as exc:
        raise RuntimeError(f"could not run '{smi}': {exc}")

    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip().splitlines()
        raise RuntimeError(f"'{smi}' exited {proc.returncode}: "
                           f"{detail[-1] if detail else 'no output'}")

    devices: list[GpuDevice] = []
    for line in proc.stdout.splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 5 or not parts[0]:
            continue
        try:
            devices.append(GpuDevice(
                index=int(parts[0]), name=parts[1],
                total_mib=float(parts[2]), used_mib=float(parts[3]),
                compute_cap=parts[4]))
        except ValueError:
            continue                      # an unparseable row is not fatal
    if not devices:
        raise RuntimeError(f"'{smi}' reported no CUDA devices.")
    return devices


def preflight(required_mib: float = DEFAULT_REQUIRED_MIB,
              smi: str = "nvidia-smi") -> PreflightResult:
    """Decide whether it is safe to launch. Never raises for GPU problems."""
    try:
        devices = query_devices(smi)
    except RuntimeError as exc:
        return PreflightResult(False, str(exc), [], required_mib)

    best = max(devices, key=lambda d: d.free_mib)
    if best.free_mib < required_mib:
        return PreflightResult(
            False,
            f"GPU '{best.name}' has only {best.free_mib:.0f} MiB free of "
            f"{best.total_mib:.0f} MiB, but {required_mib:.0f} MiB is needed. "
            f"The card is shared — wait for the other job or pick another device.",
            devices, required_mib)

    return PreflightResult(
        True,
        f"{best.name} ({best.compute_cap}) — {best.free_mib:.0f} MiB free "
        f"of {best.total_mib:.0f} MiB.",
        devices, required_mib)


def require_gpu(required_mib: float = DEFAULT_REQUIRED_MIB,
                smi: str = "nvidia-smi") -> PreflightResult:
    """Like preflight(), but raise GpuUnavailable instead of returning ok=False."""
    result = preflight(required_mib, smi)
    if not result.ok:
        raise GpuUnavailable(result.reason)
    return result
