"""Tests for fep_web — no GPU, no NAMD, no network.

Run:  python3 -m unittest discover -s tests -v
      (or)  python3 tests/test_fep_web.py

The GPU tests use a *real* fake `nvidia-smi` executable so the subprocess path
is exercised for real; the process tests use an injected runner/Popen.
"""
from __future__ import annotations

import os
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fep_web import gpu, jobs, state          # noqa: E402


def _fake_smi(tmp: Path, script_body: str, name: str = "nvidia-smi") -> str:
    """Write an executable stand-in for nvidia-smi and return its path.

    `name` must differ per fake: two fakes sharing a filename means the second
    silently overwrites the first and both paths report the same GPU.
    """
    p = tmp / name
    p.write_text("#!/bin/sh\n" + script_body + "\n")
    p.chmod(p.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return str(p)


class _CP:
    """Minimal CompletedProcess stand-in for injected runners."""

    def __init__(self, rc=0, out="", err=""):
        self.returncode, self.stdout, self.stderr = rc, out, err


class _IdleRunner:
    """A runner that finds no namd3 processes."""

    def __call__(self, *a, **kw):
        return _CP(rc=1, out="")


# --------------------------------------------------------------------------
class TestGpuQuery(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_missing_nvidia_smi_is_a_clear_error(self):
        with self.assertRaises(RuntimeError) as cm:
            gpu.query_devices(smi=str(self.tmp / "definitely-not-here"))
        self.assertIn("not found", str(cm.exception))

    def test_parses_devices(self):
        smi = _fake_smi(self.tmp, 'echo "0, Tesla V100-SXM2-32GB, 16384, 4317, 7.0"')
        devs = gpu.query_devices(smi=smi)
        self.assertEqual(len(devs), 1)
        self.assertEqual(devs[0].name, "Tesla V100-SXM2-32GB")
        self.assertEqual(devs[0].total_mib, 16384)
        self.assertEqual(devs[0].used_mib, 4317)
        self.assertAlmostEqual(devs[0].free_mib, 12067)

    def test_nonzero_exit_is_reported(self):
        smi = _fake_smi(self.tmp, 'echo "driver mismatch" >&2; exit 9')
        with self.assertRaises(RuntimeError) as cm:
            gpu.query_devices(smi=smi)
        self.assertIn("exited 9", str(cm.exception))
        self.assertIn("driver mismatch", str(cm.exception))

    def test_no_devices_reported(self):
        smi = _fake_smi(self.tmp, "exit 0")
        with self.assertRaises(RuntimeError) as cm:
            gpu.query_devices(smi=smi)
        self.assertIn("no CUDA devices", str(cm.exception))

    def test_garbage_rows_are_skipped_not_fatal(self):
        smi = _fake_smi(self.tmp, 'echo "not,a,row"; echo "0, V100, 16384, 100, 7.0"')
        devs = gpu.query_devices(smi=smi)
        self.assertEqual(len(devs), 1)


class TestPreflight(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_ok_when_enough_free(self):
        smi = _fake_smi(self.tmp, 'echo "0, Tesla V100, 16384, 500, 7.0"')
        r = gpu.preflight(required_mib=1500, smi=smi)
        self.assertTrue(r.ok)
        self.assertAlmostEqual(r.best_free_mib, 15884)

    def test_refuses_when_free_memory_below_threshold(self):
        smi = _fake_smi(self.tmp, 'echo "0, Tesla V100, 16384, 15500, 7.0"')
        r = gpu.preflight(required_mib=1500, smi=smi)
        self.assertFalse(r.ok)
        self.assertIn("884 MiB free", r.reason)
        self.assertIn("shared", r.reason)

    def test_refuses_when_no_smi(self):
        r = gpu.preflight(smi=str(self.tmp / "nope"))
        self.assertFalse(r.ok)
        self.assertIn("not found", r.reason)
        self.assertEqual(r.devices, [])

    def test_picks_the_freest_device_of_several(self):
        smi = _fake_smi(self.tmp,
                        'echo "0, V100, 16384, 15000, 7.0"; '
                        'echo "1, V100, 16384, 100, 7.0"')
        r = gpu.preflight(required_mib=1500, smi=smi)
        self.assertTrue(r.ok)
        self.assertAlmostEqual(r.best_free_mib, 16284)

    def test_require_gpu_raises(self):
        smi = _fake_smi(self.tmp, 'echo "0, V100, 16384, 16300, 7.0"')
        with self.assertRaises(gpu.GpuUnavailable):
            gpu.require_gpu(required_mib=1500, smi=smi)


# --------------------------------------------------------------------------
def _make_system(tmp: Path, *, windows_done=0, legs=("complex", "solvent"),
                 controller: bool = True) -> Path:
    """A directory shaped like 6I5I_DUAL_FEP/."""
    sysdir = tmp / "SYS"
    for leg in legs:
        d = sysdir / leg
        d.mkdir(parents=True, exist_ok=True)
        for i in range(windows_done):
            (d / f".md_forward_w{i:02d}.done").touch()
    if controller:
        (sysdir / jobs.CONTROLLER).write_text("#!/bin/bash\necho fake controller\n")
    return sysdir


class TestState(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.runner = _IdleRunner()

    def test_counts_done_markers(self):
        s = _make_system(self.tmp, windows_done=6)
        st = state.read_state(s, runner=self.runner)
        fwd = st.legs["complex"].stages["md_forward"]
        self.assertEqual(fwd.done, 6)
        self.assertEqual(fwd.total, 15)
        self.assertFalse(fwd.complete)
        self.assertAlmostEqual(fwd.fraction, 6 / 15)

    def test_complete_when_all_windows_and_stages_done(self):
        s = _make_system(self.tmp, windows_done=15)
        for leg in ("complex", "solvent"):
            for stg in ("nvt_equil", "npt_equil", "md_backward"):
                (s / leg / f".{stg}.done").touch()
            for i in range(15):
                (s / leg / f".md_backward_w{i:02d}.done").touch()
        st = state.read_state(s, runner=self.runner)
        self.assertTrue(st.complete)
        self.assertFalse(st.running)

    def test_missing_dir_is_reported_not_raised(self):
        st = state.read_state(self.tmp / "nope", runner=self.runner)
        self.assertIsNotNone(st.error)
        self.assertEqual(st.legs, {})

    def test_stale_pidfile_does_not_look_running(self):
        s = _make_system(self.tmp, windows_done=1)
        (s / "controller_run.pid").write_text("999999999")   # cannot exist
        st = state.read_state(s, runner=self.runner)
        self.assertIsNone(st.controller_pid)
        self.assertFalse(st.running)

    def test_live_pidfile_is_detected(self):
        s = _make_system(self.tmp, windows_done=1)
        (s / "controller_run.pid").write_text(str(os.getpid()))
        st = state.read_state(s, runner=self.runner)
        self.assertEqual(st.controller_pid, os.getpid())
        self.assertTrue(st.running)

    def test_salvaged_single_run_windows_are_counted(self):
        """Windows completed by a crashed pre-controller run have no markers.

        Mirrors the real 6I5I complex/forward layout: the pre-controller single
        run completed w00-w05, then per-window runs covered w06-w14. The two
        sets are disjoint, so a marker-only count reports 9/15 for a leg that is
        actually finished.
        """
        s = _make_system(self.tmp, controller=False)
        for i in range(6, 15):                          # markers w06..w14 only
            (s / "complex" / f".md_forward_w{i:02d}.done").touch()
        fep = s / "complex" / "md_forward.fepout"
        # a single run that completed windows 0..5 (boundary rows at 250000..1500000)
        rows = []
        for i in range(6):
            for step in range(i * 250000 + 500, (i + 1) * 250000 + 1, 500):
                rows.append(f"FepEnergy: {step} 0 0 0 0 0.0 0.0 300.0 0.0")
        fep.write_text("\n".join(rows) + "\n")

        st = state.read_state(s, runner=self.runner)
        fwd = st.legs["complex"].stages["md_forward"]
        self.assertEqual(fwd.done, 15, "6 salvaged + 9 marked should be complete")
        self.assertTrue(fwd.complete)
        self.assertIn("salvaged", fwd.note)

    def test_partial_single_run_window_is_not_counted(self):
        """A window that did not reach its boundary step is incomplete."""
        s = _make_system(self.tmp, controller=False)
        fep = s / "complex" / "md_forward.fepout"
        rows = [f"FepEnergy: {step} 0 0 0 0 0.0 0.0 300.0 0.0"
                for step in range(500, 250001, 500)]          # window 0 complete
        rows += [f"FepEnergy: {step} 0 0 0 0 0.0 0.0 300.0 0.0"
                 for step in range(250500, 400001, 500)]      # window 1 truncated
        fep.write_text("\n".join(rows) + "\n")
        st = state.read_state(s, runner=self.runner)
        self.assertEqual(st.legs["complex"].stages["md_forward"].done, 1)

    def test_stage_finished_without_marker_is_detected(self):
        s = _make_system(self.tmp, controller=False)
        leg = s / "complex"
        (leg / "npt_equil.coor").write_text("x")
        (leg / "npt_equil.log").write_text("...\nWallClock: 1.0\nEnd of program\n")
        st = state.read_state(s, runner=self.runner)
        npt = st.legs["complex"].stages["npt_equil"]
        self.assertTrue(npt.complete)
        self.assertIn("without a marker", npt.note)

    def test_output_without_clean_log_is_not_claimed_done(self):
        s = _make_system(self.tmp, controller=False)
        leg = s / "complex"
        (leg / "npt_equil.coor").write_text("x")
        (leg / "npt_equil.log").write_text("FATAL ERROR: something\n")
        st = state.read_state(s, runner=self.runner)
        self.assertFalse(st.legs["complex"].stages["npt_equil"].complete)

    def test_subwindow_progress_advances_fraction(self):
        stg = state.StageState(name="md_forward", done=3, total=15, running=True,
                               rows_in_progress=state.ROWS_PER_WINDOW // 2)
        self.assertAlmostEqual(stg.fraction, 3.5 / 15)


# --------------------------------------------------------------------------
class _FakePopen:
    last_kwargs: dict = {}

    def __init__(self, *a, **kw):
        self.pid = 4242
        self.kw = kw
        _FakePopen.last_kwargs = kw


class TestJobs(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.runner = _IdleRunner()
        self.ok_smi = _fake_smi(self.tmp, 'echo "0, Tesla V100, 16384, 500, 7.0"',
                                name="smi-ok")
        self.busy_smi = _fake_smi(self.tmp,
                                  'echo "0, Tesla V100, 16384, 16300, 7.0"',
                                  name="smi-busy")

    def test_submit_refuses_on_bad_gpu(self):
        s = _make_system(self.tmp)
        with self.assertRaises(jobs.JobError) as cm:
            jobs.submit(s, smi=self.busy_smi, runner=self.runner)
        self.assertIn("GPU preflight failed", str(cm.exception))

    def test_submit_refuses_when_already_running(self):
        s = _make_system(self.tmp, windows_done=1)
        (s / "controller_run.pid").write_text(str(os.getpid()))
        with self.assertRaises(jobs.JobError) as cm:
            jobs.submit(s, smi=self.ok_smi, runner=self.runner)
        self.assertIn("already running", str(cm.exception))

    def test_submit_refuses_without_controller_script(self):
        s = self.tmp / "EMPTY"
        s.mkdir()
        with self.assertRaises(jobs.JobError) as cm:
            jobs.submit(s, smi=self.ok_smi, runner=self.runner)
        self.assertIn("not found", str(cm.exception))

    def test_submit_launches_detached_with_good_gpu(self):
        s = _make_system(self.tmp)
        res = jobs.submit(s, smi=self.ok_smi, runner=self.runner, popen=_FakePopen)
        self.assertTrue(res.started)
        self.assertEqual(res.pid, 4242)
        self.assertTrue(res.gpu.ok)
        self.assertIn("GPU ok", res.reason)
        # must be detached, or the job dies with the UI session
        self.assertTrue(_FakePopen.last_kwargs.get("start_new_session"))
        self.assertEqual(_FakePopen.last_kwargs.get("cwd"), str(s))

    def test_submit_can_skip_gpu_for_plumbing(self):
        s = _make_system(self.tmp)
        res = jobs.submit(s, smi=self.busy_smi, runner=self.runner,
                          popen=_FakePopen, allow_no_gpu=True)
        self.assertTrue(res.started)
        self.assertFalse(res.gpu.ok)
        self.assertIn("SKIPPED", res.reason)

    def test_cancel_when_idle_is_harmless(self):
        s = _make_system(self.tmp)
        self.assertIn("nothing", jobs.cancel(s, runner=self.runner))

    def test_logs_reads_tail(self):
        s = _make_system(self.tmp)
        (s / "controller_run.log").write_text("\n".join(f"line{i}" for i in range(50)))
        self.assertEqual(jobs.logs(s, tail=3), ["line47", "line48", "line49"])

    def test_logs_missing_file_is_empty(self):
        self.assertEqual(jobs.logs(_make_system(self.tmp), name="nope.log"), [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
