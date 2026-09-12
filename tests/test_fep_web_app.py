"""Tests for the Flask layer — no GPU, no NAMD, no network.

Run:  python3 -m unittest discover -s tests -v

The controller script written into each fake system is a one-line `exit 0`, so
POSTing /submit exercises the real detached-launch path (it does spawn a bash
that immediately exits) without ever starting NAMD.
"""
from __future__ import annotations

import json
import stat
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fep_web import app as webapp            # noqa: E402

SYSTEM = "TESTSYS"


def _fake_smi(tmp: Path, body: str, name: str = "nvidia-smi") -> str:
    """Write an executable stand-in for nvidia-smi.

    `name` must differ per fake: two fakes sharing a filename means the second
    silently overwrites the first, and both paths then report the same GPU.
    """
    p = tmp / name
    p.write_text("#!/bin/sh\n" + body + "\n")
    p.chmod(p.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return str(p)


def _make_root(tmp: Path, *, windows_done=0, name: str = SYSTEM) -> Path:
    root = tmp / "root"
    sysdir = root / name
    for leg in ("complex", "solvent"):
        d = sysdir / leg
        d.mkdir(parents=True, exist_ok=True)
        for i in range(windows_done):
            (d / f".md_forward_w{i:02d}.done").touch()
    (sysdir / "fep_run.py").write_text("# marker: makes this a prepared system\n")
    (sysdir / "analyze_fep.py").write_text("# presence is what readiness() checks\n")
    (sysdir / "respawn_controller.sh").write_text("#!/bin/bash\nexit 0\n")
    (sysdir / "controller_run.log").write_text("line-a\nline-b\nline-c\n")
    return root


class AppCase(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.ok_smi = _fake_smi(self.tmp, 'echo "0, Tesla V100, 16384, 500, 7.0"',
                                name="smi-ok")
        self.busy_smi = _fake_smi(self.tmp,
                                  'echo "0, Tesla V100, 16384, 16300, 7.0"',
                                  name="smi-busy")

    def client(self, root, **kw):
        kw.setdefault("smi", self.ok_smi)
        app = webapp.create_app(root, **kw)
        app.config["TESTING"] = True
        return app.test_client()


class TestDiscovery(unittest.TestCase):
    def test_finds_only_dirs_with_fep_run_py(self):
        tmp = Path(tempfile.mkdtemp())
        root = _make_root(tmp)
        (root / "not-a-system").mkdir()
        found = webapp.discover_systems(root)
        self.assertEqual(list(found), [SYSTEM])

    def test_missing_root_is_empty_not_an_error(self):
        self.assertEqual(webapp.discover_systems(Path("/nope/nope")), {})


class TestPages(AppCase):
    def test_index_with_no_systems(self):
        empty = self.tmp / "empty"
        empty.mkdir()
        r = self.client(empty).get("/")
        self.assertEqual(r.status_code, 200)
        self.assertIn(b"No prepared systems", r.data)

    def test_index_lists_a_system(self):
        r = self.client(_make_root(self.tmp, windows_done=3)).get("/")
        self.assertEqual(r.status_code, 200)
        self.assertIn(SYSTEM.encode(), r.data)

    def test_system_detail_renders_stages(self):
        r = self.client(_make_root(self.tmp, windows_done=3)).get(f"/system/{SYSTEM}")
        self.assertEqual(r.status_code, 200)
        for stage in (b"nvt_equil", b"npt_equil", b"md_forward", b"md_backward"):
            self.assertIn(stage, r.data)

    def test_unknown_system_is_404(self):
        r = self.client(_make_root(self.tmp)).get("/system/NOPE")
        self.assertEqual(r.status_code, 404)


class TestPathSafety(AppCase):
    def test_traversal_is_refused(self):
        root = _make_root(self.tmp)
        client = self.client(root)
        for attempt in ("..", "../..", "%2e%2e", "....//", "/etc"):
            r = client.get(f"/system/{attempt}")
            self.assertIn(r.status_code, (404, 308),
                          f"{attempt!r} returned {r.status_code}")

    def test_absolute_path_is_refused(self):
        r = self.client(_make_root(self.tmp)).get("/system//etc/passwd")
        self.assertIn(r.status_code, (404, 308))


class TestAuth(AppCase):
    def test_token_required_when_set(self):
        client = self.client(_make_root(self.tmp), token="s3cret")
        self.assertEqual(client.get("/").status_code, 401)

    def test_token_accepted_via_query(self):
        client = self.client(_make_root(self.tmp), token="s3cret")
        self.assertEqual(client.get("/?token=s3cret").status_code, 200)

    def test_token_accepted_via_header(self):
        client = self.client(_make_root(self.tmp), token="s3cret")
        r = client.get("/", headers={"X-FEP-Token": "s3cret"})
        self.assertEqual(r.status_code, 200)

    def test_wrong_token_rejected(self):
        client = self.client(_make_root(self.tmp), token="s3cret")
        self.assertEqual(client.get("/?token=nope").status_code, 401)

    def test_no_token_configured_is_open(self):
        self.assertEqual(self.client(_make_root(self.tmp)).get("/").status_code, 200)

    def test_token_survives_navigation_via_session(self):
        """In-page links cannot carry the token, so auth must stick.

        Without this, clicking from the dashboard into a system would 401 and
        the UI would be unusable behind a token.
        """
        client = self.client(_make_root(self.tmp), token="s3cret")
        self.assertEqual(client.get("/?token=s3cret").status_code, 200)
        # subsequent requests carry only the session cookie
        self.assertEqual(client.get("/").status_code, 200)
        self.assertEqual(client.get(f"/system/{SYSTEM}").status_code, 200)
        self.assertEqual(client.get("/api/gpu").status_code, 200)

    def test_session_is_not_shared_with_unauthenticated_client(self):
        root = _make_root(self.tmp)
        good = self.client(root, token="s3cret")
        good.get("/?token=s3cret")
        stranger = self.client(root, token="s3cret")
        self.assertEqual(stranger.get("/").status_code, 401)


class TestApi(AppCase):
    def test_gpu_endpoint(self):
        d = json.loads(self.client(_make_root(self.tmp)).get("/api/gpu").data)
        self.assertTrue(d["ok"])
        self.assertEqual(d["devices"][0]["name"], "Tesla V100")

    def test_systems_endpoint_shape(self):
        d = json.loads(self.client(_make_root(self.tmp, windows_done=2))
                       .get("/api/systems").data)
        self.assertIn(SYSTEM, d)
        self.assertIn("legs", d[SYSTEM])
        self.assertEqual(d[SYSTEM]["legs"]["complex"]["md_forward"]["total"], 15)
        self.assertEqual(d[SYSTEM]["legs"]["complex"]["md_forward"]["done"], 2)

    def test_log_endpoint(self):
        d = json.loads(self.client(_make_root(self.tmp))
                       .get(f"/api/system/{SYSTEM}/log?tail=2").data)
        self.assertEqual(d["lines"], ["line-b", "line-c"])


class TestResult(AppCase):
    def test_result_explains_what_is_missing(self):
        r = self.client(_make_root(self.tmp)).get(f"/system/{SYSTEM}/result")
        self.assertEqual(r.status_code, 409)
        d = json.loads(r.data)
        self.assertIn("not assembled yet", d["error"])

    def test_detail_page_says_not_ready(self):
        r = self.client(_make_root(self.tmp)).get(f"/system/{SYSTEM}")
        self.assertIn(b"not assembled yet", r.data)


class TestControls(AppCase):
    def test_submit_redirects_and_logs(self):
        root = _make_root(self.tmp)
        client = self.client(root)
        r = client.post(f"/system/{SYSTEM}/submit")
        self.assertEqual(r.status_code, 302)
        log = (root / SYSTEM / "controller_run.log").read_text()
        self.assertIn("launched by fep_web", log)

    def test_submit_blocked_when_gpu_is_busy(self):
        root = _make_root(self.tmp)
        client = self.client(root, smi=self.busy_smi)
        r = client.post(f"/system/{SYSTEM}/submit", follow_redirects=True)
        self.assertIn(b"GPU preflight failed", r.data)
        self.assertNotIn("launched by fep_web",
                         (root / SYSTEM / "controller_run.log").read_text())

    def test_cancel_when_idle_is_reported_not_raised(self):
        r = self.client(_make_root(self.tmp)).post(f"/system/{SYSTEM}/cancel",
                                                   follow_redirects=True)
        self.assertEqual(r.status_code, 200)
        self.assertIn(b"nothing was running", r.data)


if __name__ == "__main__":
    unittest.main(verbosity=2)
