"""Flask UI over fep_web.

Session-scoped by design: you start it when the GPU box is up and it lives only
as long as that session. This is why there is no user database and no permanent
URL handling -- the requirement is "see and drive the job while I'm working",
not "a always-on service".

All job logic lives in `fep_web.jobs` / `state` / `gpu`; this module only renders
and routes, so the same backend can back a CLI or a different UI later.

Run:
    python3 -m fep_web.app --root . --port 8080 --host 0.0.0.0

Binding 0.0.0.0 exposes it to the network. If you reach it through a port proxy,
set a shared secret with --token (or FEP_WEB_TOKEN) -- there is no other auth.
"""
from __future__ import annotations

import argparse
import os
import secrets
from pathlib import Path

from flask import (Flask, abort, flash, jsonify, redirect, render_template,
                   request, url_for)

from . import jobs, results
from .gpu import preflight
from .state import WINDOWS, read_state

DEFAULT_PORT = 8080          # the port AI Studio's service proxy documents
MARKER = "fep_run.py"        # a directory with this is a prepared system


def discover_systems(root: Path) -> dict[str, Path]:
    """Prepared system directories, keyed by name. Sorted for stable output."""
    root = Path(root)
    if not root.is_dir():
        return {}
    return {d.name: d for d in sorted(root.iterdir())
            if d.is_dir() and (d / MARKER).exists()}


def create_app(root: Path, token: str | None = None, smi: str = "nvidia-smi") -> Flask:
    root = Path(root).resolve()
    app = Flask(__name__)
    app.secret_key = secrets.token_hex(16)      # flash messages only; ephemeral
    app.config["SYSTEMS_ROOT"] = root
    app.config["TOKEN"] = token
    app.config["SMI"] = smi

    # -- helpers ----------------------------------------------------------
    def systems() -> dict[str, Path]:
        return discover_systems(app.config["SYSTEMS_ROOT"])

    def system_or_404(name: str) -> Path:
        # Validate against the discovered set rather than touching the name,
        # so a caller cannot escape the systems root via ../ or a symlink.
        found = systems().get(name)
        if found is None:
            abort(404, description=f"no system named {name!r}")
        return found

    def gpu():
        return preflight(smi=app.config["SMI"])

    # -- auth -------------------------------------------------------------
    @app.before_request
    def _check_token():
        want = app.config["TOKEN"]
        if not want:
            return None
        supplied = (request.headers.get("X-FEP-Token")
                    or request.args.get("token"))
        if supplied != want:
            abort(401, description="missing or wrong token")
        return None

    # -- views ------------------------------------------------------------
    @app.route("/")
    def index():
        rows = [{"name": name, "path": path, **_summarise(path)}
                for name, path in systems().items()]
        return render_template("index.html", gpu=gpu(), rows=rows)

    @app.route("/system/<name>")
    def system_detail(name):
        path = system_or_404(name)
        return render_template("system.html", gpu=gpu(), name=name, path=path,
                               st=read_state(path), windows=WINDOWS,
                               readiness=results.readiness(path),
                               log_lines=jobs.logs(path, tail=200))

    @app.route("/system/<name>/submit", methods=["POST"])
    def system_submit(name):
        path = system_or_404(name)
        try:
            res = jobs.submit(path, smi=app.config["SMI"])
            flash(res.reason, "ok")
        except jobs.JobError as exc:
            flash(str(exc), "error")
        return redirect(url_for("system_detail", name=name))

    @app.route("/system/<name>/cancel", methods=["POST"])
    def system_cancel(name):
        path = system_or_404(name)
        try:
            flash(jobs.cancel(path), "ok")
        except jobs.JobError as exc:
            flash(str(exc), "error")
        return redirect(url_for("system_detail", name=name))

    @app.route("/system/<name>/result")
    def system_result(name):
        path = system_or_404(name)
        try:
            return jsonify(results.bar_result(path))
        except results.AnalysisError as exc:
            return jsonify({"error": str(exc)}), 409

    # -- json API (used by the pages' polling) ----------------------------
    @app.route("/api/gpu")
    def api_gpu():
        g = gpu()
        return jsonify({"ok": g.ok, "reason": g.reason,
                        "best_free_mib": g.best_free_mib,
                        "devices": [d.__dict__ for d in g.devices]})

    @app.route("/api/systems")
    def api_systems():
        return jsonify({name: _state_dict(read_state(p))
                        for name, p in systems().items()})

    @app.route("/api/system/<name>")
    def api_system(name):
        return jsonify(_state_dict(read_state(system_or_404(name))))

    @app.route("/api/system/<name>/log")
    def api_log(name):
        system_or_404(name)
        tail = request.args.get("tail", default=120, type=int)
        return jsonify({"lines": jobs.logs(systems()[name], tail=tail)})

    return app


def _summarise(path: Path) -> dict:
    """Roll a system's stages up into one progress fraction for the dashboard."""
    st = read_state(path)
    stages = [s for ls in st.legs.values() for s in ls.stages.values()]
    done = sum(s.done for s in stages)
    total = sum(s.total for s in stages)
    return {
        "state": st,
        "done": done,
        "total": total,
        "fraction": (done / total) if total else 0.0,
    }


def _state_dict(st) -> dict:
    """SystemState -> JSON-serialisable dict."""
    return {
        "path": str(st.path),
        "error": st.error,
        "running": st.running,
        "complete": st.complete,
        "controller_pid": st.controller_pid,
        "namd_pids": st.namd_pids,
        "active_window": st.active_window,
        "legs": {
            leg: {name: {"done": s.done, "total": s.total,
                         "complete": s.complete, "running": s.running,
                         "fraction": round(s.fraction, 4), "note": s.note}
                  for name, s in ls.stages.items()}
            for leg, ls in st.legs.items()
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="NAMD RBFE web UI")
    ap.add_argument("--root", type=Path, default=Path.cwd(),
                    help="directory containing the system folders")
    ap.add_argument("--host", default="0.0.0.0")
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--token", default=os.environ.get("FEP_WEB_TOKEN"),
                    help="shared secret; required if you expose this beyond "
                         "localhost (env FEP_WEB_TOKEN also works)")
    ap.add_argument("--smi", default="nvidia-smi")
    a = ap.parse_args(argv)

    app = create_app(a.root, token=a.token, smi=a.smi)
    found = discover_systems(Path(a.root).resolve())
    print(f"NAMD RBFE UI — {len(found)} system(s) under {a.root}: "
          f"{', '.join(found) or '(none)'}")
    print(f"listening on http://{a.host}:{a.port}"
          + ("  [token required]" if a.token else "  [NO AUTH]"))
    app.run(host=a.host, port=a.port, debug=False, threaded=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
