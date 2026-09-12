#!/usr/bin/env python3
"""AI Studio deployment entry point for the NAMD RBFE web UI.

AI Studio's service deployment ("未创建过服务 ... 在 Python 文件中编写服务代码，
点击顶部部署按钮进行部署") expects a Python file *in the project* holding the
service code. This is that file. Open it in the AI Studio editor and click 部署.

It exposes the WSGI callable as a module-level `app` (for platforms that import
it) and also runs a server directly when executed.

    python3 serve_ui.py                 # token auto-generated, shown at startup
    FEP_WEB_TOKEN=xyz python3 serve_ui.py
    python3 serve_ui.py --port 9000

Port/host follow Baidu's documented requirements for service deployment:
bind 0.0.0.0 and listen on 8080. The app is also prefix-tolerant, so it works
whether or not the /api_serving/<port> prefix is stripped by the proxy.
"""
from __future__ import annotations

import argparse
import os
import secrets
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from fep_web.app import create_app, discover_systems   # noqa: E402

TOKEN_FILE = Path.home() / ".fep_web_token"
DEFAULT_HOST = "0.0.0.0"
DEFAULT_PORT = 8080


def resolve_token(explicit: str | None = None) -> str:
    """env var -> stored file -> freshly generated (and stored).

    Persisted so the URL you bookmark keeps working across redeploys; /tmp is
    not a safe place for it because the container is recycled.
    """
    if explicit:
        return explicit
    env = os.environ.get("FEP_WEB_TOKEN")
    if env:
        return env
    try:
        existing = TOKEN_FILE.read_text().strip()
        if existing:
            return existing
    except OSError:
        pass
    token = secrets.token_urlsafe(12)
    try:
        TOKEN_FILE.write_text(token + "\n")
        TOKEN_FILE.chmod(0o600)
    except OSError:
        pass                      # read-only home: the token still works this run
    return token


TOKEN = resolve_token()

# Module-level WSGI callable, for a platform that imports this module.
app = create_app(ROOT, token=TOKEN)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="NAMD RBFE UI (AI Studio entry point)")
    ap.add_argument("--host", default=DEFAULT_HOST)
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--root", type=Path, default=ROOT,
                    help="directory containing the system folders")
    ap.add_argument("--token", default=None)
    ap.add_argument("--no-token", action="store_true",
                    help="disable auth (only if the proxy already authenticates)")
    a = ap.parse_args(argv)

    global app
    token = None if a.no_token else resolve_token(a.token)
    if token != TOKEN:
        app = create_app(a.root, token=token)

    found = discover_systems(Path(a.root).resolve())
    print(f"NAMD RBFE UI — {len(found)} system(s) under {a.root}: "
          f"{', '.join(found) or '(none)'}", flush=True)
    if token:
        print(f"token: {token}", flush=True)
        print(f"open:  http://<host>:{a.port}/?token={token}", flush=True)
    else:
        print("AUTH DISABLED", flush=True)
    app.run(host=a.host, port=a.port, debug=False, threaded=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
