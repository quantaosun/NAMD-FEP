"""jupyter-server-proxy registration for the NAMD RBFE web UI.

INSTALL:  bash deploy/install_proxy_config.sh
          (copies this file to ~/.jupyter/jupyter_server_config.py)

This is the file that makes the UI reachable in a browser. It is kept *inside
the repo* on purpose: the previous copy existed only at
~/.jupyter/jupyter_server_config.py, so the 2026-09-13 container recycle
destroyed it and took the deployment down with it. Nothing about the working
setup was recoverable from git.

Browser URL:  <your WebIDE URL>/.../<server-id>/fep/
              NOT /proxy/fep/ -- this registers a *named* server, which
              jupyter-server-proxy mounts at ujoin(base_url, name).
              See RUNBOOK section 10.

Why a proxy server rather than a tunnel or the AI Studio deploy button: see
RUNBOOK section 10. Short version -- outbound 7844 (cloudflared) and 22 (SSH
tunnels) are both blocked, and the deploy button ships the app to a GPU-less
container that cannot see /home/aistudio at all, so it could not drive a job
even if it did deploy.

Traps encoded below, each of which has already bitten:

  * ``command`` must be a LIST -- the proxy execs it directly, with no shell.
  * There is NO ``cwd`` option (jupyter_server_proxy/config.py,
    ``make_server_process``), so the repo is made importable via PYTHONPATH.
  * ``{port}`` is substituted by ``str.format`` (config.py,
    ``process_args``/``_render_template``). The proxy picks a free port and
    waits for the process to bind it. Braces are single, not doubled -- the
    ``{{port}}`` in the upstream docstring is docstring escaping, not syntax.
  * Bind 127.0.0.1, not 0.0.0.0. The proxy is then the only way in, and it sits
    behind JupyterHub's auth. Binding 0.0.0.0 would put the UI on the LAN with
    no auth of its own.
  * FEP_WEB_TOKEN is forced empty. The app reads an empty token as "no auth"
    (fep_web/app.py: ``if not want: return None``). Without this line a token
    inherited from the Jupyter server's environment would silently start
    demanding a ``?token=`` the browser does not have.
"""

c = get_config()  # noqa: F821  -- injected by jupyter into config namespaces

REPO = "/home/aistudio/work/NAMD-FEP"
PYTHON = "/opt/conda/envs/webide/bin/python3.7"

c.ServerProxy.servers = {
    "fep": {
        "command": [
            PYTHON, "-m", "fep_web.app",
            "--root", REPO,
            "--host", "127.0.0.1",
            "--port", "{port}",
        ],
        "environment": {
            "PYTHONPATH": REPO,
            "FEP_WEB_TOKEN": "",
        },
        # No launcher icon. The UI is reached by URL, and staying out of the
        # Lab launcher discourages a second, hand-started instance -- which
        # squats the port the proxy wants and turns /fep/ into a bare 500.
        "launcher_entry": {"enabled": False},
        # Seconds to wait for the app to bind its port. This is a STARTUP
        # readiness timeout, not a request timeout -- per-request is hardcoded
        # at 300s (handlers.py:506). The app starts in about a second.
        "timeout": 30,
    }
}
