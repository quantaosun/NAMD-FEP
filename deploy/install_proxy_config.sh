#!/usr/bin/env bash
# Install / restore the jupyter-server-proxy registration for the NAMD RBFE UI.
#
#   bash deploy/install_proxy_config.sh          # install
#   bash deploy/install_proxy_config.sh --check  # verify, change nothing
#
# RUN THIS AFTER EVERY CONTAINER RECYCLE. ~/.jupyter/ is not persisted: the
# 2026-09-13 recycle replaced it with the platform default and silently removed
# the named proxy server, so /fep/ stopped resolving. The app itself was fine
# the whole time -- only this registration was lost.
#
# The file it installs is inert until the Jupyter server restarts (handlers are
# registered at startup), so installing is safe while a session is live; the
# route appears on the next WebIDE restart.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
SRC="$HERE/jupyter_server_config.py"
TARGET="$HOME/.jupyter/jupyter_server_config.py"
PYTHON="/opt/conda/envs/webide/bin/python3.7"

CHECK_ONLY=0
[[ "${1:-}" == "--check" ]] && CHECK_ONLY=1

ok()   { printf '  \033[32m✓\033[0m %s\n' "$1"; }
bad()  { printf '  \033[31m✗\033[0m %s\n' "$1"; }
warn() { printf '  \033[33m!\033[0m %s\n' "$1"; }
die()  { bad "$1"; exit 1; }

echo "NAMD RBFE UI — jupyter-server-proxy registration"
echo "  repo:   $REPO"
echo "  target: $TARGET"
echo

# -- preconditions ------------------------------------------------------------
[[ -x "$PYTHON" ]] || die "interpreter not found: $PYTHON"
ok "interpreter $PYTHON"

# The proxy runs the app with this interpreter, so IT -- not the shell's
# python3 -- must be able to import both the extension and the app.
if "$PYTHON" -c "import jupyter_server_proxy" 2>/dev/null; then
  ver="$("$PYTHON" -c 'import jupyter_server_proxy as m; print(m.__version__)' 2>/dev/null || echo '?')"
  ok "jupyter_server_proxy $ver (via $PYTHON)"
else
  die "jupyter_server_proxy is not importable by $PYTHON — the WebIDE needs it installed"
fi

if "$PYTHON" -c "import sys; sys.path.insert(0, '$REPO'); import fep_web.app" 2>/dev/null; then
  ok "fep_web.app imports (flask present)"
else
  die "fep_web.app does not import under $PYTHON — run: $PYTHON -c 'import sys; sys.path.insert(0,\"$REPO\"); import fep_web.app'"
fi

if "$PYTHON" -m jupyter serverextension list 2>/dev/null | grep -q 'jupyter_server_proxy.*enabled'; then
  ok "extension enabled"
else
  warn "jupyter_server_proxy does not appear as 'enabled' — check the WebIDE config dirs"
fi

# -- install ------------------------------------------------------------------
if (( CHECK_ONLY )); then
  if [[ -f "$TARGET" ]]; then ok "installed: $TARGET"; else bad "MISSING: $TARGET"; fi
  exit 0
fi

mkdir -p "$(dirname "$TARGET")"

if [[ -f "$TARGET" ]] && ! cmp -s <(sed "s#^REPO = .*#REPO = \"$REPO\"#" "$SRC") "$TARGET"; then
  backup="$TARGET.bak.$(date +%Y%m%d-%H%M%S)"
  cp -p "$TARGET" "$backup"
  warn "existing config backed up to $backup"
fi

# The committed file carries the canonical repo path; rewrite it if this
# checkout lives elsewhere so PYTHONPATH and --root stay correct.
sed "s#^REPO = .*#REPO = \"$REPO\"#" "$SRC" > "$TARGET"
chmod 644 "$TARGET"
ok "wrote $TARGET"

# Fail loudly here rather than with a 404 later: a config that does not parse
# takes the whole Jupyter server down with it on restart.
if "$PYTHON" - "$TARGET" <<'PY'
import builtins, runpy, sys, types
# config files call get_config(), which jupyter injects into their namespace.
# runpy executes in a fresh globals dict, so stub it via builtins instead.
builtins.get_config = lambda: types.SimpleNamespace(
    ServerProxy=types.SimpleNamespace(servers=None))
ns = runpy.run_path(sys.argv[1])
servers = ns["c"].ServerProxy.servers
assert "fep" in servers, servers
cmd = servers["fep"]["command"]
assert isinstance(cmd, list) and "{port}" in cmd, cmd
assert servers["fep"]["environment"]["FEP_WEB_TOKEN"] == ""
print("    command: " + " ".join(cmd))
PY
then
  ok "config parses; server 'fep' well-formed"
else
  die "installed config failed validation — restoring is advised"
fi

# -- what next ----------------------------------------------------------------
cat <<EOF

Next: RESTART the WebIDE / Jupyter server. The proxy reads this config at
startup, so the route does not exist until then.

After the restart, open (note: no /proxy/):

    <your WebIDE URL>/bj-cpu-01/user/696173/10639452/fep/

i.e. take the URL in your browser's address bar and replace its last path
segment with 'fep/'. If it 404s, confirm nothing is squatting port 8080
(ps -eo pid,cmd | grep -E 'serve_ui|fep_web.app') -- a hand-started server
makes the proxy's spawn fail, which surfaces as a bare 500.
EOF
