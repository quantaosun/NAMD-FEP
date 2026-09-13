# RUNBOOK — operating the NAMD RBFE workflow

Operational guide for the current work. The root `README.md` is the original
2021-era public tutorial (LigParGen + Feprepare) and does **not** describe this
pipeline. `NAMD_RBFE_Guide_zh.md` is the detailed Chinese walkthrough and is kept
in sync with this file.

---

## 1. Repo map — what is current and what is history

| path | status |
|---|---|
| **`6I5I_DUAL_FEP/`** | **CURRENT.** The system that produced the reported ΔΔG. Prep, configs, per-window results, controller, analysis. |
| `toppar/` | **CURRENT.** CHARMM36 parameters. Includes the project-specific `par_water_ions_clean.prm`. |
| `fep_pipeline/` | ⚠️ **NOT WORKING.** Never ran end-to-end; `run` path is broken (see §7). Only `alignment.py` is usable. |
| `T4L_RBFE/` | 📁 historical. NAMD tutorial system, built with the Feprepare/CHARMM-GUI/LigParGen web services. |
| `6I5I_RBFE/` | 📁 historical. 6I5I attempt via acpype + GAFF. |
| `6I5I_FEP/` | 📁 historical. 6I5I via acpype/GAFF + RDKit atom mapping (`atom_map.json`). Source of the ligand PDBs/SDFs still cited in `DDG_preliminary.md`. |
| `6I5I_deck/` | 📁 presentation slides (reveal.js). |
| `*.namd` at root, `NAMD-FEP-tutorial.pdf`, `NAMD-FEP_local.ipynb` | 📁 historical tutorial material. |

Anything marked historical is kept for provenance. **Do not build on it.**

---

## 2. Quick start — the one entry point

Everything runs from `6I5I_DUAL_FEP/`.

```bash
cd 6I5I_DUAL_FEP

# 1. Check the GPU first (a job started without one wastes the setup time)
nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv

# 2. Start / resume the whole job. Idempotent — safe to re-run any number of times.
bash respawn_controller.sh     # detached, survives the box's ~4.5 h job kills
```

`respawn_controller.sh` → `run_checkpointed.sh` → `fep_run.py`. It is
marker-guarded: every λ-window writes a `.done` marker, so re-running resumes at
the interrupted window instead of replaying.

**Watch progress:**
```bash
tail -f controller_run.log          # controller decisions
bash checkpoint_status.sh &         # one status line every 5 min
python3 fep_run.py windows complex md_forward   # JSON: how much is done
```

### ⚠️ `run_all.sh` is NOT the entry point

It is the original non-resumable launcher: it replays every stage from
`nvt_equil` and overwrites `md_forward`/`md_backward`. It now **refuses to run**
when production state exists, and exits 1. Use `run_checkpointed.sh`.
(`run_all.sh --force` overrides — only for a genuinely fresh start.)

---

## 3. The full pipeline, in order

```bash
cd 6I5I_DUAL_FEP

# --- Preparation (only needed once; regenerates hybrid/, ionized.*, *.namd) ---
python3 prepare_hybrid.py      # inputs/{ref,mut}.{rtf,prm,pdb} -> hybrid/hybrid.{rtf,prm,pdb}
python3 build_system.py        # VMD psfgen + solvate + ionize -> complex/, solvent/
python3 write_fep_inputs.py    # ionized.fep B-factors + fep.tcl + the 8 .namd configs

# --- Equilibration (per leg: nvt then npt) ---
# --- Production: one NAMD process per lambda window, forward then backward ---
#     all handled by the controller; do not launch these by hand
```

Each leg is `nvt_equil → npt_equil → md_forward (0→1) → md_backward (1→0)`.
15 λ-windows per direction, 50 000 equilibration + 250 000 production steps each.

> **Regenerating preparation files desyncs the committed results.** `hybrid/*`
> are checked-in *build products* of `prepare_hybrid.py`. If you re-run step 1,
> `hybrid/` no longer matches the run that produced the reported ΔΔG.

---

## 4. Analysis

```bash
cd 6I5I_DUAL_FEP

python3 audit_fep.py              # full audit: parser validation, BAR, hysteresis,
                                  # stationarity, overlap, block-bootstrap error bars
python3 audit_fep.py --skip-boot  # fast
python3 audit_fep.py --scan       # ΔΔG vs equilibration-trim sensitivity

python3 analyze_fep.py bar \
  complex/md_forward_combined.fepout complex/md_backward_combined.fepout \
  solvent/md_forward_combined.fepout solvent/md_backward_combined.fepout
```

**Current result:** ΔΔG(12H → desmethyl) = **−0.106 ± 0.120 kcal/mol**,
95% CI **[−0.355, +0.118]** — i.e. zero within error. See
`DDG_preliminary.md` §7 for the corrections that produced this number and the
known limitations.

> ⚠️ **Every per-window `.fepout` opens with 99 pre-equilibration rows** out of 500
> (`alchEquilSteps 50000` ÷ `alchOutFreq 500`). NAMD's own `dE_avg`/`dG` columns
> exclude them, and so must any analysis. Both scripts above trim by default;
> `--no-trim` reproduces the old (biased) behaviour.
>
> ⚠️ **Do not use `fep_pipeline.analysis`** — deleted 2026-09-12. Its greedy regex
> captured the cumulative `net change until now` instead of the per-window value.

---

## 5. Reassembling a combined fepout

`<stage>_combined.fepout` is the analysis input, built from the per-window files.
Rebuild with:

```bash
python3 fep_run.py assemble complex md_forward    # and md_backward, solvent/*
```

`assemble` is deliberately strict — it refuses unless every window is both
**complete (500 rows)** and **marked `.done`**, because a killed window leaves a
partial `.fepout` that would otherwise be assembled as a silently truncated
window. It preserves NAMD's own `STEPS OF EQUILIBRATION` and
`Free energy change ...` comment lines; those are what make the output
verifiable against NAMD's internal estimator.

---

## 6. Traps that have already bitten

| trap | what happens | guard now |
|---|---|---|
| `--dry` marking windows done | every later real run skips them; the controller declares the leg complete | removed — `--dry` writes no marker |
| assembling a partial window | silent truncation, biased ΔΔG | `assemble` requires the `.done` marker + 500 rows |
| `fep_run.py run <leg> fwd --start N` with no `--from` | window N is seeded from `npt_equil` instead of w(N−1) | pass `--from window_snapshots/md_forward_w05` explicitly |
| running `run_all.sh` on a finished job | ~14 h of duplicated GPU work, clobbered `.fepout` | refuses unless `--force` |
| κ `temperature` in a restart-chained config | `FATAL ERROR: Cannot specify both an initial temperature and a velocity file` | only `nvt_equil` sets it |
| a second `namd3` starting | splits the single GPU — makes everything *slower* | `wait_gpu_free` gate before every launch |

---

## 7. GPU notes (this box)

- **V100, sm_70, 16384 MiB usable — and shared.** Other processes may already hold
  memory. Check `nvidia-smi` before starting, not just that a GPU exists.
- Run NAMD with **`+p1`**, never `+p8`. In GPU-resident mode throughput scales
  *inversely* with PE count (`+p1` 75.8 ns/day vs `+p8` 13.2). `run_all.sh` and
  `run_checkpointed.sh` both set this; don't "optimise" it back.
- One GPU ⇒ **strictly sequential legs**. `wait_gpu_free` enforces this.
- NAMD peaks at ~590 MiB, so memory is not the constraint — latency is.
- Judge throughput from `TIMING:` sec/step, **not** the cumulative `PERFORMANCE:`
  line (a leading minimization drags that average down for hours).

Binary: `/home/aistudio/NAMD_3.0.3_Source/Linux-x86_64-g++.gpuresident/namd3`
VMD: `/home/aistudio/vmd-env/bin/vmd` (not on `PATH` — hardcoded in `build_system.py`)

---

## 8. `fep_pipeline/` — why it is marked broken

A 2643-line package advertised as a Feprepare replacement. It has **never run
end-to-end** (no output, no tests, no history). Three independent breaks in `run`:

1. `system_builder.py:60-62` passes `.prm` files and the hybrid **PDB** to
   psfgen's `topology` and emits no `parameters` lines — it cannot build a PSF.
2. `prepare_fep.py:339` references `alchFile ionized_{leg}.fep`; nothing generates it.
3. `fep_file.build_system_fep` — the function that would — is never called.

Salvageable: `alignment.py` (pure-stdlib Kabsch + B-factor marking) and the
`namd_config.py` templates (which still carry the restart-chain `temperature` bug).
Its `analysis.py` was deleted for returning wrong numbers.

---

## 9. `fep_web/` — job control as a library (foundation for the UI)

A host-agnostic layer over the scripts above, so the same logic can back a CLI,
a Flask app, or the existing shell entry points. It knows nothing about HTTP.

```python
from pathlib import Path
from fep_web import preflight, status, submit, cancel, logs

preflight().ok                     # is there a GPU with room RIGHT NOW?
st = status(Path("6I5I_DUAL_FEP")) # per-leg/per-stage progress, running PIDs
submit(Path("6I5I_DUAL_FEP"))      # GPU-preflights, then launches detached
cancel(Path("6I5I_DUAL_FEP"))      # SIGTERM controller + namd3
logs(Path("6I5I_DUAL_FEP"), tail=50)
```

`submit()` refuses if a job is already live, and the launch is detached
(`start_new_session=True`) so it survives a closed browser/session — which is
what makes the multi-day, several-hours-a-day pattern work.

**Why it exists rather than just calling the shell scripts:** the scripts assume
a free GPU. This V100 is *shared*, so `preflight()` checks free memory against a
1500 MiB threshold (NAMD peaks at ~590 MiB). It also reads job state correctly in
two cases a naive `.done`-marker check gets wrong — see `CLAUDE.md`.

```bash
python3 -m unittest discover -s tests    # 28 tests, no GPU or NAMD required
```

### 9.1 The web UI

```bash
cd /home/aistudio/work/NAMD-FEP
python3 -m fep_web.app --root . --port 8080 --host 0.0.0.0
```

Then open `http://<host>:8080`. It lists every directory containing `fep_run.py`,
shows GPU readiness, per-stage progress, a Start/Cancel control, the controller
log (auto-refreshing while a job runs), and a **Compute ΔΔG** button that runs
the BAR analysis and reports `ddG`, `dG_complex` and `dG_solvent`.

**Session-scoped by design** — start it when the GPU box is up; it does not need
to survive a restart, because the *job* does (detached launch + `.done` markers).

⚠️ **There is no authentication unless you set one.** Binding `0.0.0.0` exposes
it to whatever can reach the port. If you reach it through a proxy:

```bash
python3 -m fep_web.app --root . --port 8080 --token "$(openssl rand -hex 8)"
# then: http://<host>:8080/?token=<that>   (or send X-FEP-Token)
```

The GPU check gates the Start button: it is disabled and the reason shown when
the card has no room, so you cannot accidentally queue a job onto a contended
GPU.

---

## 10. Exposing the UI — SOLVED 2026-09-12 (jupyter-server-proxy)

**The UI is reachable in a browser at `<your WebIDE URL>/…/<server-id>/fep/`**
— no tunnel, no deploy button, no public IP needed. Verified working end to end
(HTTP 200, correct `<title>NAMD RBFE — systems</title>`, assets correctly
prefixed).

### The URL

JupyterHub's service prefix is `/bj-cpu-01/user/696173/10639452/`, so on this
container the path is:

```
https://<ai-studio-webide-host>/bj-cpu-01/user/696173/10639452/fep/
```

Take whatever URL is in your browser address bar when the WebIDE/Jupyter is
open and **replace its last path segment with `fep/`** (e.g. `…/10639452/lab`
→ `…/10639452/fep/`). The public host is terminated by AI Studio's ingress and
is **not discoverable from inside the box** — every hostname in the logs
(`jupyter-696173-10639452:8888`, `127.0.0.1:8888`) is internal, so read the host
out of your own browser rather than guessing it.

⚠️ **Not `/proxy/fep/`.** A *named* server is mounted at
`ujoin(base_url, sp.name, …)` → `/<base>/fep/`. `/proxy/<name>/` is for
jupyter-server-proxy's *unnamed* (arbitrary-port) form and 404s here.
`/fep/` with no base prefix also 404s.

### How it works

`~/.jupyter/jupyter_server_config.py` registers a named proxy server `fep` that
runs `python3.7 -m fep_web.app --host 127.0.0.1 --port 8080`. The proxy **starts
the app on demand** (so nothing to launch by hand — and nothing to leak: the app
binds loopback only, making the authenticated Jupyter proxy the sole way in).
Auth is JupyterHub's own; the app runs with **no `--token`**, i.e. there is no
second auth layer to configure.

Two constraints, both already encoded in that file:

* there is **no `cwd` option** (jupyter_server_proxy `config.py:115`
  `make_server_process`) — the repo is made importable via `PYTHONPATH`;
* `command` must be a **list**, and `{port}` is substituted by the proxy.

**Do not also run `serve_ui.py` / `python3 -m fep_web.app` by hand on 8080.**
That squats the port the proxy wants; the proxy's spawn then fails and `/fep/`
returns **500** (with a generic Jupyter error page that does not mention the
port). If `/fep/` 500s, check for a stray process on 8080 *first*.

### ⚠️ `~/.jupyter/` does NOT survive a container recycle

**This already happened.** The 2026-09-13 11:05 recycle replaced `~/.jupyter/`
with the platform default and removed `jupyter_server_config.py`, silently
taking the deployment down: `/fep/` stopped resolving (it now 302s to the
JupyterHub login, and the route behind it is gone). `jupyter_server_proxy`
itself survived — only the *registration* was lost. Nothing about the working
setup was recoverable from git.

So the config now lives **in the repo** and is restored by one command:

```bash
bash deploy/install_proxy_config.sh           # install / restore
bash deploy/install_proxy_config.sh --check   # verify, change nothing
```

Run it after every recycle. It validates the interpreter, the extension, and
that `fep_web.app` imports *under the proxy's own interpreter* (not the shell's
`python3`), then writes `~/.jupyter/jupyter_server_config.py` and re-parses the
result — a config that does not parse would take the whole Jupyter server down
on restart, so it is checked before you get there.

The file it installs is **inert until the Jupyter server restarts** (handlers
are registered at startup), so installing while a session is live is safe — the
route simply appears after the next WebIDE restart.

> Verify a deploy **after** the restart, not before: an unauthenticated `curl`
> to `/fep/` returns **302** either way, because JupyterHub's login redirect
> runs before route resolution. A 302 is therefore *not* evidence the route
> exists. Check the Jupyter log for a `404 GET .../fep/` after a real browser
> request instead.

### What was ruled out (keep, so it is not re-derived)

| outbound | result |
|---|---|
| TCP **443** | ✅ works (GitHub, Cloudflare API, gh-proxy) |
| TCP **80** | ✅ works |
| TCP 22 (SSH) | ❌ blocked — no SSH tunnels (localhost.run, serveo, pinggy) |
| UDP 7844 | ❌ `operation not permitted` |
| TCP 7844 | ❌ `i/o timeout` |

**cloudflared cannot work here**: it uses port **7844 for both QUIC and
HTTP/2**, so `--protocol http2` does not rescue it. The quick tunnel creates a
hostname (`*.trycloudflare.com`) and then never connects. ngrok would work
(dials out on 443) but needs an account authtoken — unnecessary now.

### What AI Studio supports

`jupyter_server_proxy` **is** installed in the WebIDE env — `3.2.4`, at
`/home/aistudio/external-libraries/lib/python3.7/site-packages/`, and
**enabled** (`external-libraries/etc/jupyter`). An earlier revision of this
section claimed it was absent; that was wrong. (It is *not* present in the
py3.10 `python35-paddle120-env`, so use the 3.7 interpreter for anything the
proxy runs.)

AI Studio's Codelab ships **`codelab_gradio_extension`** and
**`codelab_streamlit_extension`** (confirmed in `~/.codelab-jupyter.log`) — the
deploy button supports **Gradio and Streamlit only, not Flask**. The Gradio
extension discovers files ending in **`.gradio.py`**
(`codelab_gradio_extension/handlers.py`, `util.py`).

Its deploy handler takes a `file` argument resolved against
`jupyter_root_dir` (= `/home/aistudio`) — read from
`codelab_gradio_extension/handlers.py:check_file`.

⚠️ The same extension contains a **`DeployNotAvailableHandler`**, so deployment
is gated by project type. Baidu's docs say it requires **BML Codelab with
PaddlePaddle 2.5.2+**. Whether this project qualifies is **unconfirmed** — it is
the open question.

### What exists now

| file | purpose |
|---|---|
| **`deploy/install_proxy_config.sh`** | **restore the `/fep/` route after a recycle — run this first** |
| **`deploy/jupyter_server_config.py`** | the named-server registration, version-controlled so a recycle cannot lose it |
| `fep_ui.gradio.py` (repo root) | the Gradio app the deploy button wants (unconfirmed, unneeded) |
| `/home/aistudio/fep_ui.gradio.py` | symlink to the above, so it is visible from the project root |
| `serve_ui.py` (repo root) | Flask launcher, local use only — **not** deployable |
| `fep_web/app.py` | the Flask UI (prefix-tolerant) |

`fep_ui.gradio.py` uses only the API common to Gradio 3.x and 5.x and is
verified launching under **both** (3.19.1 on py3.7, 5.27.1 on py3.10), because
which interpreter the deploy picks is not knowable in advance.

### If `/fep/` ever breaks

Check these in order — the first two are 500s, the last two are 404s:

1. **Something else is on 8080.** `ps -eo pid,cmd | grep -E 'serve_ui|fep_web.app'`
   — a hand-started server makes the proxy's spawn fail with a bare 500.
2. **The proxy's child died.** The app also runs as a *child of the Jupyter
   server*, so the Jupyter log is where its traceback is, not `~/fep_ui.log`.
3. **Missing base prefix** → 404. You want `/<base>/fep/`, not `/fep/` and not
   `/proxy/fep/`.
4. **Proxy config not loaded** → 404 on the correct path. Confirm with
   `jupyter serverextension list` (expect `jupyter_server_proxy 3.2.4 enabled`);
   if it is missing, the WebIDE needs a restart, since the extension list is
   read at server start.

### The 部署 button — still unconfirmed, but no longer needed

Kept for reference only. The deploy handler packages a folder, uploads it to
BOS and runs it as a *separate cloud application* — with no GPU, no `namd3` and
no `/home/aistudio`, so a deployed UI could not see or drive the job anyway.
That is why §10 uses jupyter-server-proxy instead. `codelab_gradio_extension`
is still gated by `DeployNotAvailableHandler` (project type), unconfirmed.

### Local fallback (on the box only, NOT for browser access)

```bash
cd /home/aistudio/work/NAMD-FEP
python3 -m fep_web.app --root . --port 8080 --host 0.0.0.0 --token <secret>
```

Useful for debugging the app directly, but it **conflicts with the proxy** (see
above) — stop it before expecting `/fep/` to work. Binding `0.0.0.0` also
exposes the port to the LAN, which is why the proxy config deliberately uses
`127.0.0.1`.
