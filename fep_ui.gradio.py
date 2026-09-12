#!/usr/bin/env python3
"""NAMD RBFE console — AI Studio deployment app.

The filename matters: AI Studio's Codelab deploy button
(`codelab_gradio_extension`) discovers files ending in ``.gradio.py``. That is
why this is Gradio and not the Flask app in `fep_web/app.py` — the platform's
deploy path supports Gradio and Streamlit only.

All the logic lives in the host-agnostic `fep_web` backend, so this file only
renders and routes. It deliberately uses the small API subset common to Gradio
3.x and 5.x (no `gr.Timer`, no `every=`, no Dataframe-schema assumptions):
the same file may run under either, depending on which interpreter the
extension picks.

Deploy: open this file in the AI Studio editor and click 部署.
"""
from __future__ import annotations

import sys
from pathlib import Path

import gradio as gr

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from fep_web import jobs, preflight, results          # noqa: E402
from fep_web.app import discover_systems              # noqa: E402
from fep_web.state import read_state                  # noqa: E402

PORT = 8080


# ---------------------------------------------------------------- rendering
def _systems() -> dict:
    return discover_systems(ROOT)


def _gpu_md() -> str:
    g = preflight()
    gate = "**Ready**" if g.ok else "**Blocked**"
    lines = [f"### GPU — {gate}", "", g.reason, ""]
    if g.devices:
        lines += ["| device | name | free / total | cc |",
                  "|---|---|---|---|"]
        for d in g.devices:
            lines.append(f"| {d.index} | {d.name} | "
                         f"{d.free_mib:.0f} / {d.total_mib:.0f} MiB | "
                         f"{d.compute_cap} |")
    else:
        lines.append("_No CUDA device visible — GPU-resident NAMD cannot run here._")
    return "\n".join(lines)


def _status_md(system: str) -> str:
    if not system:
        return "_No system selected._"
    st = read_state(_systems()[system])
    head = ("**running**" if st.running else
            "**complete**" if st.complete else "**idle**")
    out = [f"### {system} — {head}", ""]
    if st.controller_pid:
        out.append(f"controller pid `{st.controller_pid}`"
                   + (f", namd3 `{', '.join(map(str, st.namd_pids))}`"
                      if st.namd_pids else ""))
        out.append("")
        if st.active_window:
            out.append(f"active window: `{st.active_window}`")
            out.append("")
    out += ["| leg | stage | done | state | note |",
            "|---|---|---|---|---|"]
    for leg, ls in st.legs.items():
        for name, s in ls.stages.items():
            state = ("complete" if s.complete else
                     "running" if s.running else
                     f"partial" if s.done else "not started")
            out.append(f"| {leg} | {name} | {s.done}/{s.total} | {state} | "
                       f"{s.note or ''} |")
    done = sum(s.done for ls in st.legs.values() for s in ls.stages.values())
    total = sum(s.total for ls in st.legs.values() for s in ls.stages.values())
    out += ["", f"**overall {done}/{total} stages**"]
    return "\n".join(out)


def _choices():
    return sorted(_systems())


# ---------------------------------------------------------------- actions
def refresh(system):
    return _gpu_md(), _status_md(system), "\n".join(jobs.logs(ROOT / system, tail=80)) \
        if system else ""


def start(system):
    if not system:
        return "select a system first", _status_md(system)
    try:
        res = jobs.submit(_systems()[system])
        return f"✅ {res.reason}", _status_md(system)
    except jobs.JobError as exc:
        return f"⛔ {exc}", _status_md(system)


def stop(system):
    if not system:
        return "select a system first", _status_md(system)
    try:
        return f"🛑 {jobs.cancel(_systems()[system])}", _status_md(system)
    except jobs.JobError as exc:
        return f"⛔ {exc}", _status_md(system)


def analyse(system):
    if not system:
        return "select a system first"
    try:
        d = results.bar_result(_systems()[system])
    except results.AnalysisError as exc:
        return f"⛔ {exc}"
    return (f"### ΔΔG = {d['ddG']:+.3f} kcal/mol\n\n"
            f"| leg | 0 → 1 |\n|---|---|\n"
            f"| complex | {d['dG_complex']:+.3f} kcal/mol |\n"
            f"| solvent | {d['dG_solvent']:+.3f} kcal/mol |\n\n"
            f"_Per-window equilibration samples are discarded automatically. "
            f"Run `audit_fep.py` for error bars and convergence checks._")


# ---------------------------------------------------------------- layout
def build() -> gr.Blocks:
    names = _choices()
    with gr.Blocks(title="NAMD RBFE") as demo:
        gr.Markdown("# NAMD RBFE console\n"
                    "Relative binding free energy — job control and results.")

        with gr.Row():
            system = gr.Dropdown(choices=names,
                                 value=names[0] if names else None,
                                 label="System", scale=3)
            refresh_btn = gr.Button("Refresh", scale=1)

        gpu_md = gr.Markdown(_gpu_md())

        with gr.Row():
            start_btn = gr.Button("▶  Start / resume", variant="primary")
            stop_btn = gr.Button("■  Cancel")

        message = gr.Markdown("")
        status_md = gr.Markdown(_status_md(names[0]) if names else "")

        with gr.Row():
            analyse_btn = gr.Button("Compute ΔΔG")
        result_md = gr.Markdown("")

        gr.Markdown("### Controller log")
        log_box = gr.Textbox(lines=18, max_lines=18, show_label=False,
                             interactive=False,
                             placeholder="(no output yet)")

        refresh_btn.click(refresh, [system], [gpu_md, status_md, log_box])
        start_btn.click(start, [system], [message, status_md])
        stop_btn.click(stop, [system], [message, status_md])
        analyse_btn.click(analyse, [system], [result_md])
        system.change(refresh, [system], [gpu_md, status_md, log_box])

    return demo


if __name__ == "__main__":
    app = build()
    # server_name/port are what the AI Studio proxy expects.
    app.launch(server_name="0.0.0.0", server_port=PORT)
