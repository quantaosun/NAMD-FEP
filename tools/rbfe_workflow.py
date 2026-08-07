#!/usr/bin/env python3
"""Small, dependency-free helpers for running independent NAMD FEP windows."""

from __future__ import annotations

import argparse
import json
import shlex
from pathlib import Path


def lambda_windows(count: int) -> list[tuple[float, float]]:
    if count < 2:
        raise ValueError("at least two lambda windows are required")
    step = 1.0 / (count - 1)
    return [(round(i * step, 10), round((i + 1) * step, 10)) for i in range(count - 1)]


def validate_inputs(files: list[Path]) -> None:
    missing = [str(path) for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing required input: " + ", ".join(missing))


def render_window(template: str, lambda_value: float, next_lambda: float) -> str:
    replacements = {
        "@ALCH_LAMBDA@": f"{lambda_value:.10g}",
        "@ALCH_LAMBDA2@": f"{next_lambda:.10g}",
    }
    for marker, value in replacements.items():
        template = template.replace(marker, value)
    if "@" in template:
        raise ValueError("template contains an unresolved @...@ marker")
    return template


def namd_command(namd: Path, config: Path, cpus: int, device: int) -> list[str]:
    if cpus < 1:
        raise ValueError("cpus must be positive")
    return [
        str(namd),
        f"+p{cpus}",
        "+setcpuaffinity",
        "--CUDASOAintegrate",
        "on",
        "+devices",
        str(device),
        str(config),
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate")
    validate.add_argument("files", nargs="+", type=Path)

    generate = subparsers.add_parser("generate")
    generate.add_argument("--template", type=Path, required=True)
    generate.add_argument("--output-dir", type=Path, required=True)
    generate.add_argument("--windows", type=int, default=16)

    command = subparsers.add_parser("command")
    command.add_argument("--namd", type=Path, required=True)
    command.add_argument("--config", type=Path, required=True)
    command.add_argument("--cpus", type=int, default=4)
    command.add_argument("--device", type=int, default=0)

    args = parser.parse_args()
    if args.command == "validate":
        validate_inputs(args.files)
        print(json.dumps({"valid": True, "files": [str(p) for p in args.files]}))
    elif args.command == "generate":
        validate_inputs([args.template])
        template = args.template.read_text()
        args.output_dir.mkdir(parents=True, exist_ok=True)
        windows = lambda_windows(args.windows)
        for index, (lam, lam2) in enumerate(windows):
            output = args.output_dir / f"window_{index:03d}.namd"
            output.write_text(render_window(template, lam, lam2))
        print(f"generated {len(windows)} independent windows in {args.output_dir}")
    else:
        print(shlex.join(namd_command(args.namd, args.config, args.cpus, args.device)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
