#!/usr/bin/env python3
"""Top-level CLI for the scDaisychain XCI package."""

from __future__ import annotations

import argparse
import sys

from .pipeline import add_run_arguments, run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="scDaisychain",
        description="Run or inspect the scDaisychain/XCI BAM-to-matrix pipeline.",
    )
    sub = parser.add_subparsers(dest="command")

    run = sub.add_parser("run", help="Run the full pipeline.")
    add_run_arguments(run)

    return parser


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()

    if not argv:
        parser.print_help()
        return 0

    args = parser.parse_args(argv)

    if args.command == "run":
        return run_pipeline(args)

    parser.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
