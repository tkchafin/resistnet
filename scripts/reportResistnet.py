#!/usr/bin/env python3
"""Builds (or rebuilds) the interactive report of a ResistNet run."""
import argparse

from resistnet.report import build_report


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Builds the interactive HTML report of a runResistnet.py run "
            "(DIR/NAME/NAME.report.html) from its output files"
        )
    )
    parser.add_argument(
        "-o", "--out", required=True,
        help="The run's output directory and prefix, as given to "
             "runResistnet.py -o (DIR/NAME)"
    )
    args = parser.parse_args()
    build_report(args.out)


if __name__ == "__main__":
    main()
