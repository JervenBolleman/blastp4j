#!/usr/bin/env python3
"""Prints JMH CSV results side by side as a Markdown table: one row per benchmark and parameter set, one column per
run.

Usage: jmh-table.py LABEL=FILE.csv...

Each cell is the score and its 99.9% error. Benchmarks missing from a run (it failed, or did not exist in that code)
have an empty cell.
"""
import csv
import sys


def main(args):
    runs = [a.split("=", 1) for a in args]
    rows = {}  # (benchmark, params) -> unit, in order of first appearance
    cells = {}  # (label, benchmark, params) -> text
    for label, path in runs:
        try:
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    name = r["Benchmark"].rsplit(".", 2)
                    name = ".".join(name[-2:])
                    params = ", ".join(f"{k[len('Param: '):]}={v}" for k, v in r.items()
                                       if k.startswith("Param: ") and v)
                    key = (name, params)
                    rows.setdefault(key, r["Unit"])
                    error = r["Score Error (99.9%)"]
                    cells[(label, *key)] = fmt(r["Score"]) + (" ± " + fmt(error) if error not in ("", "NaN") else "")
        except FileNotFoundError:
            pass
    labels = [label for label, _ in runs]
    print("| benchmark | parameters | unit | " + " | ".join(labels) + " |")
    print("|---|---|---|" + "--:|" * len(labels))
    for (name, params), unit in rows.items():
        print(f"| {name} | {params} | {unit} | " + " | ".join(cells.get((l, name, params), "") for l in labels) + " |")


def fmt(s):
    x = float(s)
    return f"{x:.1f}" if abs(x) >= 10 else f"{x:.3g}"


if __name__ == "__main__":
    main(sys.argv[1:])
