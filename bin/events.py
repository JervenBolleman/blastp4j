#!/usr/bin/env python3
"""Splits the hardware events of a run by compiled method, as machcode.py does for its cycles.

Usage: events.py LOG SAMPLES [TOP]

LOG is the JVM's output with -XX:CompileCommand=print for some methods (see machcode.py). SAMPLES is
"perf script -F event,period,ip" of a perf record of the same run with several events, one sample per line:
its period, its event and the sampled instruction address (hex). Each sample counts for its period, so the
totals are the events' counts. Prints, for the TOP (default 15) compiled methods with the most cycles, each
event's share of its total and its count per 1,000 cycles in that method; then the 10 instructions with the
most of each event other than cycles.
"""
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from machcode import methods  # noqa: E402

CYCLES = ("cycles", "cpu-clock")


def read(samples):
    """The events' counts by instruction address: {event: Counter(address: count)}."""
    counts = collections.defaultdict(collections.Counter)
    with open(samples) as f:
        for line in f:
            fields = line.split()
            event = next((x[:-1] for x in fields if x.endswith(":")), None)
            if event is None or len(fields) < 3:
                continue
            try:
                period, address = int(fields[0]), int(fields[-1], 16)
            except ValueError:
                continue
            counts[event.split(":")[0]][address] += period
    return counts


def main(log, samples, top=15):
    counts = read(samples)
    cycles = next((e for e in counts if e in CYCLES), None)
    if cycles is None:
        print("no cycles in the samples")
        return
    events = [cycles] + sorted(e for e in counts if e != cycles)
    totals = {e: sum(counts[e].values()) for e in events}
    rows = []
    for m in methods(log):
        inside = {e: sum(n for a, n in counts[e].items() if m["start"] <= a < m["end"]) for e in events}
        if inside[cycles]:
            rows.append((m, inside))
    rows.sort(key=lambda r: -r[1][cycles])
    print("totals: " + ", ".join(f"{e} {totals[e] / 1e9:.1f}G" for e in events))
    print("share of each event's total (per 1,000 cycles in the method):")
    print(f"{'':48s}" + "".join(f"{e[:24]:>26s}" for e in events))
    for m, inside in rows[:top]:
        cells = [f"{100 * inside[cycles] / totals[cycles]:25.2f}%"]
        for e in events[1:]:
            rate = 1000 * inside[e] / inside[cycles]
            cells.append(f"{100 * inside[e] / totals[e]:15.2f}% ({rate:6.1f})")
        print(f"{m['name'][-48:]:48s}" + "".join(cells))
    covered = {e: sum(r[1][e] for r in rows) for e in events}
    print(f"{'elsewhere':48s}" + "".join(f"{100 * (1 - covered[e] / totals[e]):25.2f}%" for e in events))
    names = {}
    for m, _ in rows:
        names[(m["start"], m["end"])] = m["name"]
    for e in events[1:]:
        print(f"\nthe 10 instructions with the most {e}:")
        for a, n in counts[e].most_common(10):
            name = next((v for (s, t), v in names.items() if s <= a < t), "elsewhere")
            print(f"  {100 * n / totals[e]:6.2f}%  {a:#x}  {name}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 15)
