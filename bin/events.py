#!/usr/bin/env python3
"""Splits the hardware events of a run by compiled method, as machcode.py does for its cycles.

Usage: events.py LOG SAMPLES [STAT [TOP]]

LOG is the JVM's output with -XX:CompileCommand=print for some methods (see machcode.py). SAMPLES is
"perf script -F event,period,ip" of a perf record of the same run with several events, one sample per line:
its period, its event and the sampled instruction address (hex). Each sample counts for its period, so the
totals are the events' counts. Prints, for the TOP (default 15) compiled methods with the most cycles, each
event's share of its total and its count per 1,000 cycles in that method; then the 10 instructions with the
most of each event other than cycles.

With more events than the core has counters, perf multiplexes them: each event is counted only part of the
time, so the shares of one event are right but the counts of two events are not comparable. STAT, a
"perf stat -x," output of a run of the same program with the same events, then gives each event's count per
cycle, and the counts per 1,000 cycles are scaled to it.

Code freed by the JVM can be reused for another method at the same addresses; methods whose address ranges
overlap share one row, and the samples in it are counted once.
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


def stat(path):
    """The events' counts in a perf stat -x, output: {event: count}."""
    counts = {}
    with open(path) as f:
        for line in f:
            fields = line.strip().split(",")
            if len(fields) > 2 and fields[0].isdigit():
                counts[fields[2].split(":")[0]] = int(fields[0])
    return counts


def main(log, samples, stat_path=None, top=15):
    counts = read(samples)
    cycles = next((e for e in counts if e in CYCLES), None)
    if cycles is None:
        print("no cycles in the samples")
        return
    events = [cycles] + sorted(e for e in counts if e != cycles)
    totals = {e: sum(counts[e].values()) for e in events}
    # The factor of each event's count per cycle in the samples to that in perf stat.
    scale = {e: 1.0 for e in events}
    if stat_path:
        counted = stat(stat_path)
        for e in events[1:]:
            if counted.get(e) and counted.get(cycles) and totals[e]:
                scale[e] = (counted[e] / counted[cycles]) / (totals[e] / totals[cycles])
    # Overlapping address ranges (code reused after the JVM freed it) merged into one.
    ranges = []
    for m in sorted(methods(log), key=lambda m: m["start"]):
        last = ranges[-1] if ranges else None
        if last and m["start"] < last["end"]:
            names = last["name"].split(" / ")
            name = last["name"] if m["name"] in names else last["name"] + " / " + m["name"]
            ranges[-1] = dict(last, end=max(last["end"], m["end"]), name=name)
        else:
            ranges.append(m)
    rows = []
    for m in ranges:
        inside = {e: sum(n for a, n in counts[e].items() if m["start"] <= a < m["end"]) for e in events}
        if inside[cycles]:
            rows.append((m, inside))
    rows.sort(key=lambda r: -r[1][cycles])
    print("totals: " + ", ".join(f"{e} {totals[e] / 1e9:.1f}G" for e in events))
    print("share of each event's total (per 1,000 cycles in the method"
          + (", scaled to perf stat's counts per cycle):" if stat_path else "):"))
    print(f"{'':48s}" + "".join(f"{e[:24]:>26s}" for e in events))
    for m, inside in rows[:top]:
        cells = [f"{100 * inside[cycles] / totals[cycles]:25.2f}%"]
        for e in events[1:]:
            rate = 1000 * scale[e] * inside[e] / inside[cycles]
            cells.append(f"{100 * inside[e] / totals[e]:15.2f}% ({rate:6.1f})")
        print(f"{m['name'][-48:]:48s}" + "".join(cells))
    bounds = [(m["start"], m["end"]) for m, _ in rows]
    covered = {e: sum(n for a, n in counts[e].items() if any(s <= a < t for s, t in bounds)) for e in events}
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
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None,
         int(sys.argv[4]) if len(sys.argv) > 4 else 15)
