#!/usr/bin/env python3
"""Measures the noise of the architecture comparison from the runs it has already made.

Usage: run-noise.py [REF] [MIN_REPEATS]

REF (default origin/reports) is the reports branch, whose runs/<id>/summary.json this reads; the branch must be
fetched. MIN_REPEATS (default 3) is how many runs a commit needs before its repeats are reported.

Three measurements, each per architecture:

1. The same code on both sides of a job. A job runs the code and its reference in the same interleaved rounds, so
   where the two sides hold the same program their difference is pure noise. Three sources of such pairs: a
   calibration run (`code` and `reference` the same commit), a commit whose diff touches no file the programs are
   built from (the workflow, the report script, the comparison script), and, for fastblastp alone, a commit that
   touches only the sources blastp uses and blastp-fast does not (BatchWordFinder, WordFinder). The spread of those
   differences is what one paired run can resolve. It is printed for three ways of summarising a side's rounds
   (their median, their mean, their minimum), because that choice changes the answer.

2. Where that noise lives: the differences of the single rounds within a job, against the differences of the job
   means. If a job's mean is no tighter than its rounds allow (round-to-round / sqrt(rounds)), the noise is
   round-to-round and more rounds buy their square root; what is left over attaches to the side for the whole job
   and more rounds cannot touch it.

3. One commit measured in several runs. Hosted runners differ from run to run and the x86-64 ones differ in CPU
   model, so this is grouped by model: the spread of BLAST+'s wall time, of blastp4j's, and of their ratio.
"""

import json
import math
import statistics as st
import subprocess
import sys
from collections import defaultdict

# The sources blastp uses and blastp-fast does not, so that a commit touching only these leaves fastblastp's
# program identical on both sides of a job (blastp-fast seeds with FastBatchWordFinder and CompressedTable).
BLASTP_ONLY = {"src/main/java/swiss/sib/blastp4j/seed/BatchWordFinder.java",
               "src/main/java/swiss/sib/blastp4j/seed/WordFinder.java"}
# Paths that change what is measured without changing either program.
HARNESS = (".github/", "bin/compare-blastp.sh", "bin/events.py", "bin/machcode.py", "README.md", "docs/",
           "reference/", "for_ncbi")


def summaries(ref):
    """Every run's summary.json on the reports branch, oldest first, as rows with the run's fields folded in."""
    names = subprocess.run(["git", "ls-tree", "--name-only", ref, "runs/"],
                           capture_output=True, text=True, check=True).stdout.split()
    out = []
    for name in sorted(names):
        path = f"{ref}:{name}/summary.json"
        got = subprocess.run(["git", "show", path], capture_output=True, text=True)
        if got.returncode:
            continue
        run = json.loads(got.stdout)
        for row in run.get("rows", []):
            row = dict(row)
            row["run_id"] = run.get("run_id", name.split("/")[-1])
            row["date"] = run.get("date", "")
            out.append(row)
    return out


def changed(a, b):
    """The files a..b changes, or None when either commit is not in this clone."""
    got = subprocess.run(["git", "diff", "--name-only", a, b], capture_output=True, text=True)
    return got.stdout.split() if got.returncode == 0 else None


def same_code(a, b):
    """Whether a..b leaves (the whole program, fastblastp's program) identical."""
    files = changed(a, b)
    if files is None:
        return False, False
    build = [f for f in files if not f.startswith(HARNESS)]
    return not build, all(f in BLASTP_ONLY for f in build)


def rounds(row, key):
    return [x for x in (row.get("rounds") or {}).get(key) or [] if x]


def pairs(rows, which):
    """The jobs whose two sides of WHICH ("full" or "fast") hold the same program, as (code rounds, ref rounds)."""
    seen, out = {}, []
    for row in rows:
        a, b = row.get("ref_commit"), row.get("jvm_commit")
        if not a or not b:
            continue
        if a == b:
            whole = fast = True
        else:
            if (a, b) not in seen:
                seen[(a, b)] = same_code(a, b)
            whole, fast = seen[(a, b)]
        if not (whole if which == "full" else whole or fast):
            continue
        code = rounds(row, "jvm" if which == "full" else "fast")
        ref = rounds(row, "ref" if which == "full" else "fastref")
        if len(code) >= 2 and len(ref) >= 2 and len(code) == len(ref):
            out.append((row, code, ref))
    return out


def rms(v):
    return math.sqrt(st.mean([x * x for x in v])) if v else float("nan")


def report(rows, repeats):
    arches = sorted({r.get("architecture") for r in rows if r.get("architecture") in ("aarch64", "x86_64")})

    print("1. The same code on both sides of a job: the difference is noise.")
    print(f"   {'architecture':12} {'what':5} {'jobs':>4} {'mean':>8} {'sd of medians':>14} "
          f"{'sd of means':>12} {'sd of minima':>13}")
    for arch in arches:
        for which in ("full", "fast"):
            got = [(c, f) for row, c, f in pairs(rows, which) if row.get("architecture") == arch]
            if len(got) < 3:
                print(f"   {arch:12} {which:5} {len(got):4d}  too few")
                continue
            def by(f):
                return [(f(c) - f(r)) / f(r) * 100 for c, r in got]
            print(f"   {arch:12} {which:5} {len(got):4d} {st.mean(by(st.mean)):+7.2f}% "
                  f"{st.stdev(by(st.median)):13.2f}% {st.stdev(by(st.mean)):11.2f}% {st.stdev(by(min)):12.2f}%")

    print("\n2. Where the noise lives.")
    print(f"   {'architecture':12} {'what':5} {'jobs':>4} {'round to round':>15} {'job means':>10} "
          f"{'rounds alone':>13} {'job level':>10}")
    for arch in arches:
        for which in ("full", "fast"):
            got = [(c, f) for row, c, f in pairs(rows, which) if row.get("architecture") == arch]
            if len(got) < 3:
                continue
            within, means, k = [], [], []
            for c, r in got:
                d = [(x - y) / y * 100 for x, y in zip(c, r)]
                within.append(st.stdev(d))
                means.append(st.mean(d))
                k.append(len(d))
            per_round, per_job = rms(within), st.stdev(means)
            alone = per_round / math.sqrt(st.mean(k))
            left = math.sqrt(max(0.0, per_job ** 2 - alone ** 2))
            print(f"   {arch:12} {which:5} {len(got):4d} {per_round:14.2f}% {per_job:9.2f}% "
                  f"{alone:12.2f}% {left:9.2f}%")

    print("\n3. One commit measured in several runs, by CPU model.")
    groups = defaultdict(list)
    for row in rows:
        if row.get("jvm_commit") and row.get("rounds"):
            groups[(row["jvm_commit"], row.get("architecture"), row.get("cpu"))].append(row)
    print(f"   {'commit':8} {'architecture':12} {'cpu':34} {'runs':>4} {'BLAST+':>8} {'blastp4j':>9} {'ratio':>7}")
    for (commit, arch, cpu), got in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        if len(got) < repeats:
            continue
        def cv(key):
            v = [st.mean(rounds(r, key)) for r in got if rounds(r, key)]
            return st.stdev(v) / st.mean(v) * 100 if len(v) > 1 else float("nan")
        ratios = [st.mean(rounds(r, "jvm")) / st.mean(rounds(r, "blast"))
                  for r in got if rounds(r, "jvm") and rounds(r, "blast")]
        rcv = st.stdev(ratios) / st.mean(ratios) * 100 if len(ratios) > 1 else float("nan")
        print(f"   {commit:8} {arch or '?':12} {(cpu or '?')[:34]:34} {len(got):4d} "
              f"{cv('blast'):7.2f}% {cv('jvm'):8.2f}% {rcv:6.2f}%")
    print("   (each figure is a coefficient of variation over the runs: BLAST+'s wall time, blastp4j's, "
          "and blastp4j/BLAST+.)")


if __name__ == "__main__":
    ref = sys.argv[1] if len(sys.argv) > 1 else "origin/reports"
    least = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    report(summaries(ref), least)
