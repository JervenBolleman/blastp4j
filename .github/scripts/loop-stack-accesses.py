#!/usr/bin/env python3
"""Counts stack accesses (register spills and reloads) in one compiled function: statically, and weighted by samples.

Reads HotSpot -XX:+PrintAssembly output (with hsdis) or objdump output, x86-64 (AT&T syntax) or AArch64. Stack
accesses are memory operands based on the stack pointer (and on x29, the AArch64 frame pointer). Instructions are
"in a loop" when their basic block is in a cycle of the control-flow graph. With --samples (one sampled instruction
address per line, from "perf script -F ip"), the samples that fall in the function are counted per instruction, which
weighs each instruction by how often it ran.

Usage:
  loop-stack-accesses.py --header
  loop-stack-accesses.py --hotspot FILE --method NAME --label TEXT [--samples FILE]
  loop-stack-accesses.py --objdump FILE --label TEXT [--samples FILE]
For HotSpot, the static counts are of the last non-OSR C2 compilation of the method; the samples are matched against
all its C2 compilations.
"""
import argparse
import collections
import re
import sys

INSN = re.compile(r"^\s*(?:0x)?([0-9a-f]+):\s+(\S+)\s*(.*)$")
TARGET = re.compile(r"(?:^|[\s,#])(?:0x)?([0-9a-f]+)$")


def hotspot_blocks(path, method):
    """The listings of the C2 compilations of method, from the verified entry point: [(is_osr, lines)]."""
    blocks, current, started = [], None, False
    for line in open(path, errors="replace"):
        if line.startswith("Compiled method ("):
            # "Compiled method (c2) 3031 1409       4       pkg.Class::method (710 bytes)"; an OSR one has a "%".
            fields = line.split()
            current = None
            if fields[2] == "(c2)" and any(f.endswith("::" + method) for f in fields):
                current = []
                blocks.append(("%" in fields, current))
            started = False
        elif current is not None:
            if line.startswith("[Verified Entry Point]"):
                started = True
            elif started and line.startswith(("[Exception Handler]", "[Stub Code]", "[/Disassembly]", "-----")):
                current = None
            elif started:
                current.append(line)
    if not blocks:
        sys.exit(f"no C2 compilation of {method} in {path}")
    return blocks


def objdump_lines(path):
    """The instructions of the first function in an objdump listing."""
    out, inside = [], False
    for line in open(path, errors="replace"):
        if re.match(r"^[0-9a-f]+ <.*>:$", line.strip()):
            if inside:
                break
            inside = True
        elif inside:
            if not line.strip():
                if out:
                    break
                continue
            out.append(line)
    if not out:
        sys.exit(f"no function in {path}")
    return out


def parse(lines):
    insns = []
    for line in lines:
        line = line.split(";")[0]  # HotSpot comments
        line = line.split("//")[0]  # AArch64 objdump comments
        line = re.sub(r"\s+#\s.*$", "", line)  # x86 objdump comments
        line = re.sub(r"\s*<[^>]*>", "", line)  # symbolic branch targets
        m = INSN.match(line)
        if m:
            insns.append((int(m.group(1), 16), m.group(2).lower(), m.group(3).strip()))
    return insns


def branch_kind(mnemonic, ops):
    """'jump' (unconditional, direct), 'cond' (conditional, direct), 'end' (no fall-through) or None."""
    if mnemonic.startswith(("ret", "hlt", "ud2", "brk", "udf")) or mnemonic in ("br",):
        return "end"
    if mnemonic.startswith("jmp"):
        return "end" if "*" in ops or "%" in ops else "jump"
    if mnemonic == "b":
        return "jump"
    if mnemonic.startswith("j") or mnemonic.startswith(("b.", "cbz", "cbnz", "tbz", "tbnz")):
        return "cond"
    return None


def target(ops):
    m = TARGET.search(ops)
    return int(m.group(1), 16) if m else None


def in_loops(insns):
    """The addresses of the instructions whose basic block lies on a cycle of the control-flow graph."""
    addrs = [a for a, _, _ in insns]
    index = {a: i for i, a in enumerate(addrs)}
    leaders = {0}
    for i, (a, m, ops) in enumerate(insns):
        kind = branch_kind(m, ops)
        if kind:
            if i + 1 < len(insns):
                leaders.add(i + 1)
            t = target(ops) if kind != "end" else None
            if t in index:
                leaders.add(index[t])
    starts = sorted(leaders)
    block_of = {}
    for b, s in enumerate(starts):
        end = starts[b + 1] if b + 1 < len(starts) else len(insns)
        for i in range(s, end):
            block_of[i] = b
    succ = collections.defaultdict(list)
    for b, s in enumerate(starts):
        last = (starts[b + 1] if b + 1 < len(starts) else len(insns)) - 1
        a, m, ops = insns[last]
        kind = branch_kind(m, ops)
        if kind in ("jump", "cond"):
            t = target(ops)
            if t in index:
                succ[b].append(block_of[index[t]])
        if kind in (None, "cond") and last + 1 < len(insns):
            succ[b].append(b + 1)
    # Tarjan's strongly connected components, iterative.
    n = len(starts)
    num, low, on, stack, cyclic, counter = [-1] * n, [0] * n, [False] * n, [], set(), [0]
    for root in range(n):
        if num[root] >= 0:
            continue
        work = [(root, 0)]
        while work:
            v, k = work.pop()
            if k == 0:
                num[v] = low[v] = counter[0]
                counter[0] += 1
                stack.append(v)
                on[v] = True
            if k < len(succ[v]):
                work.append((v, k + 1))
                w = succ[v][k]
                if num[w] < 0:
                    work.append((w, 0))
                elif on[w]:
                    low[v] = min(low[v], num[w])
                continue
            if low[v] == num[v]:
                comp = []
                while True:
                    w = stack.pop()
                    on[w] = False
                    comp.append(w)
                    if w == v:
                        break
                if len(comp) > 1 or v in succ[v]:
                    cyclic.update(comp)
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[v])
    return {addrs[i] for i in range(len(insns)) if block_of[i] in cyclic}


def stack_access(mnemonic, ops):
    """'load', 'store' or None."""
    if "(%rsp" in ops:  # x86-64, AT&T: the destination is the last operand
        if mnemonic.startswith("push"):
            return "store"
        if mnemonic.startswith("pop"):
            return "load"
        last = ops.rsplit(",", 1)[-1]
        return "store" if "(%rsp" in last and mnemonic.lstrip("v").startswith("mov") else "load"
    if re.search(r"\[(sp|x29)\b", ops):  # AArch64
        if mnemonic.startswith("st"):
            return "store"
        if mnemonic.startswith("ld"):
            return "load"
    return None


def static_counts(insns):
    loop = in_loops(insns)
    c = collections.Counter()
    for a, m, ops in insns:
        kind = stack_access(m, ops)
        c["insns"] += 1
        c[kind] += kind is not None
        if a in loop:
            c["loop insns"] += 1
            c["loop " + str(kind)] += kind is not None
    return c


def sample_counts(listings, path):
    """Samples in the given instruction listings, and those on stack loads and stores."""
    kinds = {}
    for insns in listings:
        for a, m, ops in insns:
            kinds[a] = stack_access(m, ops) or "other"
    c = collections.Counter()
    for line in open(path, errors="replace"):
        line = line.strip()
        if line:
            try:
                kind = kinds.get(int(line, 16))
            except ValueError:
                continue
            if kind:
                c["total"] += 1
                c[kind] += 1
    return c


def pct(part, whole):
    return f"{100 * part / whole:.1f}%" if whole else "-"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hotspot")
    p.add_argument("--objdump")
    p.add_argument("--method")
    p.add_argument("--samples")
    p.add_argument("--label", default="")
    p.add_argument("--header", action="store_true")
    a = p.parse_args()
    if a.header:
        print("| code | instructions | stack loads / stores | in loops: instructions | stack loads / stores "
              "| samples | on stack loads | on stack stores |")
        print("|---|--:|--:|--:|--:|--:|--:|--:|")
        return
    if a.hotspot:
        blocks = hotspot_blocks(a.hotspot, a.method)
        normal = [lines for osr, lines in blocks if not osr]
        insns = parse(normal[-1] if normal else blocks[-1][1])
        listings = [parse(lines) for _, lines in blocks]
    else:
        insns = parse(objdump_lines(a.objdump))
        listings = [insns]
    c = static_counts(insns)
    row = (f"| {a.label} | {c['insns']} | {c['load']} / {c['store']} | {c['loop insns']} "
           f"| {c['loop load']} / {c['loop store']} |")
    if a.samples:
        s = sample_counts(listings, a.samples)
        row += f" {s['total']} | {pct(s['load'], s['total'])} | {pct(s['store'], s['total'])} |"
    else:
        row += " - | - | - |"
    print(row)


if __name__ == "__main__":
    main()
