#!/usr/bin/env python3
"""Annotates the machine code the JVM printed with CPU samples of the same process.

Usage: machcode.py LOG IPS [TOP]

LOG is the JVM's output with -XX:CompileCommand=print for some methods. With the hsdis library the JVM prints the
disassembly; without it, the code as bytes ([MachCode]), which is then disassembled here with objdump (for this
machine's architecture). IPS is "perf script -F ip" of a perf record of the same run, one sampled instruction address
(hex) per line. Each instruction is shown with its samples: first the TOP (default 40) instructions with the
most samples, then every compiled method that has samples, in full (after the first line starting with ==). Shares are of all samples of the run.
"""
import collections
import os
import platform
import re
import subprocess
import sys
import tempfile

HEADER = re.compile(r"^Compiled method \((c1|c2)\)\s+\d+\s+\d+\s+(%?)\s*\d\s+(\S+)")
INSN = re.compile(r"^\s+0x([0-9a-f]+):\s+([a-z][^;]*?)\s*(?:;.*)?$")
CODE = re.compile(r"^\s+0x([0-9a-f]+): ((?:[0-9a-f]{4} ?\|? ?)+)\s*$")
COMMENT = re.compile(r"^\s+0x([0-9a-f]+): +;(.*)")


def methods(log):
    """The compiled methods: name, compiler, their instructions (hsdis) or code bytes, and comments by address."""
    found = []
    current = None
    in_code = False
    with open(log, errors="replace") as f:
        for line in f:
            m = HEADER.match(line)
            if m:
                current = {"name": m.group(3) + (" (OSR)" if m.group(2) else ""), "compiler": m.group(1),
                           "code": bytearray(), "start": None, "insns": [],
                           "comments": collections.defaultdict(list)}
                found.append(current)
                in_code = in_disassembly = False
                continue
            if current is None:
                continue
            if line.startswith("[MachCode]"):
                in_code = True
            elif line.startswith("[/MachCode]"):
                in_code = False
            elif line.startswith("[Disassembly]"):
                in_disassembly = True
            elif line.startswith("[/Disassembly]"):
                in_disassembly = False
            elif in_disassembly:
                m = COMMENT.match(line)
                if m:
                    current["comments"][int(m.group(1), 16)].append(m.group(2).strip())
                    continue
                m = INSN.match(line)
                if m:
                    current["insns"].append((int(m.group(1), 16), re.sub(r"\s+", " ", m.group(2))))
            elif in_code:
                m = CODE.match(line)
                if m:
                    address = int(m.group(1), 16)
                    data = bytes.fromhex(m.group(2).replace("|", "").replace(" ", ""))
                    if current["start"] is None:
                        current["start"] = address
                    # Fill a gap (code not printed) with zeros, so that addresses stay right.
                    end = current["start"] + len(current["code"])
                    if address > end:
                        current["code"] += bytes(address - end)
                    if address >= end:
                        current["code"] += data
                    continue
                m = COMMENT.match(line)
                if m:
                    current["comments"][int(m.group(1), 16)].append(m.group(2).strip())
    for m in found:
        if m["insns"]:
            m["start"] = m["insns"][0][0]
            m["end"] = m["insns"][-1][0] + 16
        elif m["start"] is not None:
            m["end"] = m["start"] + len(m["code"])
    return [m for m in found if m["insns"] or m["code"]]


def disassemble(start, code):
    """objdump's disassembly of the code at address start: (address, text) per instruction."""
    machine = platform.machine()
    arch = ["-m", "aarch64"] if machine in ("aarch64", "arm64") else ["-m", "i386:x86-64", "-M", "intel"]
    with tempfile.NamedTemporaryFile(delete=False) as f:
        f.write(code)
    try:
        out = subprocess.run(["objdump", "-D", "-b", "binary", *arch, "--adjust-vma=" + hex(start), f.name],
                             capture_output=True, text=True, check=True).stdout
    finally:
        os.unlink(f.name)
    result = []
    for line in out.splitlines():
        m = re.match(r"^\s*([0-9a-f]+):\t(?:[0-9a-f]{2,8} ?)+\s*\t(.*)$", line)
        if m:
            result.append((int(m.group(1), 16), m.group(2).strip()))
    return result


def main(log, ips, top=40):
    samples = collections.Counter()
    with open(ips) as f:
        for line in f:
            line = line.strip()
            if line:
                samples[int(line, 16)] += 1
    total = sum(samples.values())
    if not total:
        print("no samples")
        return
    listings = []
    hot = []
    for m in methods(log):
        start, end = m["start"], m["end"]
        inside = sum(n for a, n in samples.items() if start <= a < end)
        if not inside:
            continue
        lines = m["insns"] or disassemble(start, bytes(m["code"]))
        listings.append((inside, m, lines))
        for a, text in lines:
            if samples[a]:
                hot.append((samples[a], a, m["name"], text))
    attributed = sum(inside for inside, _, _ in listings)
    print(f"{total} samples, {100 * attributed / total:.2f}% in printed compiled methods "
          f"({100 * (total - attributed) / total:.2f}% elsewhere: interpreter, stubs, native, unprinted code)")
    print("the 15 printed compiled methods with the most:")
    for inside, m, _ in sorted(listings, key=lambda x: -x[0])[:15]:
        print(f"  {100 * inside / total:6.2f}%  {m['compiler']} {m['name']}")
    print(f"\nthe {top} instructions with the most samples:")
    for n, a, name, text in sorted(hot, reverse=True)[:top]:
        print(f"  {100 * n / total:6.2f}%  {a:#x}  {text:48s} {name}")
    for inside, m, lines in sorted(listings, key=lambda x: -x[0]):
        print(f"\n== {m['compiler']} {m['name']}: {100 * inside / total:.2f}% of the samples")
        for a, text in lines:
            for c in m["comments"].get(a, []):
                print(f"{'':17s};{c}")
            n = samples[a]
            share = f"{100 * n / total:6.2f}%" if n else ""
            print(f"{share:>8s}  {a:#x}  {text}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 40)
