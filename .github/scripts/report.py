#!/usr/bin/env python3
"""Builds the HTML report of an architecture comparison run (.github/workflows/arch-compare.yml).

Usage: report.py OUT_DIR JOB_DIR...

Each JOB_DIR is one job's report directory:
  machine.txt      key: value lines (architecture, CPU model, cores, memory, runner image, JDK, BLAST+, commits)
  report.md        the job's Markdown sections, as written to the step summary
  full-search.txt  the output of bin/compare-blastp.sh
  jvm-profile.html the async-profiler flame graph (optional)

Writes OUT_DIR/index.html, with a summary table of the full searches on top, and copies each flame graph to
OUT_DIR/<architecture>/jvm-profile.html. The Markdown is the small subset the workflow writes (headings, paragraphs,
fenced code, pipe tables, inline code), so no Markdown library is needed.
"""
import html
import os
import re
import shutil
import sys

RUN_URL = os.environ.get("RUN_URL", "")


def read(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            return f.read()
    except FileNotFoundError:
        return ""


def machine(text):
    m = {}
    for line in text.splitlines():
        key, sep, value = line.partition(":")
        if sep:
            m[key.strip()] = value.strip()
    return m


def inline(text):
    parts = re.split(r"(`[^`]*`)", text)
    return "".join(f"<code>{html.escape(p[1:-1])}</code>" if p.startswith("`") and p.endswith("`") and len(p) > 1
                   else html.escape(p) for p in parts)


def cells(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [c.strip() for c in line.split("|")]


def markdown(text, shift=1):
    """The workflow's Markdown subset as HTML; headings moved down by SHIFT levels."""
    out, lines, i = [], text.splitlines(), 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("```"):
            j = i + 1
            while j < len(lines) and not lines[j].startswith("```"):
                j += 1
            out.append("<pre>" + html.escape("\n".join(lines[i + 1:j])) + "</pre>")
            i = j + 1
        elif m := re.match(r"(#{1,6})\s+(.*)", line):
            level = min(6, len(m.group(1)) + shift)
            out.append(f"<h{level}>{inline(m.group(2))}</h{level}>")
            i += 1
        elif line.lstrip().startswith("|") and i + 1 < len(lines) and re.match(r"^\s*\|?[\s:|-]+\|?\s*$", lines[i + 1]):
            head = cells(line)
            right = [c.endswith(":") for c in cells(lines[i + 1])]
            j = i + 2
            rows = []
            while j < len(lines) and lines[j].lstrip().startswith("|"):
                rows.append(cells(lines[j]))
                j += 1

            def row(tag, cs):
                return "<tr>" + "".join(
                    f'<{tag}{" class=num" if k < len(right) and right[k] else ""}>{inline(c)}</{tag}>'
                    for k, c in enumerate(cs)) + "</tr>"

            out.append('<div class="scroll"><table><thead>' + row("th", head) + "</thead><tbody>"
                       + "".join(row("td", r) for r in rows) + "</tbody></table></div>")
            i = j
        elif line.strip():
            j = i + 1
            while (j < len(lines) and lines[j].strip() and not lines[j].startswith("```")
                   and not lines[j].startswith("#") and not lines[j].lstrip().startswith("|")):
                j += 1
            out.append("<p>" + inline(" ".join(lines[i:j])) + "</p>")
            i = j
        else:
            i += 1
    return "\n".join(out)


def full_search(text):
    """Median wall times and counters per system from bin/compare-blastp.sh's output."""
    wall, counters = {}, {}
    for line in text.splitlines():
        if m := re.match(r"(\w+) median wall: ([\d.]+) s", line):
            wall[m.group(1)] = float(m.group(2))
        elif m := re.match(r"(\w+): (instructions:u=.*)", line):
            counters[m.group(1)] = {k: float(v) for k, v in re.findall(r"([\w-]+):u=([\d.]+)G", m.group(2))}
    identical = "DIFFERENT" not in text and "(identical)" in text
    return wall, counters, identical


def ratio(a, b):
    return f"{a / b:.2f}×" if a and b else "–"


def main():
    out, jobs = sys.argv[1], sys.argv[2:]
    os.makedirs(out, exist_ok=True)
    summary, sections, commits = [], [], set()
    for job in sorted(jobs):
        m = machine(read(os.path.join(job, "machine.txt")))
        arch = m.get("Architecture", os.path.basename(job))
        wall, counters, identical = full_search(read(os.path.join(job, "full-search.txt")))
        blast, jvm = counters.get("blast", {}), counters.get("jvm", {})
        commits.add(m.get("Code", ""))
        profile = ""
        if os.path.exists(os.path.join(job, "jvm-profile.html")):
            os.makedirs(os.path.join(out, arch), exist_ok=True)
            shutil.copy(os.path.join(job, "jvm-profile.html"), os.path.join(out, arch, "jvm-profile.html"))
            profile = f'<a href="{html.escape(arch)}/jvm-profile.html">flame graph</a>'
        summary.append(
            "<tr>"
            f"<td>{html.escape(arch)}</td><td>{html.escape(m.get('CPU model', '?'))}</td>"
            f"<td class=num>{html.escape(m.get('CPUs', '?'))}</td>"
            f"<td class=num>{wall.get('blast', '–')}</td><td class=num>{wall.get('jvm', '–')}</td>"
            f"<td class=num><b>{ratio(wall.get('jvm'), wall.get('blast'))}</b></td>"
            f"<td class=num>{ratio(jvm.get('cycles'), blast.get('cycles'))}</td>"
            f"<td class=num>{ratio(jvm.get('instructions'), blast.get('instructions'))}</td>"
            + (f"<td class=num>{wall.get('ref', '–')}</td>" if "ref" in wall else "")
            + f"<td>{'identical' if identical else '<b class=bad>different</b>'}</td><td>{profile}</td></tr>")
        facts = "".join(f"<dt>{html.escape(k)}</dt><dd>{inline(v)}</dd>" for k, v in m.items())
        sections.append(
            f'<section><h2 id="{html.escape(arch)}">{html.escape(arch)}: {html.escape(m.get("CPU model", ""))}</h2>'
            f'<dl class="facts">{facts}</dl>'
            + markdown(read(os.path.join(job, "report.md")), shift=1)
            + "</section>")
    has_ref = any("ref median wall" in read(os.path.join(j, "full-search.txt")) for j in jobs)
    commit = ", ".join(sorted(c for c in commits if c)) or "?"
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>blastp4j vs BLAST+</title>
<style>
:root {{ --bg: #fdfdfc; --fg: #1d1d1b; --muted: #66645e; --line: #dddad3; --head: #f1efe9; --accent: #2f5d8a; --bad: #b3261e; }}
@media (prefers-color-scheme: dark) {{
  :root {{ --bg: #17181a; --fg: #e6e4df; --muted: #a09d96; --line: #34363a; --head: #222428; --accent: #8db7e0; --bad: #f2867d; }}
}}
body {{ background: var(--bg); color: var(--fg); font: 15px/1.5 system-ui, sans-serif; margin: 0 auto; max-width: 1100px; padding: 24px 16px 64px; }}
h1 {{ font-size: 1.6rem; margin: 0 0 4px; }}
h2 {{ font-size: 1.3rem; margin-top: 2.5rem; border-bottom: 1px solid var(--line); padding-bottom: 4px; }}
h3, h4 {{ font-size: 1.05rem; margin-top: 1.6rem; }}
a {{ color: var(--accent); }}
.muted {{ color: var(--muted); }}
.bad {{ color: var(--bad); }}
.scroll {{ overflow-x: auto; }}
table {{ border-collapse: collapse; font-size: 13px; margin: 8px 0; }}
th, td {{ border: 1px solid var(--line); padding: 4px 8px; text-align: left; vertical-align: top; }}
th {{ background: var(--head); }}
.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
pre {{ background: var(--head); border: 1px solid var(--line); padding: 8px 10px; overflow-x: auto; font-size: 12px; }}
code {{ font-size: 0.9em; }}
dl.facts {{ display: grid; grid-template-columns: max-content 1fr; gap: 2px 16px; font-size: 13px; }}
dl.facts dt {{ color: var(--muted); }}
dl.facts dd {{ margin: 0; overflow-wrap: anywhere; }}
</style>
</head>
<body>
<h1>blastp4j against BLAST+</h1>
<p class="muted">Code {inline(commit)}{f' · <a href="{html.escape(RUN_URL)}">workflow run</a>' if RUN_URL else ''}.
Hosted runners are virtual machines and differ from run to run: compare the ratios within one row, not times across runs.</p>
<div class="scroll"><table>
<thead><tr><th>architecture</th><th>CPU model</th><th class=num>CPUs</th><th class=num>BLAST+ wall s</th>
<th class=num>blastp4j wall s</th><th class=num>wall</th><th class=num>cycles</th><th class=num>instructions</th>
{'<th class=num>reference wall s</th>' if has_ref else ''}<th>output</th><th></th></tr></thead>
<tbody>{''.join(summary)}</tbody>
</table></div>
<p class="muted">Wall times are medians of the full-search rounds; the ratios are blastp4j over BLAST+.</p>
{''.join(sections)}
</body>
</html>
"""
    with open(os.path.join(out, "index.html"), "w", encoding="utf-8") as f:
        f.write(page)


if __name__ == "__main__":
    main()
