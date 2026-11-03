#!/usr/bin/env python3
"""Builds the HTML reports of the architecture comparison runs (.github/workflows/arch-compare.yml).

Usage: report.py run OUT_DIR JOB_DIR...
       report.py index SITE_DIR

"run" builds one run's report. Each JOB_DIR is one job's report directory:
  machine.txt      key: value lines (architecture, CPU model, cores, memory, runner image, JDK, BLAST+, commits)
  report.md        the job's Markdown sections, as written to the step summary
  full-search.txt  the output of bin/compare-blastp.sh
  jvm-profile.html the async-profiler flame graph (optional)
It writes OUT_DIR/index.html, with a summary table of the full searches on top, OUT_DIR/summary.json with that table's
numbers, and each flame graph as OUT_DIR/<architecture>/jvm-profile.html. The Markdown is the small subset the
workflow writes (headings, paragraphs, fenced code, pipe tables, inline code), so no Markdown library is needed.

"index" writes SITE_DIR/index.html: every run in SITE_DIR/runs/*/summary.json, newest first.

Settings (environment, for "run"): RUN_ID, RUN_URL, REPO_URL (for links to the commits), RUN_DATE (default now),
RUN_EVENT (what started the run: schedule or workflow_dispatch).
"""
import datetime
import json
import html
import os
import re
import shutil
import sys

RUN_URL = os.environ.get("RUN_URL", "")
REPO_URL = os.environ.get("REPO_URL", "")


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
    return a / b if a and b else None


def times(r):
    return f"{r:.2f}×" if r else "–"


STYLE = """
:root { --bg: #fdfdfc; --fg: #1d1d1b; --muted: #66645e; --line: #dddad3; --head: #f1efe9; --accent: #2f5d8a; --bad: #b3261e; }
@media (prefers-color-scheme: dark) {
  :root { --bg: #17181a; --fg: #e6e4df; --muted: #a09d96; --line: #34363a; --head: #222428; --accent: #8db7e0; --bad: #f2867d; }
}
body { background: var(--bg); color: var(--fg); font: 15px/1.5 system-ui, sans-serif; margin: 0 auto; max-width: 1100px; padding: 24px 16px 64px; }
h1 { font-size: 1.6rem; margin: 0 0 4px; }
h2 { font-size: 1.3rem; margin-top: 2.5rem; border-bottom: 1px solid var(--line); padding-bottom: 4px; }
h3, h4 { font-size: 1.05rem; margin-top: 1.6rem; }
a { color: var(--accent); }
.muted { color: var(--muted); }
.bad { color: var(--bad); }
.scroll { overflow-x: auto; }
table { border-collapse: collapse; font-size: 13px; margin: 8px 0; }
th, td { border: 1px solid var(--line); padding: 4px 8px; text-align: left; vertical-align: top; }
th { background: var(--head); }
.num { text-align: right; font-variant-numeric: tabular-nums; }
pre { background: var(--head); border: 1px solid var(--line); padding: 8px 10px; overflow-x: auto; font-size: 12px; }
code { font-size: 0.9em; }
dl.facts { display: grid; grid-template-columns: max-content 1fr; gap: 2px 16px; font-size: 13px; }
dl.facts dt { color: var(--muted); }
dl.facts dd { margin: 0; overflow-wrap: anywhere; }
"""


def page(title, body):
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>{STYLE}</style>
</head>
<body>
{body}
</body>
</html>
"""


def write(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def commit_link(code):
    """The Code line of machine.txt ("`abc1234` (subject)") with the hash linked to the commit."""
    m = re.match(r"`([0-9a-f]+)`\s*(.*)", code)
    if not m:
        return inline(code)
    sha = html.escape(m.group(1))
    link = f'<a href="{html.escape(REPO_URL)}/commit/{sha}"><code>{sha}</code></a>' if REPO_URL else f"<code>{sha}</code>"
    return f"{link} {html.escape(m.group(2))}"


def run(out, jobs):
    os.makedirs(out, exist_ok=True)
    rows, summary, sections, commits = [], [], [], set()
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
            profile = f"{arch}/jvm-profile.html"
        row = {
            "architecture": arch, "cpu": m.get("CPU model", "?"), "cpus": m.get("CPUs", "?"),
            "blast_wall": wall.get("blast"), "jvm_wall": wall.get("jvm"), "ref_wall": wall.get("ref"),
            "wall": ratio(wall.get("jvm"), wall.get("blast")),
            "cycles": ratio(jvm.get("cycles"), blast.get("cycles")),
            "instructions": ratio(jvm.get("instructions"), blast.get("instructions")),
            "identical": identical, "profile": profile,
        }
        rows.append(row)
        link = f'<a href="{html.escape(profile)}">flame graph</a>' if profile else ""
        summary.append(
            "<tr>"
            f"<td>{html.escape(arch)}</td><td>{html.escape(row['cpu'])}</td>"
            f"<td class=num>{html.escape(row['cpus'])}</td>"
            f"<td class=num>{row['blast_wall'] or '–'}</td><td class=num>{row['jvm_wall'] or '–'}</td>"
            f"<td class=num><b>{times(row['wall'])}</b></td>"
            f"<td class=num>{times(row['cycles'])}</td><td class=num>{times(row['instructions'])}</td>"
            + (f"<td class=num>{row['ref_wall'] or '–'}</td>" if "ref" in wall else "")
            + f"<td>{'identical' if identical else '<b class=bad>different</b>'}</td>"
            f"<td>{link}</td></tr>")
        facts = "".join(f"<dt>{html.escape(k)}</dt><dd>{inline(v)}</dd>" for k, v in m.items())
        sections.append(
            f'<section><h2 id="{html.escape(arch)}">{html.escape(arch)}: {html.escape(row["cpu"])}</h2>'
            f'<dl class="facts">{facts}</dl>'
            + markdown(read(os.path.join(job, "report.md")), shift=1)
            + "</section>")
    has_ref = any(r["ref_wall"] is not None for r in rows)
    code = sorted(c for c in commits if c)
    date = os.environ.get("RUN_DATE") or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    body = f"""<p class="muted"><a href="../../">All runs</a></p>
<h1>blastp4j against BLAST+</h1>
<p class="muted">{date} · code {", ".join(commit_link(c) for c in code) or "?"}{f' · <a href="{html.escape(RUN_URL)}">workflow run</a>' if RUN_URL else ''}.
Hosted runners are virtual machines and differ from run to run: compare the ratios within one row, not times across runs.</p>
<div class="scroll"><table>
<thead><tr><th>architecture</th><th>CPU model</th><th class=num>CPUs</th><th class=num>BLAST+ wall s</th>
<th class=num>blastp4j wall s</th><th class=num>wall</th><th class=num>cycles</th><th class=num>instructions</th>
{'<th class=num>reference wall s</th>' if has_ref else ''}<th>output</th><th></th></tr></thead>
<tbody>{''.join(summary)}</tbody>
</table></div>
<p class="muted">Wall times are medians of the full-search rounds; the ratios are blastp4j over BLAST+.</p>
{''.join(sections)}"""
    write(os.path.join(out, "index.html"), page("blastp4j vs BLAST+, " + date, body))
    write(os.path.join(out, "summary.json"), json.dumps(
        {"run_id": os.environ.get("RUN_ID", os.path.basename(os.path.normpath(out))), "run_url": RUN_URL,
         "date": date, "event": os.environ.get("RUN_EVENT", ""), "code": code, "rows": rows}, indent=1))


def index(site):
    runs = []
    for d in os.listdir(os.path.join(site, "runs")) if os.path.isdir(os.path.join(site, "runs")) else []:
        text = read(os.path.join(site, "runs", d, "summary.json"))
        if text:
            runs.append((d, json.loads(text)))
    runs.sort(key=lambda r: r[1].get("date", ""), reverse=True)
    trs = []
    for d, s in runs:
        rows = s.get("rows") or [{}]
        for k, r in enumerate(rows):
            tr = "<tr>"
            if k == 0:
                n = len(rows)
                event = {"schedule": "daily", "workflow_dispatch": "manual"}.get(s.get("event", ""), "")
                tr += (f'<td rowspan={n}><a href="runs/{html.escape(d)}/">{html.escape(s.get("date", d))}</a></td>'
                       f'<td rowspan={n}>{event}</td>'
                       f'<td rowspan={n}>{", ".join(commit_link(c) for c in s.get("code", []))}</td>')
            tr += (f"<td>{html.escape(r.get('architecture', '?'))}</td><td>{html.escape(r.get('cpu', '?'))}</td>"
                   f"<td class=num>{html.escape(str(r.get('cpus', '?')))}</td>"
                   f"<td class=num><b>{times(r.get('wall'))}</b></td><td class=num>{times(r.get('cycles'))}</td>"
                   f"<td class=num>{times(r.get('instructions'))}</td>"
                   f"<td>{'identical' if r.get('identical') else '<b class=bad>different</b>'}</td></tr>")
            trs.append(tr)
    body = f"""<h1>blastp4j against BLAST+</h1>
<p class="muted">Every run of the architecture comparison workflow, newest first. The ratios are blastp4j over BLAST+ in
the full search, measured in the same job. Hosted runners are virtual machines on varying CPUs, so compare runs on the
same CPU model.</p>
<div class="scroll"><table>
<thead><tr><th>run</th><th>trigger</th><th>code</th><th>architecture</th><th>CPU model</th><th class=num>CPUs</th><th class=num>wall</th>
<th class=num>cycles</th><th class=num>instructions</th><th>output</th></tr></thead>
<tbody>{''.join(trs)}</tbody>
</table></div>"""
    write(os.path.join(site, "index.html"), page("blastp4j vs BLAST+", body))


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3:])
    elif len(sys.argv) == 3 and sys.argv[1] == "index":
        index(sys.argv[2])
    else:
        sys.exit(__doc__)
