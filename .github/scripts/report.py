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
numbers (plus, per job, the measured commits and the wall time of each round per system), and each flame graph as
OUT_DIR/<architecture>/jvm-profile.html. The Markdown is the small subset the workflow writes (headings, paragraphs, fenced code, pipe tables, inline code), so no Markdown library is needed.

"index" writes SITE_DIR/index.html: every run in SITE_DIR/runs/*/summary.json, newest first, and above them violin
plots of the full search's rounds per measured commit (drawn in the browser with Plotly.js; only the newest commit is
open). The rounds come from summary.json, or, for runs from before they were kept there, from the run's index.html.

Settings (environment, for "run"): RUN_ID, RUN_URL, REPO_URL (for links to the commits), RUN_DATE (default now),
RUN_EVENT (what started the run: schedule, push or workflow_dispatch).
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


def search_rounds(text):
    """The wall time of each round per system ({"blast": [s, ...], "jvm": [...], "ref": [...]}, None for a missing
    round) and the measured commits ({"jvm": hash, "ref": hash}) from bin/compare-blastp.sh's output."""
    walls, commits = {}, {}
    for line in text.splitlines():
        line = line.strip()
        if m := re.match(r"round (\d+) (\w+): ([\d.]+)\s", line + " "):
            walls.setdefault(m.group(2), {})[int(m.group(1))] = float(m.group(3))
        elif re.fullmatch(r"\w+: [0-9a-f]{7,40}(, \w+: [0-9a-f]{7,40})*", line):
            commits.update(re.findall(r"(\w+): ([0-9a-f]+)", line))
    n = max((max(w) for w in walls.values()), default=0)
    return {system: [w.get(k) for k in range(1, n + 1)] for system, w in walls.items()}, commits


def code_hash(code):
    """The hash of a Code line of machine.txt ("`abc1234` (subject)")."""
    m = re.match(r"`([0-9a-f]+)`", code or "")
    return m.group(1) if m else ""


def ratio(a, b):
    return a / b if a and b else None


def times(r):
    return f"{r:.2f}×" if r else "–"


STYLE = """
:root { --bg: #fdfdfc; --fg: #1d1d1b; --muted: #66645e; --line: #dddad3; --head: #f1efe9; --accent: #2f5d8a; --bad: #b3261e;
  --arm: #2a78d6; --x86: #eb6834; }
@media (prefers-color-scheme: dark) {
  :root { --bg: #17181a; --fg: #e6e4df; --muted: #a09d96; --line: #34363a; --head: #222428; --accent: #8db7e0; --bad: #f2867d;
    --arm: #3987e5; --x86: #d95926; }
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
details.commit { border: 1px solid var(--line); border-radius: 4px; margin: 8px 0; padding: 0 12px; }
details.commit > summary { cursor: pointer; padding: 8px 0; }
details.commit[open] > summary { border-bottom: 1px solid var(--line); }
.plots { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 420px), 1fr)); gap: 8px; }
.plot { min-height: 360px; min-width: 0; overflow: hidden; }
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
        search = read(os.path.join(job, "full-search.txt"))
        wall, counters, identical = full_search(search)
        blast, jvm = counters.get("blast", {}), counters.get("jvm", {})
        rounds, tested = search_rounds(search)
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
            "jvm_commit": tested.get("jvm") or code_hash(m.get("Code", "")), "ref_commit": tested.get("ref", ""),
            "rounds": rounds,
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


def backfill(page_html):
    """{architecture: (rounds, jvm commit, reference commit)} from a run page whose summary.json predates the rounds:
    each job's section has the full-search output in a <pre> block. Read-only."""
    out = {}
    for arch, body in re.findall(r'<section><h2 id="([^"]*)">(.*?)</section>', page_html, re.S):
        text = "\n".join(html.unescape(p) for p in re.findall(r"<pre>(.*?)</pre>", body, re.S))
        rounds, commits = search_rounds(text)
        code = re.search(r"<dt>Code</dt><dd><code>([0-9a-f]+)</code>", body)
        out[html.unescape(arch)] = (rounds, commits.get("jvm") or (code.group(1) if code else ""), commits.get("ref", ""))
    return out


def measured(rounds, system):
    return any(x is not None for x in rounds.get(system) or [])


def per_commit(site, runs):
    """The runs grouped by the blastp4j commit they measured, the newest commit (by its latest run) first. Each group
    has the hash, the Code line, the runs and, per job with rounds, the wall times of each round per system. Also the
    number of runs whose rounds were read from their pages (summary.json without rounds)."""
    groups, backfilled = {}, set()
    for d, s in runs:
        pages = None
        for r in s.get("rows") or []:
            rounds, jvm, ref = r.get("rounds"), r.get("jvm_commit", ""), r.get("ref_commit", "")
            if rounds is None:
                if pages is None:
                    pages = backfill(read(os.path.join(site, "runs", d, "index.html")))
                rounds, jvm_page, ref_page = pages.get(r.get("architecture"), ({}, "", ""))
                jvm, ref = jvm or jvm_page, ref or ref_page
                if rounds:
                    backfilled.add(d)
            jvm = (jvm or next((code_hash(c) for c in s.get("code", []) if code_hash(c)), ""))[:7]
            if not jvm:
                continue
            g = groups.setdefault(jvm, {"hash": jvm, "code": "", "date": "", "runs": {}, "jobs": []})
            g["code"] = g["code"] or next((c for c in s.get("code", []) if code_hash(c).startswith(jvm)), "")
            g["date"] = max(g["date"], s.get("date", ""))
            run = g["runs"].setdefault(d, {"id": d, "date": s.get("date", d), "measured": False})
            if measured(rounds, "jvm") and (measured(rounds, "blast") or measured(rounds, "ref")):
                run["measured"] = True
                job = {"run": d, "arch": r.get("architecture", "?"), "cpu": r.get("cpu", "?"), "ref_commit": ref[:7]}
                job.update({k: rounds[k] for k in ("blast", "jvm", "ref") if k in rounds})
                g["jobs"].append(job)
    ordered = sorted(groups.values(), key=lambda g: g["date"], reverse=True)
    for g in ordered:
        g["runs"] = sorted(g["runs"].values(), key=lambda r: r["date"], reverse=True)
    return ordered, len(backfilled)


PLOTLY = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"

# Draws the violin plots of a commit when its <details> is opened (the newest commit's on load), with the colors of
# the page's CSS variables.
PLOT_SCRIPT = """
(function () {
  const data = JSON.parse(document.getElementById("perf-data").textContent);
  const NAMES = {blast: "BLAST+", ref: "reference"};
  const METRICS = [
    ["blast", "blastp4j / BLAST+, wall time per round"],
    ["ref", "blastp4j / reference commit, wall time per round"],
  ];
  const ARCH = {aarch64: "--arm", x86_64: "--x86"};
  function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function hex(color) { return /^#[0-9a-f]{6}$/i.test(color) ? color : "#888888"; }

  function traces(commit, den) {
    const by = {};
    for (const j of commit.jobs) {
      const a = j.jvm || [], b = j[den] || [];
      for (let i = 0; i < Math.min(a.length, b.length); i++) {
        if (a[i] == null || !b[i]) continue;
        const t = by[j.arch] || (by[j.arch] = {y: [], text: []});
        t.y.push(a[i] / b[i]);
        t.text.push("run " + j.run + ", round " + (i + 1) + "<br>" + j.cpu + "<br>blastp4j " + a[i].toFixed(2)
          + " s, " + NAMES[den] + (den === "ref" && j.ref_commit ? " " + j.ref_commit : "") + " " + b[i].toFixed(2) + " s");
      }
    }
    const fg = css("--fg"), bg = css("--bg");
    return Object.keys(by).sort().map(function (arch) {
      const color = hex(css(ARCH[arch] || "--muted"));
      return {
        type: "violin", orientation: "v", name: arch, x: by[arch].y.map(function () { return arch; }), y: by[arch].y,
        text: by[arch].text, hoveron: "points", hovertemplate: "%{text}<br>ratio %{y:.3f}<extra></extra>",
        points: "all", pointpos: 0, jitter: 0.6, spanmode: "hard", width: 0.8,
        line: {color: color, width: 2}, fillcolor: color + "33",
        marker: {color: color, size: 8, line: {color: bg, width: 1}},
        box: {visible: true, width: 0.12, fillcolor: bg, line: {color: fg, width: 1}},
      };
    });
  }

  function layout(title, den, n) {
    const fg = css("--fg"), bg = css("--bg"), line = css("--line");
    return {
      title: {text: title, font: {size: 14}, x: 0, xanchor: "left"},
      paper_bgcolor: bg, plot_bgcolor: bg, font: {color: fg, family: "system-ui, sans-serif", size: 12},
      showlegend: false, height: 360, margin: {l: 60, r: 12, t: 40, b: 36},
      xaxis: {type: "category", range: [-0.5, n - 0.5], linecolor: line, tickcolor: line},
      yaxis: {title: {text: "below 1: blastp4j faster than " + NAMES[den]}, gridcolor: line, linecolor: line,
              zeroline: false, tickformat: ".2f"},
      shapes: [{type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: 1, y1: 1, layer: "below",
                line: {color: css("--muted"), width: 1.5, dash: "dash"}}],
      hoverlabel: {bgcolor: css("--head"), bordercolor: line, font: {color: fg}},
    };
  }

  function render(details) {
    const box = details.querySelector(".plots");
    if (typeof Plotly === "undefined") {
      box.textContent = "The plots need Plotly.js from cdn.jsdelivr.net, which did not load.";
      return;
    }
    for (const old of box.querySelectorAll(".plot")) Plotly.purge(old);
    box.textContent = "";
    const commit = data[Number(details.dataset.k)];
    const plots = [];
    for (const [den, title] of METRICS) {
      const t = traces(commit, den);
      if (!t.length) continue;
      const div = document.createElement("div");
      div.className = "plot";
      box.appendChild(div);
      plots.push([div, t, layout(title, den, t.length)]);
    }
    // Draw once all the plots are in the grid, so that each gets its column's width.
    for (const [div, t, l] of plots) {
      Plotly.newPlot(div, t, l, {responsive: true, displaylogo: false, modeBarButtonsToRemove: ["select2d", "lasso2d"]});
    }
    details.dataset.rendered = "1";
  }

  const all = document.querySelectorAll("details.commit");
  for (const d of all) {
    d.addEventListener("toggle", function () { if (d.open && !d.dataset.rendered) render(d); });
    if (d.open && !d.dataset.rendered) render(d);
  }
  // Draw again with the other colors when the color scheme changes.
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", function () {
    for (const d of all) {
      if (!d.dataset.rendered) continue;
      if (d.open) render(d); else delete d.dataset.rendered;
    }
  });
})();
"""


def commit_sections(groups):
    """The index's per-commit sections, only the newest open, and the plots' data as JSON for a <script> element."""
    out, data = [], []
    for g in groups:
        if not g["jobs"]:
            continue
        k = len(data)
        data.append({"hash": g["hash"], "jobs": g["jobs"]})
        n = sum(1 for r in g["runs"] if r["measured"])
        links = ", ".join(f'<a href="runs/{html.escape(r["id"])}/">{html.escape(r["date"])}</a>'
                          + ("" if r["measured"] else " (no rounds)") for r in g["runs"])
        refs = sorted({j["ref_commit"] for j in g["jobs"] if j["ref_commit"] and "ref" in j})
        title = commit_link(g["code"]) if g["code"] else f"<code>{html.escape(g['hash'])}</code>"
        out.append(
            f'<details class="commit" data-k="{k}"{" open" if k == 0 else ""}>'
            f'<summary>{title} <span class="muted">· {n} run{"" if n == 1 else "s"} · latest '
            f'{html.escape(g["date"])}</span></summary>'
            f'<p class="muted">Runs: {links}.'
            + (f' Reference: {", ".join(f"<code>{html.escape(r)}</code>" for r in refs)}.' if refs else "")
            + '</p><div class="plots"></div></details>')
    return "\n".join(out), json.dumps(data, separators=(",", ":")).replace("</", "<\\/")


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
                event = {"schedule": "daily", "push": "push", "workflow_dispatch": "manual"}.get(s.get("event", ""), "")
                tr += (f'<td rowspan={n}><a href="runs/{html.escape(d)}/">{html.escape(s.get("date", d))}</a></td>'
                       f'<td rowspan={n}>{event}</td>'
                       f'<td rowspan={n}>{", ".join(commit_link(c) for c in s.get("code", []))}</td>')
            tr += (f"<td>{html.escape(r.get('architecture', '?'))}</td><td>{html.escape(r.get('cpu', '?'))}</td>"
                   f"<td class=num>{html.escape(str(r.get('cpus', '?')))}</td>"
                   f"<td class=num><b>{times(r.get('wall'))}</b></td><td class=num>{times(r.get('cycles'))}</td>"
                   f"<td class=num>{times(r.get('instructions'))}</td>"
                   f"<td>{'identical' if r.get('identical') else '<b class=bad>different</b>'}</td></tr>")
            trs.append(tr)
    groups, backfilled = per_commit(site, runs)
    sections, data = commit_sections(groups)
    body = f"""<h1>blastp4j against BLAST+</h1>
<p class="muted">Every run of the architecture comparison workflow. The ratios are blastp4j over BLAST+ in the full
search, measured in the same job. Hosted runners are virtual machines on varying CPUs, so compare runs on the same CPU
model.</p>
<h2>Per commit</h2>
<p class="muted">Each point is one round of the full search: blastp4j's wall time over BLAST+'s, or over the reference
commit's, in the same round of the same job, which keeps the different runner CPUs comparable. All rounds of all runs
of a commit, newest commit first; open a commit to draw its plots.</p>
{sections or '<p class="muted">No runs with rounds yet.</p>'}
<h2>Runs</h2>
<p class="muted">Newest first.</p>
<div class="scroll"><table>
<thead><tr><th>run</th><th>trigger</th><th>code</th><th>architecture</th><th>CPU model</th><th class=num>CPUs</th><th class=num>wall</th>
<th class=num>cycles</th><th class=num>instructions</th><th>output</th></tr></thead>
<tbody>{''.join(trs)}</tbody>
</table></div>
<script type="application/json" id="perf-data">{data}</script>
<script src="{PLOTLY}" charset="utf-8"></script>
<script>{PLOT_SCRIPT}</script>"""
    write(os.path.join(site, "index.html"), page("blastp4j vs BLAST+", body))
    print(f"index: {len(runs)} runs; {sum(1 for g in groups if g['jobs'])} commits with rounds of {len(groups)}; "
          f"the rounds of {backfilled} runs read from their pages", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3:])
    elif len(sys.argv) == 3 and sys.argv[1] == "index":
        index(sys.argv[2])
    else:
        sys.exit(__doc__)
