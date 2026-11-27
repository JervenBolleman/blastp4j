#!/usr/bin/env python3
"""Builds the HTML reports of the architecture comparison runs (.github/workflows/arch-compare.yml).

Usage: report.py run OUT_DIR JOB_DIR...
       report.py index SITE_DIR

"run" builds one run's report. Each JOB_DIR is one job's report directory:
  machine.txt      key: value lines (architecture, CPU model, cores, memory, runner image, JDK, BLAST+, commits)
  report.md        the job's Markdown sections, as written to the step summary
  full-search.txt  the output of bin/compare-blastp.sh
  <system>-profile.html the async-profiler flame graphs of jvm, almost and fast (optional)
It writes OUT_DIR/index.html, with a summary table of the full searches on top, OUT_DIR/summary.json with that table's
numbers (plus, per job, the measured commits and the wall time of each round per system), and each flame graph as
OUT_DIR/<architecture>/<system>-profile.html. The Markdown is the small subset the workflow writes (headings, paragraphs, fenced code, pipe tables, inline code), so no Markdown library is needed.

"index" writes SITE_DIR/index.html: every run in SITE_DIR/runs/*/summary.json, newest first, above them violin
plots of the full search's rounds per measured commit (drawn in the browser with Plotly.js; only the newest commit is
open), and above those one chart of all the runs over time, six lines: each of blastp4j, almostblastp and fastblastp
over the BLAST+ task it must match, on each architecture. The rounds come from summary.json, or, for runs from before
they were kept there, from the run's index.html.

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


# The database and the queries changed with the push of 2026-11-22 20:30 UTC (arch-compare.yml): before it the
# Phaseolus vulgaris proteome and its first 500 proteins, after it the human one and 750 evenly spaced proteins, which
# is about three times the instructions. A ratio to BLAST+ on one workload is not the same quantity as on the other, so
# the history chart draws the boundary. Runs from before machine.txt named the proteome are told apart by their date;
# the first run of the human workload is 2026-11-22 21:50 UTC.
HUMAN_FROM = "2026-11-22 21:00"


def workload(date, named=""):
    """The proteome a run measured: the name its machine.txt gave, or the one its date implies."""
    return named or ("UP000005640_9606" if date >= HUMAN_FROM else "UP000000226_3885")


# bin/jvm.options' [avx512] stanza gives 512-bit vectors to AMD family 26 and later (Zen 5, Zen 6) and to Intel parts
# with avx512f (Ice Lake-SP and later Xeons); every other x86-64 runner keeps 256, because 512 was slower on Zen 4 and
# on the older Xeons. Such a job therefore measures different code from the other x86-64 jobs, so the history chart
# marks its points. The signal is the options the job ran with, which machine.txt carries from bin/jvm-options; jobs
# from before it did are told apart by their CPU model, as the workload is by its date, and the only model the runners
# had offered by then was Zen 5's. Named for the width and not for a CPU: the stanza was Zen 5 alone until 2026-11-26,
# and summaries written before that call the key zen5, which read_jobs still accepts.
WIDE512_MODELS = ("9V45",)


def wide512(cpu, options=""):
    """Whether a job ran with 512-bit vectors: what its own options say, or what its CPU model implies."""
    if "MaxVectorSize" in options:
        return "MaxVectorSize=64" in options
    return any(model in cpu for model in WIDE512_MODELS)


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
#history { min-height: 460px; min-width: 0; overflow: hidden; margin-top: 8px; }
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


def subject(code, limit=70):
    """The subject of a Code line of machine.txt ("`abc1234` (subject)"), escaped for the charts' hover text."""
    m = re.match(r"`[0-9a-f]+`\s*\((.*)\)\s*$", (code or "").strip())
    text = m.group(1) if m else ""
    return html.escape(text if len(text) <= limit else text[:limit - 1] + "…", quote=False)


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
        profiles = {}
        for system in ("jvm", "almost", "fast"):
            name = f"{system}-profile.html"
            if os.path.exists(os.path.join(job, name)):
                os.makedirs(os.path.join(out, arch), exist_ok=True)
                shutil.copy(os.path.join(job, name), os.path.join(out, arch, name))
                profiles[system] = f"{arch}/{name}"
        row = {
            "architecture": arch, "cpu": m.get("CPU model", "?"), "cpus": m.get("CPUs", "?"),
            "wide512": wide512(m.get("CPU model", ""), m.get("JVM options", "")),
            "blast_wall": wall.get("blast"), "jvm_wall": wall.get("jvm"), "ref_wall": wall.get("ref"),
            "blastfast_wall": wall.get("blastfast"), "fast_wall": wall.get("fast"),
            "fast": ratio(wall.get("fast"), wall.get("blastfast")),
            "wall": ratio(wall.get("jvm"), wall.get("blast")),
            "cycles": ratio(jvm.get("cycles"), blast.get("cycles")),
            "instructions": ratio(jvm.get("instructions"), blast.get("instructions")),
            "identical": identical, "profile": profiles.get("jvm", ""),
            "proteome": m.get("Proteome", ""),
            "jvm_commit": tested.get("jvm") or code_hash(m.get("Code", "")), "ref_commit": tested.get("ref", ""),
            "rounds": rounds,
        }
        rows.append(row)
        link = " ".join(f'<a href="{html.escape(p)}">{"flame graph" if s == "jvm" else s}</a>' for s, p in profiles.items())
        summary.append(
            "<tr>"
            f"<td>{html.escape(arch)}</td><td>{html.escape(row['cpu'])}</td>"
            f"<td class=num>{html.escape(row['cpus'])}</td>"
            f"<td class=num>{row['blast_wall'] or '–'}</td><td class=num>{row['jvm_wall'] or '–'}</td>"
            f"<td class=num><b>{times(row['wall'])}</b></td>"
            f"<td class=num>{times(row['cycles'])}</td><td class=num>{times(row['instructions'])}</td>"
            + (f"<td class=num>{row['ref_wall'] or '–'}</td>" if "ref" in wall else "")
            + (f"<td class=num>{row['blastfast_wall'] or '–'}</td><td class=num>{row['fast_wall'] or '–'}</td>"
               f"<td class=num><b>{times(row['fast'])}</b></td>" if "fast" in wall else "")
            + f"<td>{'identical' if identical else '<b class=bad>different</b>'}</td>"
            f"<td>{link}</td></tr>")
        facts = "".join(f"<dt>{html.escape(k)}</dt><dd>{inline(v)}</dd>" for k, v in m.items())
        sections.append(
            f'<section><h2 id="{html.escape(arch)}">{html.escape(arch)}: {html.escape(row["cpu"])}</h2>'
            f'<dl class="facts">{facts}</dl>'
            + markdown(read(os.path.join(job, "report.md")), shift=1)
            + "</section>")
    has_ref = any(r["ref_wall"] is not None for r in rows)
    has_fast = any(r["fast_wall"] is not None for r in rows)
    code = sorted(c for c in commits if c)
    date = os.environ.get("RUN_DATE") or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    body = f"""<p class="muted"><a href="../../">All runs</a></p>
<h1>blastp4j against BLAST+</h1>
<p class="muted">{date} · code {", ".join(commit_link(c) for c in code) or "?"}{f' · <a href="{html.escape(RUN_URL)}">workflow run</a>' if RUN_URL else ''}.
Hosted runners are virtual machines and differ from run to run: compare the ratios within one row, not times across runs.</p>
<div class="scroll"><table>
<thead><tr><th>architecture</th><th>CPU model</th><th class=num>CPUs</th><th class=num>BLAST+ wall s</th>
<th class=num>blastp4j wall s</th><th class=num>wall</th><th class=num>cycles</th><th class=num>instructions</th>
{'<th class=num>reference wall s</th>' if has_ref else ''}{'<th class=num>BLAST+ blastp-fast wall s</th><th class=num>fastblastp wall s</th><th class=num>fast wall</th>' if has_fast else ''}<th>output</th><th></th></tr></thead>
<tbody>{''.join(summary)}</tbody>
</table></div>
<p class="muted">Wall times are medians of the full-search rounds; the ratios are blastp4j over BLAST+, and fastblastp
over BLAST+ with <code>-task blastp-fast</code>. Each output is compared with BLAST+'s for its task.</p>
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
                job = {"run": d, "date": s.get("date", d), "arch": r.get("architecture", "?"),
                       "cpu": r.get("cpu", "?"), "ref_commit": ref[:7],
                       # "zen5" is what summaries written before 2026-11-27 call it.
                       "wide512": bool(r.get("wide512", r.get("zen5", wide512(r.get("cpu", ""))))),
                       "workload": workload(s.get("date", ""), r.get("proteome", ""))}
                job.update({k: rounds[k] for k in ("blast", "jvm", "ref", "almost", "blastfast", "fast", "fastref") if k in rounds})
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
  const NAMES = {blast: "BLAST+", ref: "reference", jvm: "blastp4j", almost: "almostblastp",
    blastfast: "BLAST+ blastp-fast", fast: "fastblastp", fastref: "reference fastblastp"};
  // Numerator, denominator, title.
  const METRICS = [
    ["jvm", "blast", "blastp4j / BLAST+, wall time per round"],
    ["jvm", "ref", "blastp4j / reference commit, wall time per round"],
    ["almost", "blast", "almostblastp / BLAST+, wall time per round"],
    ["fast", "blastfast", "fastblastp / BLAST+ blastp-fast, wall time per round"],
    ["fast", "fastref", "fastblastp / reference commit's, wall time per round"],
  ];
  const ARCH = {aarch64: "--arm", x86_64: "--x86"};
  function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function hex(color) { return /^#[0-9a-f]{6}$/i.test(color) ? color : "#888888"; }

  function traces(commit, num, den) {
    const by = {};
    for (const j of commit.jobs) {
      const a = j[num] || [], b = j[den] || [];
      for (let i = 0; i < Math.min(a.length, b.length); i++) {
        if (a[i] == null || !b[i]) continue;
        const t = by[j.arch] || (by[j.arch] = {y: [], text: []});
        t.y.push(a[i] / b[i]);
        t.text.push("run " + j.run + ", round " + (i + 1) + "<br>" + j.cpu + "<br>" + NAMES[num] + " " + a[i].toFixed(2)
          + " s, " + NAMES[den] + ((den === "ref" || den === "fastref") && j.ref_commit ? " " + j.ref_commit : "") + " " + b[i].toFixed(2) + " s");
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

  function layout(title, num, den, t) {
    const fg = css("--fg"), bg = css("--bg"), line = css("--line");
    // The ratio axis starts at 0, so that heights read as proportions.
    const ys = [].concat(...t.map(function (x) { return x.y; }));
    const hi = Math.max(1, ...ys) * 1.05;
    const n = t.length;
    return {
      title: {text: title, font: {size: 14}, x: 0, xanchor: "left"},
      paper_bgcolor: bg, plot_bgcolor: bg, font: {color: fg, family: "system-ui, sans-serif", size: 12},
      showlegend: false, height: 360, margin: {l: 60, r: 12, t: 40, b: 36},
      xaxis: {type: "category", range: [-0.5, n - 0.5], linecolor: line, tickcolor: line},
      yaxis: {title: {text: "below 1: " + NAMES[num] + " faster than " + NAMES[den]}, gridcolor: line, linecolor: line,
              zeroline: false, range: [0, hi]},
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
    for (const [num, den, title] of METRICS) {
      const t = traces(commit, num, den);
      if (!t.length) continue;
      const div = document.createElement("div");
      div.className = "plot";
      box.appendChild(div);
      plots.push([div, t, layout(title, num, den, t)]);
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


# Draws the history chart: one line per architecture and program, a point per job, over the runs' dates. It reads the
# same data as the violin plots, which carries each job's date, CPU, workload and rounds.
HISTORY_SCRIPT = """
(function () {
  const el = document.getElementById("history");
  const data = JSON.parse(document.getElementById("perf-data").textContent);
  const NAMES = {blast: "BLAST+", jvm: "blastp4j", almost: "almostblastp",
    blastfast: "BLAST+ blastp-fast", fast: "fastblastp"};
  // Each program over the BLAST+ task whose output it must match, and the dash that tells the three apart within one
  // architecture's color.
  const SERIES = [["jvm", "blast", "solid"], ["almost", "blast", "dash"], ["fast", "blastfast", "dot"]];
  const ARCH = {aarch64: "--arm", x86_64: "--x86"};
  function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function hex(color) { return /^#[0-9a-f]{6}$/i.test(color) ? color : "#888888"; }
  function when(x) { return new Date(x.replace(" ", "T") + ":00Z"); }
  function iso(d) { return d.toISOString().replace("T", " ").slice(0, 16); }
  function median(xs) {
    const v = xs.slice().sort(function (a, b) { return a - b; });
    const h = v.length >> 1;
    return v.length % 2 ? v[h] : (v[h - 1] + v[h]) / 2;
  }

  // One point per job and metric: the median of its rounds' ratios, with the bar over all of them. Both times of a
  // ratio are from the same round of the same job, so the runner's CPU cancels out of it.
  function points(num, den) {
    const by = {};
    for (const commit of data) {
      for (const j of commit.jobs) {
        if (!(j.arch in ARCH)) continue;
        const a = j[num] || [], b = j[den] || [], rs = [];
        for (let i = 0; i < Math.min(a.length, b.length); i++) {
          if (a[i] != null && b[i]) rs.push(a[i] / b[i]);
        }
        if (!rs.length || !j.date) continue;
        const lo = Math.min.apply(null, rs), hi = Math.max.apply(null, rs);
        (by[j.arch] || (by[j.arch] = [])).push({
          x: j.date.replace(" UTC", ""), y: median(rs), lo: lo, hi: hi, run: j.run, workload: j.workload || "",
          wide512: !!j.wide512,
          text: "<b>" + commit.hash + "</b> " + (commit.subject || "") + "<br>" + j.cpu
            + (j.wide512 ? " · 512-bit vectors" : "") + " · " + j.arch
            + " · run " + j.run + " · " + j.date + "<br>" + NAMES[num] + " over " + NAMES[den] + ", "
            + rs.length + (rs.length === 1 ? " round " : " rounds ") + lo.toFixed(3) + "–" + hi.toFixed(3)
            + "<br>" + (j.workload || "") + " · click for this run's report",
        });
      }
    }
    for (const arch in by) by[arch].sort(function (p, q) { return p.x < q.x ? -1 : p.x > q.x ? 1 : 0; });
    return by;
  }

  function build() {
    const bg = css("--bg"), out = [];
    let anyWide512 = false;
    for (const [num, den, dash] of SERIES) {
      const by = points(num, den);
      for (const arch of Object.keys(by).sort()) {
        const color = hex(css(ARCH[arch])), ps = by[arch];
        anyWide512 = anyWide512 || ps.some(function (p) { return p.wide512; });
        out.push({points: ps, trace: {
          type: "scatter", mode: "lines+markers", name: arch + " " + NAMES[num], legendgroup: arch,
          x: ps.map(function (p) { return p.x; }), y: ps.map(function (p) { return p.y; }),
          text: ps.map(function (p) { return p.text; }), customdata: ps.map(function (p) { return p.run; }),
          hovertemplate: "%{text}<br>ratio <b>%{y:.3f}</b><extra></extra>",
          line: {color: color, width: 2, dash: dash},
          // A hollow diamond where the job ran with 512-bit vectors: the color is the architecture and
          // the dash the program, so the symbol is what is left to say the CPU changed the code. Hollow rather than
          // filled because these points sit in dense clusters, where one filled shape among others of the same color
          // is not something an eye picks out; the open one draws the page through it. It needs the series color on
          // its outline, where a filled marker takes the background as a ring to separate it from its neighbours.
          marker: {color: color, size: ps.map(function (p) { return p.wide512 ? 11 : 8; }),
            symbol: ps.map(function (p) { return p.wide512 ? "diamond-open" : "circle"; }),
            line: {width: ps.map(function (p) { return p.wide512 ? 2 : 1; }),
              color: ps.map(function (p) { return p.wide512 ? color : bg; })}},
          error_y: {type: "data", symmetric: false, color: color, thickness: 1, width: 3,
            array: ps.map(function (p) { return p.hi - p.y; }),
            arrayminus: ps.map(function (p) { return p.y - p.lo; })},
        }});
      }
    }
    // A key for the symbol. The legend's own entries are one per architecture and program, and their markers take the
    // symbol of whichever point happens to come first, so the one thing that says what a diamond means has to be an
    // entry of its own: a trace with no data, in the muted color, since the symbol belongs to no architecture.
    if (out.length && anyWide512) {
      const muted = css("--muted");
      out.push({points: [], trace: {
        type: "scatter", mode: "markers", name: "(512-bit vectors)", x: [null], y: [null],
        marker: {color: muted, size: 11, symbol: "diamond-open", line: {color: muted, width: 2}},
        hoverinfo: "skip",
      }});
    }
    return out;
  }

  // Where the workload changed between two runs: the ratios either side are of different searches.
  function boundaries(ps) {
    const seen = {};
    for (const p of ps) seen[p.x] = p.workload;
    const xs = Object.keys(seen).sort(), out = [];
    for (let i = 1; i < xs.length; i++) {
      if (seen[xs[i]] && seen[xs[i - 1]] && seen[xs[i]] !== seen[xs[i - 1]]) {
        out.push([new Date((when(xs[i - 1]).getTime() + when(xs[i]).getTime()) / 2), seen[xs[i]]]);
      }
    }
    return out;
  }

  function layout(ps) {
    const fg = css("--fg"), bg = css("--bg"), line = css("--line"), muted = css("--muted");
    // The range holds every bar and 1, the parity BLAST+ marks; it does not start at 0, which would flatten the lines.
    let lo = 1, hi = 1;
    for (const p of ps) { lo = Math.min(lo, p.lo); hi = Math.max(hi, p.hi); }
    const pad = (hi - lo) * 0.08 || 0.05;
    const marks = boundaries(ps);
    return {
      paper_bgcolor: bg, plot_bgcolor: bg, font: {color: fg, family: "system-ui, sans-serif", size: 12},
      height: 460, margin: {l: 60, r: 16, t: 40, b: 100}, hovermode: "closest",
      legend: {orientation: "h", x: 0, y: -0.16, yanchor: "top", font: {size: 11}, groupclick: "toggleitem"},
      xaxis: {type: "date", showgrid: false, linecolor: line, tickcolor: line,
        rangeselector: {x: 0, y: 1.12, bgcolor: css("--head"), activecolor: line, bordercolor: line,
          font: {color: fg, size: 11}, buttons: [
            {count: 2, label: "2 days", step: "day", stepmode: "backward"},
            {count: 7, label: "a week", step: "day", stepmode: "backward"},
            {step: "all", label: "all"}]}},
      yaxis: {title: {text: "wall time over BLAST+ (below 1: faster)"}, gridcolor: line, linecolor: line,
        zeroline: false, range: [lo - pad, hi + pad]},
      shapes: [{type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: 1, y1: 1, layer: "below",
        line: {color: muted, width: 1.5, dash: "dash"}}].concat(marks.map(function (m) {
          return {type: "line", xref: "x", x0: iso(m[0]), x1: iso(m[0]), yref: "paper", y0: 0, y1: 1, layer: "below",
            line: {color: muted, width: 1, dash: "dot"}};
        })),
      annotations: [{xref: "paper", x: 1, xanchor: "right", yref: "y", y: 1, yanchor: "bottom", text: "BLAST+",
        showarrow: false, font: {size: 11, color: muted}}].concat(marks.map(function (m) {
          return {xref: "x", x: iso(m[0]), xanchor: "left", yref: "paper", y: 1, yanchor: "top", text: " " + m[1],
            showarrow: false, font: {size: 10, color: muted}};
        })),
      hoverlabel: {bgcolor: css("--head"), bordercolor: line, font: {color: fg}, align: "left"},
    };
  }

  function render() {
    if (!el) return;
    if (typeof Plotly === "undefined") {
      el.textContent = "The chart needs Plotly.js from cdn.jsdelivr.net, which did not load.";
      return;
    }
    const built = build();
    if (!built.length) {
      el.textContent = "No runs with rounds yet.";
      return;
    }
    const all = [];
    for (const b of built) all.push.apply(all, b.points);
    Plotly.purge(el);
    Plotly.newPlot(el, built.map(function (b) { return b.trace; }), layout(all),
      {responsive: true, displaylogo: false, modeBarButtonsToRemove: ["select2d", "lasso2d"]});
    el.on("plotly_click", function (e) {
      const run = e.points && e.points[0] && e.points[0].customdata;
      if (run) window.open("runs/" + run + "/", "_blank", "noopener");
    });
  }

  render();
  // Draw again with the other colors when the color scheme changes.
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", render);
})();
"""


def commit_sections(groups):
    """The index's per-commit sections, only the newest open, and the plots' data as JSON for a <script> element."""
    out, data = [], []
    for g in groups:
        if not g["jobs"]:
            continue
        k = len(data)
        data.append({"hash": g["hash"], "subject": subject(g["code"]), "jobs": g["jobs"]})
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
                event = {"schedule": "scheduled", "push": "push", "workflow_dispatch": "manual"}.get(s.get("event", ""), "")
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
<h2>Over the runs</h2>
<p class="muted">Six lines: blastp4j, almostblastp and fastblastp, on each architecture, over the BLAST+ task each one
must match (<code>blastp</code> for the first two, <code>blastp -task blastp-fast</code> for fastblastp, which also
runs more queries). A point is one job: the median of its rounds' ratios, and the bar is the lowest and the highest of
them. A <b>diamond</b> instead of a circle marks a job that ran with <b>512-bit vectors</b>, which
<code>bin/jvm.options</code> gives any CPU with full-width AVX-512 units (AMD Zen 5 and later, Intel parts with
<code>avx512f</code>); the other x86-64 CPUs keep 256, so those points measure code built for a different vector
width. Numerator and denominator come from the same round of the same job, which is what makes the runs on different
runner CPUs comparable; what scatter is left between neighbouring points is the runners themselves. A dotted line marks
a change of database and queries, across which the ratios are of different searches. Click a point for that run's
report; the table below has every run, the Mac workflow's among them, which this chart leaves out.</p>
<div id="history"></div>
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
<script>{PLOT_SCRIPT}</script>
<script>{HISTORY_SCRIPT}</script>"""
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
