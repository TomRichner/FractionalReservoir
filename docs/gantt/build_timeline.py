"""Build docs/gantt/plan_timeline.html from docs/gantt/timeline.toml.

    uv run python build_timeline.py                     # "now" = the clock
    uv run python build_timeline.py --now 2026-11-12T10:45-06:00

Standard library only. The page is one self-contained file: inline CSS, JS
and SVG, no external resource, so it opens by double-click from disk.

Scheduling rule. Done tasks sit at their actual times (from git). A running
task starts at its actual start and lasts its estimate, or until a quarter
of an hour past "now" if it has already run longer. Every other task starts
at the later of "now" (rounded up to the quarter hour) and the end of the
last task it depends on, and lasts its estimate. Resources are not
modelled: two tasks whose dependencies allow it are drawn in parallel.

See Plan_Timeline.md for how to update the data.
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
import tomllib
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
Q = timedelta(minutes=15)

CODE_RE = re.compile(r"^(?P<stage>[ABCD])_(?P<step>\d+[a-z]?)_(?P<kind>[PBT])(?P<n>\d*)$")
STATUSES = {"done": "done", "running": "running", "planned": "planned", "tom": "waiting on Tom"}
GLYPH = {"done": "✓", "running": "▶", "planned": "○", "tom": "◆"}
WHO = {
    "lead": "lead agent",
    "fable": "Fable worker",
    "opus": "Opus worker",
    "workers": "worker agents",
    "tom": "Tom",
    "lead+tom": "lead + Tom",
}
KIND = {"P": "Plan", "B": "Build", "T": "Test"}
MAX_NAME = 44
MAX_BLOCK_NAME = 64

# layout, in px
PX = 10            # one quarter hour
GAP_W = 40         # a compressed idle gap
ROW_H = 24
GROUP_H = 22
BLOCK_H = 30
HEAD_H = 46
LEFT_W = 590
RIGHT_PAD = 170


class DataError(Exception):
    pass


@dataclass
class Task:
    id: str
    block: str
    group: str
    name: str
    who: str
    status: str
    what: str
    gate: str
    depends: list[str] = field(default_factory=list)
    estimate: int | None = None
    estimate_of_record: int | None = None
    estimate_of_record_text: str = ""
    estimate_of_record_source: str = ""
    took_min: int | None = None
    rule: str = ""
    guess: bool = False
    rough: bool = False
    start: datetime | None = None
    end: datetime | None = None
    approx_start: bool = False
    worker: str = ""
    issue: str = ""
    source: str = ""
    merge: str = ""
    commits: str = ""
    note: str = ""
    # filled by the scheduler
    s: datetime | None = None
    e: datetime | None = None
    est_end: datetime | None = None
    overran: bool = False


# ------------------------------------------------------------------ time helpers

def floor_q(t: datetime) -> datetime:
    t = t.replace(second=0, microsecond=0)
    return t - timedelta(minutes=t.minute % 15)


def ceil_q(t: datetime) -> datetime:
    f = floor_q(t)
    return f if f == t else f + Q


def fmt_dur(d: timedelta) -> str:
    m = int(round(d.total_seconds() / 60))
    h, m = divmod(m, 60)
    if h and m:
        return f"{h} h {m:02d} min"
    if h:
        return f"{h} h"
    return f"{m} min"


def fmt_elapsed(d: timedelta) -> str:
    m = int(round(d.total_seconds() / 60))
    h, m = divmod(m, 60)
    return f"+{h}:{m:02d}"


def parse_offset(s: str) -> timezone:
    m = re.fullmatch(r"([+-])(\d\d):(\d\d)", s)
    if not m:
        raise DataError(f"meta.utc_offset must look like -06:00, got {s!r}")
    sign = 1 if m.group(1) == "+" else -1
    return timezone(sign * timedelta(hours=int(m.group(2)), minutes=int(m.group(3))))


# ------------------------------------------------------------------ load and validate

def load(path: Path) -> tuple[dict, list[dict], list[Task]]:
    with path.open("rb") as f:
        raw = tomllib.load(f)
    meta = raw.get("meta", {})
    blocks = raw.get("block", [])
    tasks_raw = raw.get("task", [])
    errors: list[str] = []
    block_ids = [b.get("id") for b in blocks]
    for b in blocks:
        if len(b.get("name", "")) > MAX_BLOCK_NAME:
            errors.append(f"block {b.get('id')}: name longer than {MAX_BLOCK_NAME} characters")
    known = set(Task.__dataclass_fields__) - {"s", "e", "est_end", "overran"}
    tasks: list[Task] = []
    seen: set[str] = set()
    for i, t in enumerate(tasks_raw):
        tid = t.get("id", f"<task #{i + 1}>")
        unknown = set(t) - known
        if unknown:
            errors.append(f"{tid}: unknown field(s) {sorted(unknown)}")
        for req in ("id", "block", "group", "name", "who", "status", "what", "gate"):
            if req not in t:
                errors.append(f"{tid}: missing field '{req}'")
        if errors and any(e.startswith(f"{tid}: missing") for e in errors):
            continue
        m = CODE_RE.match(t["id"])
        if not m:
            errors.append(f"{tid}: code does not match <stage>_<step>_<kind>[n] "
                          f"(e.g. A_0c_T, B_3_B, A_4_B2)")
        elif m.group("stage") != t["block"]:
            errors.append(f"{tid}: code's stage letter is not its block {t['block']!r}")
        if tid in seen:
            errors.append(f"{tid}: duplicate task code")
        seen.add(tid)
        if t["block"] not in block_ids:
            errors.append(f"{tid}: unknown block {t['block']!r}")
        if t["status"] not in STATUSES:
            errors.append(f"{tid}: unknown status {t['status']!r} (one of {sorted(STATUSES)})")
        if t["who"] not in WHO:
            errors.append(f"{tid}: unknown who {t['who']!r} (one of {sorted(WHO)})")
        if len(t["name"]) > MAX_NAME:
            errors.append(f"{tid}: name longer than {MAX_NAME} characters, it would be cut")
        for k in ("start", "end"):
            if k in t and not (isinstance(t[k], datetime) and t[k].tzinfo):
                errors.append(f"{tid}: {k} must be a TOML date-time with a UTC offset")
        est = t.get("estimate")
        if est is not None:
            if not isinstance(est, int) or est <= 0 or est % 15:
                errors.append(f"{tid}: estimate {est!r} min is not a positive multiple of 15")
        st = t["status"]
        if st == "done" and ("start" not in t or "end" not in t):
            errors.append(f"{tid}: done but without actual start and end times")
        if st == "done":
            if "rule" not in t:
                errors.append(f"{tid}: done but no 'rule' saying how its duration was measured")
            tm = t.get("took_min")
            if not isinstance(tm, int) or tm <= 0:
                errors.append(f"{tid}: done but no 'took_min' (measured minutes)")
            elif isinstance(t.get("start"), datetime) and isinstance(t.get("end"), datetime):
                meas = (t["end"] - t["start"]).total_seconds() / 60
                if abs(meas - tm) > 1:
                    errors.append(f"{tid}: took_min {tm} does not match end - start "
                                  f"({meas:.0f} min)")
        eor = t.get("estimate_of_record")
        if eor is not None:
            if not isinstance(eor, int) or eor <= 0:
                errors.append(f"{tid}: estimate_of_record must be a positive whole number of minutes")
            if not t.get("estimate_of_record_source"):
                errors.append(f"{tid}: estimate_of_record needs estimate_of_record_source "
                              f"(where the estimate was recorded before the task ran)")
            if st == "planned" or (st == "tom" and "start" not in t):
                errors.append(f"{tid}: estimate_of_record is frozen when a task starts; "
                              f"a task that has not started has only 'estimate'")
        if st == "running" and "start" not in t:
            errors.append(f"{tid}: running but without an actual start time")
        if st != "done" and est is None:
            errors.append(f"{tid}: not done, so it needs an estimate (minutes)")
        if st == "planned" and ("start" in t or "end" in t):
            errors.append(f"{tid}: planned tasks are placed by the scheduler; remove start/end")
        if st == "tom" and "end" in t:
            errors.append(f"{tid}: a task waiting on Tom has no end; mark it done when it ends")
        if "start" in t and "end" in t and isinstance(t["start"], datetime) \
                and isinstance(t["end"], datetime) and t["end"] <= t["start"]:
            errors.append(f"{tid}: end is not after start")
        if not errors or not any(e.startswith(tid) for e in errors):
            tasks.append(Task(**{k: v for k, v in t.items() if k in known}))
    by_id = {t.id: t for t in tasks}
    for t in tasks:
        for d in t.depends:
            if d not in seen:
                errors.append(f"{t.id}: depends on unknown task {d!r}")
            elif d in by_id:
                ds = by_id[d].status
                if t.status == "done" and ds != "done":
                    errors.append(f"{t.id}: done, but depends on {d}, which is {ds}")
                if started(t) and ds != "done":
                    errors.append(f"{t.id}: started, but depends on {d}, which is {ds}")
    if errors:
        raise DataError("\n".join(errors))
    cycle = find_cycle(tasks)
    if cycle:
        raise DataError("dependency cycle: " + " -> ".join(cycle))
    return meta, blocks, tasks


def find_cycle(tasks: list[Task]) -> list[str] | None:
    deps = {t.id: t.depends for t in tasks}
    state: dict[str, int] = {}
    stack: list[str] = []

    def visit(n: str) -> list[str] | None:
        state[n] = 1
        stack.append(n)
        for d in deps.get(n, []):
            if state.get(d) == 1:
                return stack[stack.index(d):] + [d]
            if state.get(d) is None:
                c = visit(d)
                if c:
                    return c
        stack.pop()
        state[n] = 2
        return None

    for t in tasks:
        if state.get(t.id) is None:
            c = visit(t.id)
            if c:
                return c
    return None


# ------------------------------------------------------------------ schedule

def started(t: Task) -> bool:
    """Under way: running, or waiting on Tom with a recorded start."""
    return t.status == "running" or (t.status == "tom" and t.start is not None)


def schedule(tasks: list[Task], now: datetime) -> None:
    by_id = {t.id: t for t in tasks}
    now_q = ceil_q(now)
    done: set[str] = set()

    def place(t: Task) -> None:
        if t.id in done:
            return
        for d in t.depends:
            place(by_id[d])
        if t.status == "done":
            t.s = floor_q(t.start)
            t.e = max(ceil_q(t.end), t.s + Q)
        elif started(t):
            t.s = floor_q(t.start)
            t.est_end = t.s + timedelta(minutes=t.estimate)
            if t.est_end > now:
                t.e = t.est_end
            else:
                t.e = now_q + Q
                t.overran = True
        else:
            t.s = max([now_q] + [by_id[d].e for d in t.depends])
            t.e = t.s + timedelta(minutes=t.estimate)
        done.add(t.id)

    for t in tasks:
        place(t)


def remaining(t: Task, now: datetime) -> timedelta:
    if t.status == "done":
        return timedelta(0)
    if started(t):
        return max(t.e - ceil_q(now), timedelta(0))
    return t.e - t.s


# ------------------------------------------------------------------ axis

class Axis:
    """Maps a time to x. Idle gaps in the past are drawn compressed."""

    def __init__(self, tasks: list[Task], now: datetime, compress_hours: float):
        self.now = now
        self.now_q = ceil_q(now)
        t0 = min(t.s for t in tasks)
        t0 = t0.replace(minute=0)
        t1 = max(max(t.e for t in tasks), self.now_q + Q)
        t1 = ceil_q(t1)
        if t1.minute:
            t1 = t1.replace(minute=0) + timedelta(hours=1)
        self.t0, self.t1 = t0, t1
        # busy intervals, to find gaps in the past with nothing drawn
        iv = sorted((t.s, t.e) for t in tasks)
        merged: list[list[datetime]] = []
        for a, b in iv:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        gaps = []
        for (a0, a1), (b0, _b1) in zip(merged, merged[1:]):
            if b0 - a1 >= timedelta(hours=compress_hours) and b0 <= self.now:
                # keep an hour of context either side
                g0 = a1 + timedelta(hours=1)
                g1 = b0 - timedelta(hours=1)
                g0 = g0.replace(minute=0) + (timedelta(hours=1) if g0.minute else timedelta(0))
                g1 = g1.replace(minute=0)
                if g1 > g0:
                    gaps.append((g0, g1))
        self.gaps = gaps
        # segments: (start, end, x0, compressed)
        segs = []
        x = 0.0
        cur = t0
        for g0, g1 in gaps:
            segs.append((cur, g0, x, False))
            x += (g0 - cur) / Q * PX
            segs.append((g0, g1, x, True))
            x += GAP_W
            cur = g1
        segs.append((cur, t1, x, False))
        x += (t1 - cur) / Q * PX
        self.segs = segs
        self.width = x

    def x(self, t: datetime) -> float:
        for a, b, x0, comp in self.segs:
            if a <= t <= b:
                if comp:
                    return x0 + GAP_W * ((t - a) / (b - a))
                return x0 + (t - a) / Q * PX
        return 0.0 if t < self.t0 else self.width


# ------------------------------------------------------------------ render

def esc(s: str) -> str:
    return html.escape(s, quote=True)


def task_info(t: Task, by_id: dict[str, Task], now: datetime, tz: timezone) -> dict:
    """The data shown in the hover box."""
    m = CODE_RE.match(t.id)
    info = {
        "code": t.id,
        "name": t.name,
        "codeRead": f"stage {m.group('stage')}, step {m.group('step')}, "
                    f"{KIND[m.group('kind')]}" + (f" {m.group('n')}" if m.group("n") else ""),
        "status": STATUSES[t.status],
        "who": WHO[t.who] + (f" ({t.worker})" if t.worker else ""),
        "issue": t.issue,
        "what": t.what,
        "gate": t.gate,
        "waits": [f"{d} {by_id[d].name}" for d in t.depends],
        "merge": t.merge,
        "commits": t.commits,
        "note": t.note,
        "source": t.source,
        "rough": t.rough,
    }
    now_q = ceil_q(now)
    if t.estimate_of_record:
        info["prior"] = (f"{t.estimate_of_record_text or fmt_dur(timedelta(minutes=t.estimate_of_record))}"
                         f" (recorded before it ran: {t.estimate_of_record_source})")
    elif t.status == "done" or started(t):
        info["prior"] = "no prior estimate"
    else:
        info["prior"] = ""
    if t.status == "done":
        a, b = t.start.astimezone(tz), t.end.astimezone(tz)
        took = fmt_dur(timedelta(minutes=t.took_min))
        info["timing"] = (f"Took {took}" + (" (approximate: the start was inferred)"
                                           if t.approx_start else "")
                          + f", {a:%Y-%m-%d %H:%M} to {b:%Y-%m-%d %H:%M}. "
                          f"Measured by: {t.rule}.")
        if t.estimate_of_record:
            r = t.took_min / t.estimate_of_record
            info["timing"] += f" Actual / estimate = {r:.2f}."
    elif started(t):
        a = t.start.astimezone(tz)
        est = fmt_dur(timedelta(minutes=t.estimate))
        info["timing"] = (f"Started {a:%Y-%m-%d %H:%M}"
                          + (" (approximately)" if t.approx_start else "")
                          + f". Estimate: {est}"
                          + (" (a guess)" if t.guess else "") + ".")
        if t.overran:
            info["timing"] += (f" Still running after {fmt_dur(now - a)}: it is taking longer "
                               f"than estimated, and the bar is drawn to a quarter hour past now.")
    else:
        est = fmt_dur(timedelta(minutes=t.estimate))
        info["timing"] = (f"Estimate: {est}" + (" (a guess)" if t.guess else "") +
                          f". Placed at {fmt_elapsed(t.s - now_q)} to {fmt_elapsed(t.e - now_q)} "
                          f"of working time from now, after what it waits for.")
    return info


def build_html(meta: dict, blocks: list[dict], tasks: list[Task], now: datetime) -> str:
    tz = parse_offset(meta.get("utc_offset", "+00:00"))
    now = now.astimezone(tz)
    for t in tasks:
        for k in ("start", "end"):
            v = getattr(t, k)
            if v is not None:
                setattr(t, k, v.astimezone(tz))
        if any(b["id"] == t.block and b.get("rough") for b in blocks):
            t.rough = True
    schedule(tasks, now)
    for t in tasks:
        t.s, t.e = t.s.astimezone(tz), t.e.astimezone(tz)
    axis = Axis(tasks, now, float(meta.get("compress_idle_hours", 3)))
    by_id = {t.id: t for t in tasks}
    now_q = ceil_q(now)

    # rows
    rows: list[tuple[str, object]] = []
    for b in blocks:
        bt = [t for t in tasks if t.block == b["id"]]
        if not bt:
            continue
        rows.append(("block", b))
        group = None
        for t in bt:
            if t.group != group:
                group = t.group
                rows.append(("group", group))
            rows.append(("task", t))
    y = HEAD_H
    placed = []
    for kind, obj in rows:
        h = {"block": BLOCK_H, "group": GROUP_H, "task": ROW_H}[kind]
        placed.append((kind, obj, y, h))
        y += h
    height = y + 8
    width = axis.width + RIGHT_PAD

    # ---------------- chart svg
    c: list[str] = []
    c.append(f'<svg class="chart" width="{width:.0f}" height="{height}" '
             f'viewBox="0 0 {width:.0f} {height}" role="img" '
             f'aria-label="Gantt chart of the final-draft work">')
    c.append(DEFS)
    # projection shading
    xq = axis.x(now_q)
    c.append(f'<rect class="proj" x="{xq:.1f}" y="0" width="{axis.width - xq:.1f}" '
             f'height="{height}"/>')
    # row bands
    for kind, obj, ry, h in placed:
        if kind == "block":
            c.append(f'<rect class="band-block" x="0" y="{ry}" width="{width:.0f}" height="{h}"/>')
        elif kind == "group":
            c.append(f'<rect class="band-group" x="0" y="{ry}" width="{width:.0f}" height="{h}"/>')
    # grid: past clock times, future elapsed time
    ticks = []
    for a, b, x0, comp in axis.segs:
        if comp:
            gx = x0
            c.append(f'<rect class="gap" x="{gx:.1f}" y="18" width="{GAP_W}" height="{height - 18}"/>')
            c.append(f'<text class="gaplabel" x="{gx + GAP_W / 2:.1f}" y="41" '
                     f'text-anchor="middle">idle</text>')
            c.append(f'<text class="gaplabel" x="{gx + GAP_W / 2:.1f}" y="{HEAD_H + 14}" '
                     f'text-anchor="middle">{esc(fmt_dur(b - a))}</text>')
            continue
        t = a
        while t <= b:
            ticks.append(t)
            t += Q
    gap_starts = {g0 for g0, _g1 in axis.gaps}
    seen_x = set()
    day_marks = []
    for t in ticks:
        x = axis.x(t)
        if round(x, 1) in seen_x:
            continue
        seen_x.add(round(x, 1))
        if t < now_q:
            if t.minute == 0:
                c.append(f'<line class="g-hour" x1="{x:.1f}" y1="30" x2="{x:.1f}" y2="{height}"/>')
                if t not in gap_starts:   # its label would sit in the idle column
                    c.append(f'<text class="tick" x="{x + 2:.1f}" y="41">{t:%H}</text>')
                if t.hour == 0 or not day_marks or day_marks[-1][1] != t.date():
                    day_marks.append((x, t.date()))
            else:
                c.append(f'<line class="g-q" x1="{x:.1f}" y1="36" x2="{x:.1f}" y2="{HEAD_H}"/>')
        else:
            el = t - now_q
            mins = int(el.total_seconds() // 60)
            if mins % 60 == 0:
                c.append(f'<line class="g-hour" x1="{x:.1f}" y1="30" x2="{x:.1f}" y2="{height}"/>')
                lab = fmt_elapsed(el)
                if (mins // 60) % 2 == 0 or PX * 4 >= 44:
                    c.append(f'<text class="tick fut" x="{x + 2:.1f}" y="41">{lab}</text>')
            else:
                c.append(f'<line class="g-q" x1="{x:.1f}" y1="36" x2="{x:.1f}" y2="{HEAD_H}"/>')
    for x, d in day_marks:
        c.append(f'<line class="g-day" x1="{x:.1f}" y1="0" x2="{x:.1f}" y2="{height}"/>')
        c.append(f'<text class="day" x="{x + 3:.1f}" y="14">{d:%a %Y-%m-%d}</text>')
    c.append(f'<text class="day fut" x="{xq + 6:.1f}" y="14">projection: working time from now '
             f'(+h:mm), no pauses; not a schedule</text>')
    c.append(f'<line class="axisline" x1="0" y1="{HEAD_H}" x2="{width:.0f}" y2="{HEAD_H}"/>')

    # bars
    for kind, obj, ry, h in placed:
        if kind != "task":
            continue
        t: Task = obj
        x0, x1 = axis.x(t.s), axis.x(t.e)
        cls = f"bar st-{t.status}" + (" rough" if t.rough else "")
        by = ry + 4
        bh = h - 8
        c.append(f'<g class="task" data-id="{esc(t.id)}" tabindex="0" '
                 f'aria-label="{esc(t.id + " " + t.name + ", " + STATUSES[t.status])}">')
        c.append(f'<rect class="{cls}" x="{x0:.1f}" y="{by}" width="{max(x1 - x0, 3):.1f}" '
                 f'height="{bh}" rx="3"/>')
        if started(t) and t.overran and t.est_end:
            xe = axis.x(t.est_end)
            c.append(f'<line class="estmark" x1="{xe:.1f}" y1="{by - 2}" x2="{xe:.1f}" '
                     f'y2="{by + bh + 2}"/>')
        if t.status == "done":
            took = ("≈ " if t.approx_start else "") + fmt_dur(timedelta(minutes=t.took_min))
            if t.estimate_of_record:
                lab = (f"est. {fmt_dur(timedelta(minutes=t.estimate_of_record))}, "
                       f"took {took}")
            else:
                lab = f"took {took}"
        else:
            lab = ("≈ " if t.guess else "") + fmt_dur(timedelta(minutes=t.estimate))
            if t.guess:
                lab += " (guess)"
            if t.overran:
                lab = "running longer than the " + fmt_dur(timedelta(minutes=t.estimate)) + " estimated"
        c.append(f'<text class="barlabel" x="{x1 + 5:.1f}" y="{ry + h / 2 + 4:.1f}">'
                 f'{esc(t.id)} · {esc(lab)}</text>')
        c.append('</g>')
    # now line
    xn = axis.x(now)
    c.append(f'<line class="now" x1="{xn:.1f}" y1="0" x2="{xn:.1f}" y2="{height}"/>')
    c.append(f'<text class="nowlabel" x="{xn - 4:.1f}" y="27" text-anchor="end">now '
             f'{now:%H:%M}</text>')
    c.append('</svg>')

    # ---------------- label svg (sticky)
    L: list[str] = []
    L.append(f'<svg class="labels" width="{LEFT_W}" height="{height}" '
             f'viewBox="0 0 {LEFT_W} {height}" aria-hidden="false">')
    L.append(f'<rect class="lbg" x="0" y="0" width="{LEFT_W}" height="{height}"/>')
    L.append(f'<text class="colhead" x="8" y="41">code</text>'
             f'<text class="colhead" x="92" y="41">task</text>'
             f'<text class="colhead" x="392" y="41">who</text>'
             f'<text class="colhead" x="482" y="41">status</text>')
    L.append(f'<line class="axisline" x1="0" y1="{HEAD_H}" x2="{LEFT_W}" y2="{HEAD_H}"/>')
    for kind, obj, ry, h in placed:
        if kind == "block":
            b = obj
            rem = sum((remaining(t, now) for t in tasks if t.block == b["id"]), timedelta(0))
            L.append(f'<rect class="band-block" x="0" y="{ry}" width="{LEFT_W}" height="{h}"/>')
            L.append(f'<text class="blockname" x="8" y="{ry + 20}">{esc(b["name"])}</text>')
        elif kind == "group":
            L.append(f'<rect class="band-group" x="0" y="{ry}" width="{LEFT_W}" height="{h}"/>')
            L.append(f'<text class="groupname" x="14" y="{ry + 15}">{esc(obj)}</text>')
        else:
            t = obj
            ty = ry + h / 2 + 4
            L.append(f'<g class="task lrow" data-id="{esc(t.id)}">')
            L.append(f'<rect class="lhit" x="0" y="{ry}" width="{LEFT_W}" height="{h}"/>')
            L.append(f'<text class="code" x="8" y="{ty:.1f}">{esc(t.id)}</text>')
            L.append(f'<text class="tname" x="92" y="{ty:.1f}">{esc(t.name)}</text>')
            L.append(f'<text class="who" x="392" y="{ty:.1f}">{esc(WHO[t.who])}</text>')
            L.append(f'<text class="stat st-{t.status}-t" x="482" y="{ty:.1f}">'
                     f'{GLYPH[t.status]} {esc(STATUSES[t.status])}</text>')
            L.append('</g>')
    L.append(f'<line class="divider" x1="{LEFT_W - 0.5}" y1="0" x2="{LEFT_W - 0.5}" y2="{height}"/>')
    L.append('</svg>')

    # ---------------- summary table: remaining estimated work per block
    summ = []
    for b in blocks:
        bt = [t for t in tasks if t.block == b["id"]]
        rem = sum((remaining(t, now) for t in bt), timedelta(0))
        nd = sum(t.status == "done" for t in bt)
        guess = any(t.guess for t in bt if t.status != "done")
        summ.append(
            f'<tr><td>{esc(b["name"])}</td><td class="num">{nd} of {len(bt)}</td>'
            f'<td class="num">{esc(fmt_dur(rem)) if rem else "none"}'
            f'{" (rough, mostly guesses)" if b.get("rough") else (" (includes guesses)" if guess else "")}'
            f'</td></tr>')

    # ---------------- estimates against actuals
    pairs = [t for t in tasks if t.status == "done" and t.estimate_of_record]
    n_done = sum(t.status == "done" for t in tasks)
    if pairs:
        ratios = sorted(t.took_min / t.estimate_of_record for t in pairs)
        k = len(ratios)
        med = ratios[k // 2] if k % 2 else (ratios[k // 2 - 1] + ratios[k // 2]) / 2
        rng = (f"{ratios[0]:.2f}" if k == 1 else f"{ratios[0]:.2f} to {ratios[-1]:.2f}")
        lines = "".join(
            f"<li><code>{esc(t.id)}</code> {esc(t.name)}: estimated "
            f"{esc(t.estimate_of_record_text or fmt_dur(timedelta(minutes=t.estimate_of_record)))}, "
            f"took {esc(fmt_dur(timedelta(minutes=t.took_min)))} "
            f"(ratio {t.took_min / t.estimate_of_record:.2f})</li>" for t in pairs)
        if k >= 3:
            stat = f" Actual / estimate: median {med:.2f}, range {rng}."
        else:   # too few for a median or a range to mean anything
            stat = " Actual / estimate: " + ", ".join(f"{r:.2f}" for r in ratios) + "."
        accuracy = (f"<p>{k} of {n_done} finished tasks had an estimate recorded before they ran."
                    f"{stat}</p><ul>{lines}</ul>")
    else:
        accuracy = (f"<p>None of the {n_done} finished tasks had an estimate recorded before "
                    f"they ran.</p>")
    accuracy += ("<p class=\"meta\">The other finished tasks show their actual time alone. "
                 "Estimates are not filled in after the fact. From now on a task's estimate is "
                 "frozen when it starts, so this comparison grows with every finished task. "
                 "It is information for planning, not a score.</p>")

    info = {t.id: task_info(t, by_id, now, tz) for t in tasks}
    info_json = json.dumps(info, ensure_ascii=False, sort_keys=True, indent=None)
    info_json = info_json.replace("</", "<\\/")

    gen = f"{now:%Y-%m-%d %H:%M} (UTC{meta.get('utc_offset', '+00:00')})"
    title = meta.get("title", "Timeline")
    return PAGE.format(
        title=esc(title),
        generated=esc(gen),
        summary="\n".join(summ),
        accuracy=accuracy,
        labels="\n".join(L),
        chart="\n".join(c),
        info=info_json,
    )


DEFS = """<defs>
<pattern id="p-run" width="8" height="8" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
<rect class="p-run-bg" width="8" height="8"/><line class="p-run-fg" x1="0" y1="0" x2="0" y2="8"/></pattern>
<pattern id="p-tom" width="6" height="6" patternUnits="userSpaceOnUse">
<rect class="p-tom-bg" width="6" height="6"/><circle class="p-tom-fg" cx="3" cy="3" r="1.2"/></pattern>
</defs>"""

PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Final-draft timeline</title>
<style>
:root {{
  --bg: #fbfbfa; --fg: #1d1f22; --muted: #5d636b; --line: #d9dcdf; --line2: #eceef0;
  --band-block: #e6e9ee; --band-group: #f1f3f5; --proj: #f3f6fb; --gap: #eeeeee;
  --done: #7c8794; --run-bg: #f2c46d; --run-fg: #b8862c; --plan-bg: #ffffff; --plan: #3f6fb5;
  --tom-bg: #efe4f7; --tom-fg: #8656a8; --now: #c0392b; --tip-bg: #ffffff; --tip-b: #b9bec5;
  --code: #2e4f7a;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --bg: #17191c; --fg: #e3e5e8; --muted: #9aa1aa; --line: #3a3f45; --line2: #262a2f;
    --band-block: #2a2f36; --band-group: #21252a; --proj: #1b2230; --gap: #24272b;
    --done: #8a96a4; --run-bg: #a8792a; --run-fg: #e0b25f; --plan-bg: #1d232c; --plan: #7aa7e8;
    --tom-bg: #3a2c45; --tom-fg: #c79ae6; --now: #ff6b5b; --tip-bg: #22262b; --tip-b: #4a5058;
    --code: #9cc0f0;
  }}
}}
:root[data-theme="dark"] {{
  --bg: #17191c; --fg: #e3e5e8; --muted: #9aa1aa; --line: #3a3f45; --line2: #262a2f;
  --band-block: #2a2f36; --band-group: #21252a; --proj: #1b2230; --gap: #24272b;
  --done: #8a96a4; --run-bg: #a8792a; --run-fg: #e0b25f; --plan-bg: #1d232c; --plan: #7aa7e8;
  --tom-bg: #3a2c45; --tom-fg: #c79ae6; --now: #ff6b5b; --tip-bg: #22262b; --tip-b: #4a5058;
  --code: #9cc0f0;
}}
html, body {{ margin: 0; background: var(--bg); color: var(--fg);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }}
main {{ padding: 16px; }}
h1 {{ font-size: 20px; margin: 0 0 4px; }}
p {{ margin: 6px 0; max-width: 72rem; }}
.meta {{ color: var(--muted); }}
.cols {{ display: flex; flex-wrap: wrap; gap: 16px 32px; margin: 10px 0 14px; }}
.cols > section {{ flex: 1 1 22rem; max-width: 46rem; }}
h2 {{ font-size: 15px; margin: 6px 0; }}
table {{ border-collapse: collapse; }}
td, th {{ border-bottom: 1px solid var(--line); padding: 3px 8px 3px 0; vertical-align: top; text-align: left; }}
td.num {{ white-space: nowrap; }}
dl.legend {{ display: grid; grid-template-columns: max-content 1fr; gap: 3px 10px; margin: 0; }}
dl.legend dt {{ white-space: nowrap; }}
.sw {{ display: inline-block; width: 26px; height: 12px; vertical-align: -1px; border-radius: 2px; margin-right: 4px; }}
.sw.done {{ background: var(--done); }}
.sw.running {{ background: repeating-linear-gradient(45deg, var(--run-bg) 0 4px, var(--run-fg) 4px 6px); }}
.sw.planned {{ background: var(--plan-bg); border: 1.5px solid var(--plan); box-sizing: border-box; }}
.sw.tom {{ background: radial-gradient(var(--tom-fg) 1px, var(--tom-bg) 1.5px) 0 0 / 6px 6px; border: 1.5px solid var(--tom-fg); box-sizing: border-box; }}
.sw.rough {{ background: var(--plan-bg); border: 1.5px dashed var(--plan); box-sizing: border-box; opacity: .75; }}
.sw.now {{ width: 2px; background: var(--now); }}
code {{ color: var(--code); }}
.scroll {{ overflow-x: auto; border: 1px solid var(--line); background: var(--bg); }}
.inner {{ display: flex; width: max-content; }}
.labels {{ position: sticky; left: 0; z-index: 2; flex: none; }}
.chart {{ flex: none; display: block; }}
svg text {{ fill: var(--fg); font-size: 12px; font-family: system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }}
.lbg {{ fill: var(--bg); }}
.band-block {{ fill: var(--band-block); }}
.band-group {{ fill: var(--band-group); }}
.proj {{ fill: var(--proj); }}
.gap {{ fill: var(--gap); }}
.gaplabel {{ fill: var(--muted); font-size: 10px; }}
.g-hour {{ stroke: var(--line2); stroke-width: 1; }}
.g-q {{ stroke: var(--line); stroke-width: 1; }}
.g-day {{ stroke: var(--line); stroke-width: 1.5; }}
.axisline {{ stroke: var(--line); }}
.divider {{ stroke: var(--line); }}
.tick {{ font-size: 10px; fill: var(--muted); }}
.tick.fut {{ font-style: italic; }}
.day {{ font-size: 11px; font-weight: 600; fill: var(--muted); }}
.day.fut {{ font-style: italic; font-weight: 400; }}
.colhead {{ font-weight: 600; fill: var(--muted); font-size: 11px; }}
.blockname {{ font-weight: 700; font-size: 13px; }}
.groupname {{ font-size: 11px; fill: var(--muted); font-weight: 600; }}
.code {{ font-family: ui-monospace, Consolas, monospace; fill: var(--code); font-size: 12px; }}
.who, .stat {{ font-size: 11px; fill: var(--muted); }}
.lhit {{ fill: transparent; }}
.lrow:hover .lhit, .lrow.hl .lhit {{ fill: var(--proj); }}
.bar {{ stroke-width: 1.5; }}
.st-done {{ fill: var(--done); stroke: var(--done); }}
.st-running {{ fill: url(#p-run); stroke: var(--run-fg); }}
.st-planned {{ fill: var(--plan-bg); stroke: var(--plan); }}
.st-tom {{ fill: url(#p-tom); stroke: var(--tom-fg); }}
.rough {{ stroke-dasharray: 5 3; opacity: .78; }}
.p-run-bg {{ fill: var(--run-bg); }} .p-run-fg {{ stroke: var(--run-fg); stroke-width: 3; }}
.p-tom-bg {{ fill: var(--tom-bg); }} .p-tom-fg {{ fill: var(--tom-fg); }}
.estmark {{ stroke: var(--fg); stroke-width: 1.5; stroke-dasharray: 2 2; }}
.barlabel {{ font-size: 10.5px; fill: var(--muted); paint-order: stroke; stroke: var(--bg); stroke-width: 3px; stroke-linejoin: round; }}
.now {{ stroke: var(--now); stroke-width: 2; }}
.nowlabel {{ fill: var(--now); font-size: 11px; font-weight: 700; }}
g.task {{ cursor: default; outline: none; }}
g.task:focus .bar, g.task:hover .bar, g.task.hl .bar {{ stroke: var(--fg); stroke-width: 2.5; }}
#tip {{ position: fixed; z-index: 10; display: none; max-width: 30rem; background: var(--tip-bg);
  color: var(--fg); border: 1px solid var(--tip-b); border-radius: 6px; padding: 10px 12px;
  box-shadow: 0 4px 18px rgba(0,0,0,.18); font-size: 13px; pointer-events: none; }}
#tip h3 {{ margin: 0 0 2px; font-size: 14px; }}
#tip .read {{ color: var(--muted); font-size: 11.5px; margin-bottom: 6px; }}
#tip p {{ margin: 4px 0; }}
#tip .k {{ color: var(--muted); }}
#tip .roughnote {{ font-style: italic; color: var(--muted); }}
</style>
</head>
<body>
<main>
<h1>{title}</h1>
<p class="meta">Generated {generated} from <code>docs/gantt/timeline.toml</code> by <code>build_timeline.py</code>.</p>
<p><strong>For visibility only.</strong> This chart shows the current plan and an estimate of how long its steps take. It sets no finish date for anything. Past work sits at its real times from the git log; everything right of the red <em>now</em> line is a projection in working time from now, as if the work ran without pause (no nights, no waiting for Tom). It is regenerated after every merge and the estimates are corrected each time; blocks B, C and D will move most as they are planned.</p>
<p><strong>For agents building or testing:</strong> nobody is bound by this chart. A step takes as long as doing it properly takes; never shorten, skip, reorder or hurry work to match a bar. A step that runs longer than its estimate is not a problem, only an estimate to correct at the next update. Quality comes first and the plan is the guide.</p>
<div class="cols">
<section>
<h2>Remaining estimated work per block</h2>
<table><thead><tr><th>Block</th><th>Tasks done</th><th>Remaining, working time</th></tr></thead>
<tbody>
{summary}
</tbody></table>
<h2>Estimates against actuals</h2>
{accuracy}
</section>
<section>
<h2>How to read it</h2>
<dl class="legend">
<dt><code>A_0c_T</code></dt><dd>Task code <em>stage_step_kind</em>: block A, step 0c of the plan, a Test. Kinds: P plan, B build, T test. A trailing number (<code>A_4_B2</code>) separates several tasks of one kind in one step. Rough blocks number their steps 0, 1, 2&hellip; with 0 the planning step.</dd>
<dt><span class="sw done"></span>&#x2713; done</dt><dd>finished; drawn at its actual times</dd>
<dt><span class="sw running"></span>&#x25B6; running</dt><dd>under way now; a dotted tick marks where its estimate ended if it runs longer</dd>
<dt><span class="sw planned"></span>&#x25CB; planned</dt><dd>not started; placed after what it waits for</dd>
<dt><span class="sw tom"></span>&#x25C6; waiting on Tom</dt><dd>needs Tom: a decision, a review, or a run only Tom starts</dd>
<dt><span class="sw rough"></span>dashed</dt><dd>rough: not planned in detail, duration a guess (&asymp;)</dd>
<dt><span class="sw now"></span> now</dt><dd>when this page was generated</dd>
<dt>grey &ldquo;idle&rdquo; column</dt><dd>a past stretch with no recorded work, drawn narrow</dd>
</dl>
<p>Hover over (or tab to) a bar or a row for what the step does, who does it, what it waits for and how it is checked.</p>
</section>
</div>
<div class="scroll"><div class="inner">
{labels}
{chart}
</div></div>
</main>
<div id="tip" role="tooltip"></div>
<script id="taskdata" type="application/json">{info}</script>
<script>
(function () {{
  var info = JSON.parse(document.getElementById('taskdata').textContent);
  var tip = document.getElementById('tip');
  function el(tag, cls, text) {{
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined) e.textContent = text;
    return e;
  }}
  function line(label, text) {{
    var p = el('p'); p.appendChild(el('span', 'k', label + ': ')); p.appendChild(document.createTextNode(text));
    return p;
  }}
  function fill(d) {{
    tip.textContent = '';
    tip.appendChild(el('h3', null, d.code + '  ' + d.name));
    tip.appendChild(el('div', 'read', d.codeRead + (d.issue ? ' \\u00B7 ' + d.issue : '')));
    tip.appendChild(el('p', null, d.what));
    tip.appendChild(line('Who', d.who));
    tip.appendChild(line('Status', d.status));
    tip.appendChild(line('Time', d.timing));
    if (d.prior) tip.appendChild(line('Prior estimate', d.prior));
    tip.appendChild(line('Waits for', d.waits.length ? d.waits.join('; ') : 'nothing'));
    tip.appendChild(line('Checked by', d.gate));
    if (d.merge) tip.appendChild(line('Merged as', d.merge));
    if (d.commits) tip.appendChild(line('Commits', d.commits));
    if (d.note) tip.appendChild(line('Note', d.note));
    if (d.source) tip.appendChild(line('Times from', d.source));
    if (d.rough) tip.appendChild(el('p', 'roughnote', 'Rough: not yet planned in detail; it will change as it is planned.'));
  }}
  function place(x, y) {{
    var w = tip.offsetWidth, h = tip.offsetHeight, vw = window.innerWidth, vh = window.innerHeight;
    var left = x + 16, top = y + 16;
    if (left + w > vw - 8) left = Math.max(8, x - w - 16);
    if (top + h > vh - 8) top = Math.max(8, vh - h - 8);
    tip.style.left = left + 'px'; tip.style.top = top + 'px';
  }}
  function mark(id, on) {{
    document.querySelectorAll('g.task[data-id="' + id + '"]').forEach(function (g) {{ g.classList.toggle('hl', on); }});
  }}
  var current = null;
  document.querySelectorAll('g.task').forEach(function (g) {{
    var id = g.getAttribute('data-id');
    g.addEventListener('mouseenter', function (ev) {{ current = id; fill(info[id]); tip.style.display = 'block'; place(ev.clientX, ev.clientY); mark(id, true); }});
    g.addEventListener('mousemove', function (ev) {{ place(ev.clientX, ev.clientY); }});
    g.addEventListener('mouseleave', function () {{ tip.style.display = 'none'; mark(id, false); current = null; }});
    g.addEventListener('focus', function () {{ var r = g.getBoundingClientRect(); fill(info[id]); tip.style.display = 'block'; place(r.right, r.top); mark(id, true); }});
    g.addEventListener('blur', function () {{ tip.style.display = 'none'; mark(id, false); }});
  }});
  document.addEventListener('keydown', function (ev) {{ if (ev.key === 'Escape') tip.style.display = 'none'; }});
  // open with "now" in view: a quarter of the visible chart lies before it
  var scroller = document.querySelector('.scroll');
  var nowLine = document.querySelector('line.now');
  var leftW = document.querySelector('svg.labels').getBoundingClientRect().width;
  function scrollTo(x) {{ scroller.scrollLeft = Math.max(0, x - (scroller.clientWidth - leftW) * 0.25); }}
  scrollTo(parseFloat(nowLine.getAttribute('x1')));
  // #show=CODE opens that task's hover box (for checking the page without a mouse)
  var m = /^#show=([A-Za-z0-9_]+)$/.exec(location.hash);
  if (m && info[m[1]]) {{
    var g = document.querySelector('svg.chart g.task[data-id="' + m[1] + '"]');
    var bar = g.querySelector('rect');
    scrollTo(parseFloat(bar.getAttribute('x')));
    var r = bar.getBoundingClientRect();
    fill(info[m[1]]); tip.style.display = 'block'; place(r.left, r.bottom); mark(m[1], true);
  }}
}})();
</script>
</body>
</html>
"""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, default=HERE / "timeline.toml")
    ap.add_argument("--out", type=Path, default=HERE / "plan_timeline.html")
    ap.add_argument("--now", help="ISO time with offset, e.g. 2026-11-12T10:45-06:00 "
                                  "(default: the clock)")
    args = ap.parse_args(argv)
    try:
        meta, blocks, tasks = load(args.data)
        tz = parse_offset(meta.get("utc_offset", "+00:00"))
        if args.now:
            now = datetime.fromisoformat(args.now)
            if now.tzinfo is None:
                raise DataError("--now needs a UTC offset, e.g. 2026-11-12T10:45-06:00")
        else:
            now = datetime.now(tz).replace(second=0, microsecond=0)
        page = build_html(meta, blocks, tasks, now)
    except DataError as e:
        print("timeline data error:\n" + str(e), file=sys.stderr)
        return 1
    args.out.write_text(page, encoding="utf-8", newline="\n")
    nd = sum(t.status == "done" for t in tasks)
    print(f"wrote {args.out.name}: {len(tasks)} tasks ({nd} done), now = {now:%Y-%m-%d %H:%M}")
    for b in blocks:
        rem = sum((remaining(t, now.astimezone(tz)) for t in tasks if t.block == b["id"]),
                  timedelta(0))
        print(f"  block {b['id']}: remaining estimated work {fmt_dur(rem) if rem else 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
