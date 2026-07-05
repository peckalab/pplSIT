#!/usr/bin/env python3
"""Summarize Snakemake runtime overlap from .out logs and benchmark TSVs."""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import mean


STAMP_RE = re.compile(r"^\[(?P<stamp>[^\]]+)\]$")
RULE_RE = re.compile(r"^rule (?P<rule>[A-Za-z0-9_]+):")
JOBID_RE = re.compile(r"^\s*jobid:\s*(?P<jobid>\d+)\s*$")
FINISH_RE = re.compile(r"^Finished job (?P<jobid>\d+)\.")


@dataclass
class Job:
    source: Path
    jobid: str
    rule: str
    start: datetime
    finish: datetime | None = None

    @property
    def seconds(self) -> float | None:
        if self.finish is None:
            return None
        return (self.finish - self.start).total_seconds()


def parse_stamp(text: str) -> datetime:
    return datetime.strptime(text, "%a %b %d %H:%M:%S %Y")


def parse_log(path: Path) -> list[Job]:
    jobs: dict[str, Job] = {}
    pending_stamp: datetime | None = None
    pending_rule: tuple[str, datetime] | None = None

    for line in path.read_text(errors="replace").splitlines():
        stamp_match = STAMP_RE.match(line)
        if stamp_match:
            pending_stamp = parse_stamp(stamp_match.group("stamp"))
            continue

        finish_match = FINISH_RE.match(line)
        if finish_match and pending_stamp is not None:
            job = jobs.get(finish_match.group("jobid"))
            if job is not None:
                job.finish = pending_stamp
            continue

        rule_match = RULE_RE.match(line)
        if rule_match and pending_stamp is not None:
            pending_rule = (rule_match.group("rule"), pending_stamp)
            continue

        jobid_match = JOBID_RE.match(line)
        if jobid_match and pending_rule is not None:
            rule, start = pending_rule
            jobid = jobid_match.group("jobid")
            jobs[jobid] = Job(source=path, jobid=jobid, rule=rule, start=start)
            pending_rule = None

    return list(jobs.values())


def overlap_seconds(a: Job, b: Job) -> float:
    if a.finish is None or b.finish is None:
        return 0.0
    start = max(a.start, b.start)
    finish = min(a.finish, b.finish)
    return max(0.0, (finish - start).total_seconds())


def fmt_seconds(seconds: float | None) -> str:
    if seconds is None:
        return "unfinished"
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h:d}:{m:02d}:{s:02d}"
    return f"{m:d}:{s:02d}"


def rule_summary(jobs: list[Job]) -> list[tuple[str, int, float, float]]:
    by_rule: dict[str, list[float]] = defaultdict(list)
    for job in jobs:
        if job.seconds is not None:
            by_rule[job.rule].append(job.seconds)
    rows = []
    for rule, values in by_rule.items():
        rows.append((rule, len(values), mean(values), max(values)))
    return sorted(rows, key=lambda row: row[3], reverse=True)


def overlap_summary(jobs: list[Job]) -> list[tuple[str, str, float, int]]:
    pairs: dict[tuple[str, str], float] = defaultdict(float)
    counts: Counter[tuple[str, str]] = Counter()
    finished = [job for job in jobs if job.finish is not None]
    for i, left in enumerate(finished):
        for right in finished[i + 1 :]:
            overlap = overlap_seconds(left, right)
            if overlap <= 0:
                continue
            pair = tuple(sorted((left.rule, right.rule)))
            pairs[pair] += overlap
            counts[pair] += 1
    rows = [(a, b, sec, counts[(a, b)]) for (a, b), sec in pairs.items()]
    return sorted(rows, key=lambda row: row[2], reverse=True)


def read_benchmarks(root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in root.rglob("*.tsv"):
        try:
            with path.open(newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    row["benchmark"] = str(path)
                    rows.append(row)
        except (OSError, csv.Error):
            continue
    return rows


def benchmark_rule(path: str) -> str:
    name = Path(path).name
    if name.endswith(".tsv"):
        name = name[:-4]
    return name.split(".")[0]


def print_markdown(jobs: list[Job], benchmark_rows: list[dict[str, str]], top: int) -> None:
    complete = [job for job in jobs if job.finish is not None]
    unfinished = [job for job in jobs if job.finish is None]

    print("# Parallelism Diagnostic")
    print()
    print(f"Parsed jobs: {len(jobs)} ({len(complete)} complete, {len(unfinished)} unfinished)")
    if complete:
        start = min(job.start for job in complete)
        finish = max(job.finish for job in complete if job.finish is not None)
        print(f"Observed wall time: {fmt_seconds((finish - start).total_seconds())}")
    print()

    print("## Slowest Rules From Logs")
    print()
    print("| rule | jobs | mean | max |")
    print("| --- | ---: | ---: | ---: |")
    for rule, count, avg, max_sec in rule_summary(jobs)[:top]:
        print(f"| {rule} | {count} | {fmt_seconds(avg)} | {fmt_seconds(max_sec)} |")
    print()

    print("## Largest Rule Overlaps")
    print()
    print("| rule A | rule B | total overlap | overlaps |")
    print("| --- | --- | ---: | ---: |")
    for left, right, seconds, count in overlap_summary(jobs)[:top]:
        print(f"| {left} | {right} | {fmt_seconds(seconds)} | {count} |")
    print()

    if benchmark_rows:
        grouped: dict[str, list[float]] = defaultdict(list)
        for row in benchmark_rows:
            seconds = row.get("s")
            if not seconds:
                continue
            try:
                grouped[benchmark_rule(row["benchmark"])].append(float(seconds))
            except ValueError:
                continue

        print("## Benchmark TSVs")
        print()
        print("| benchmark rule | runs | mean | max |")
        print("| --- | ---: | ---: | ---: |")
        for rule, values in sorted(grouped.items(), key=lambda item: max(item[1]), reverse=True)[:top]:
            print(f"| {rule} | {len(values)} | {fmt_seconds(mean(values))} | {fmt_seconds(max(values))} |")
        print()

    print("## Suggested A/B Runs")
    print()
    print("Run these from the `workflow/` directory.")
    print()
    print("Current behavior:")
    print("`snakemake --configfile ../config/miguel/017388.yaml --use-conda --cores 64`")
    print()
    print("Serialize GPU and heavy I/O suspects:")
    print("`snakemake --configfile ../config/miguel/017388.yaml --use-conda --cores 64 --resources gpu=1 heavy_io=1 big_cpu=1`")
    print()
    print("Test smaller PSTH bootstrap thread claims:")
    print("`snakemake --configfile ../config/miguel/017388.yaml --use-conda --cores 64 --set-threads psth_bootstrap_profiles=16 psth_bootstrap_plots=16`")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", default=["workflow"], help="Log files or directories to scan.")
    parser.add_argument("--top", type=int, default=20, help="Number of rows per section.")
    args = parser.parse_args()

    log_paths: list[Path] = []
    benchmark_roots: list[Path] = []
    for raw in args.paths:
        path = Path(raw)
        if path.is_dir():
            log_paths.extend(path.rglob("*.out"))
            benchmark_roots.append(path)
        elif path.is_file():
            log_paths.append(path)
            benchmark_roots.append(path.parent)

    jobs: list[Job] = []
    for path in sorted(set(log_paths)):
        jobs.extend(parse_log(path))

    benchmark_rows: list[dict[str, str]] = []
    for root in sorted(set(benchmark_roots)):
        benchmark_rows.extend(read_benchmarks(root))

    print_markdown(jobs, benchmark_rows, args.top)


if __name__ == "__main__":
    main()
