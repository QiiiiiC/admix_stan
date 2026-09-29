"""Run a study stage as n parallel worker processes over disjoint slices.

A study's runner exposes stages on the command line (`simulate`, `fit`, ...)
that accept `--block-slice k/n` or `--candidate-slice k/n`: worker k takes every
n-th item starting at k.  Each result lives in its own file and is skipped once
written, so workers never share output and every stage resumes after a crash.
"""
from __future__ import annotations

import datetime
import json
from pathlib import Path
import subprocess
import sys
import time


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def in_slice(index, spec):
    """spec 'k/n' (or None): does item `index` belong to worker k of n?"""
    if not spec:
        return True
    k, n = (int(v) for v in spec.split("/"))
    return index % n == k


def prime_arviz():
    """arviz's import writes a once-per-day warning stamp through a SHARED temp
    file in ~/Library/Caches/arviz; eight workers importing it in the same second
    race on the rename and seven die.  Importing it once here, right before the
    workers launch, leaves the stamp current so they never write it."""
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import arviz  # noqa: F401
    except ImportError:
        pass


class Supervisor:
    """Launches stages of `script` against run folder `out`; status goes to out/status.json.

    progress(): optional callable returning extra fields for status.json (e.g.
    counts of finished blocks and fits), polled every `poll` seconds.
    """

    def __init__(self, script, out, progress=None, poll=10):
        self.script, self.out, self.progress, self.poll = Path(script), Path(out), progress, poll
        self.state = dict(status="running", started_utc=now())

    def save(self, **fields):
        self.state.update(fields, updated_utc=now(), **(self.progress() if self.progress else {}))
        tmp = self.out / "status.tmp.json"
        tmp.write_text(json.dumps(self.state, indent=2) + "\n")
        tmp.replace(self.out / "status.json")

    def stage(self, name, extra=(), jobs=1, slice_flag="--candidate-slice"):
        self.state.update(stage=name, stage_started_utc=now())
        if jobs > 1:
            prime_arviz()
        procs = []
        for k in range(jobs):
            args = list(extra) + ([slice_flag, f"{k}/{jobs}"] if jobs > 1 else [])
            log = self.out / (name + (f"_worker{k}" if jobs > 1 else "") + ".log")
            handle = log.open("a", buffering=1)
            handle.write(f"\nStarting {name}: {now()} {args}\n")
            cmd = [sys.executable, "-B", "-u", str(self.script), name, "--out", str(self.out), *args]
            procs.append((subprocess.Popen(cmd, cwd=self.script.parent, stdout=handle,
                                           stderr=subprocess.STDOUT), handle))
        self.state["worker_pids"] = [p.pid for p, _ in procs]
        while any(p.poll() is None for p, _ in procs):
            self.save()
            time.sleep(self.poll)
        codes = [p.returncode for p, _ in procs]
        for _, h in procs:
            h.close()
        self.state["last_returncodes"] = codes
        if any(codes):
            self.save(status="failed", finished_utc=now())
            raise SystemExit(max(codes))
        self.save()
