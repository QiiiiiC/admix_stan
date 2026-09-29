"""Durable runner: compile once, simulate once, fit in parallel workers, summarise.

Fits proceed in replicate checkpoints (every genome length gets its first N
subsamples before any length gets more), and the summary is rebuilt at each
checkpoint, so partial results are always readable.  Java and the hap-IBD jar
default to genetics_env (see methods.ibd; override with JAVA / HAPIBD_JAR or
HAPIBD_COMMAND).
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
OUT = HERE / "runs" / "default"
sys.path.insert(0, str(HERE.parents[1]))
from methods import ibd
from methods.parallel import Supervisor, now

CHECKPOINTS = (1, 5, 20, 50, 100)


def progress():
    return dict(block_caches=len(list(OUT.glob("pool_*/block_*.npz"))),
                fit_results=len(list(OUT.glob("pool_*/fits/cm_*/rep_*/*/*/*/result.json"))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subsample-limit", type=int, default=None)
    parser.add_argument("--jobs", type=int, default=8, help="parallel workers")
    parser.add_argument("--sources", nargs="+", default=["true", "hapibd"])
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, "-B", str(HERE/"run_study.py"), "plan", "--out", str(OUT),
                    "--sources", *args.sources], check=True)
    c = json.loads((OUT/"manifest.json").read_text())["config"]
    target = c["replicates"] if args.subsample_limit is None else args.subsample_limit
    if not 1 <= target <= c["replicates"]:
        parser.error("subsample limit must be between 1 and the configured replicates")
    run = Supervisor(HERE/"run_study.py", OUT, progress)
    run.save(supervisor_pid=os.getpid(), jobs=args.jobs, subsamples_per_length=target,
             planned_subsamples_per_length=c["replicates"], expected_blocks=c["blocks"]*c["pools"],
             fitted_subsample_limit=0)
    run.stage("compile")
    run.stage("simulate", ["--hapibd-command", " ".join(ibd.hapibd_command()), "--sources", *args.sources],
              jobs=args.jobs, slice_flag="--block-slice")
    for limit in sorted({x for x in CHECKPOINTS if x < target} | {target}):
        run.state["target_subsample_limit"] = limit
        run.stage("fit", ["--replicate-limit", str(limit), "--sources", *args.sources],
                  jobs=args.jobs, slice_flag="--task-slice")
        run.state["fitted_subsample_limit"] = limit
        run.stage("summarize", ["--replicate-limit", str(limit), "--sources", *args.sources])
        run.save(summarized_subsample_limit=limit)
    summary = json.loads((OUT/"summary"/"summary_status.json").read_text())
    run.save(status="completed_with_failed_fits" if summary["failed"] else "completed",
             failed_fits=summary["failed"], finished_utc=now())


if __name__ == "__main__":
    main()
