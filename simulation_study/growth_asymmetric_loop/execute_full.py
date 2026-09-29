"""Durable runner: compile once, simulate once, fit subsamples in parallel workers, summarise.

Java and the hap-IBD jar default to the copies inside genetics_env (see
methods.ibd; override with JAVA / HAPIBD_JAR or HAPIBD_COMMAND).  Workers split
the simulated blocks by --block-slice and the 29 candidate orders by
--candidate-slice, so --jobs 8 runs eight processes at once against disjoint files.
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


def progress():
    return dict(block_caches=len(list(OUT.glob("*/pool_*/block_*.npz"))),
                candidate_results=len(list(OUT.glob("*/pool_*/fits/rep_*/*/*/*/result.json"))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subsample-limit", type=int, default=None)
    parser.add_argument("--jobs", type=int, default=8, help="parallel fit workers")
    parser.add_argument("--sources", nargs="+", default=["true", "hapibd"])
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, "-B", str(HERE/"run_study.py"), "plan", "--out", str(OUT)], check=True)
    manifest = json.loads((OUT/"manifest.json").read_text()); c = manifest["config"]
    target = c["replicates"] if args.subsample_limit is None else args.subsample_limit
    if not 1 <= target <= c["replicates"]:
        parser.error("subsample limit must be between 1 and the configured replicates")
    run = Supervisor(HERE/"run_study.py", OUT, progress)
    run.save(supervisor_pid=os.getpid(), jobs=args.jobs,
             subsamples_per_scenario=target, planned_subsamples_per_scenario=c["replicates"],
             expected_blocks=2*c["blocks"]*c["pools"],
             expected_candidate_fits=len(manifest["candidate_orders"])*2*len(args.sources)*2*c["pools"]*target,
             fitted_subsample_limit=0, visualized_subsample_limit=0)
    run.stage("compile")
    run.stage("simulate", ["--hapibd-command", " ".join(ibd.hapibd_command()), "--sources", *args.sources],
              jobs=args.jobs, slice_flag="--block-slice")
    for limit in range(1, target+1):
        run.state["target_subsample_limit"] = limit
        run.stage("fit", ["--replicate-limit", str(limit), "--sources", *args.sources], jobs=args.jobs)
        run.state["fitted_subsample_limit"] = limit
        run.stage("summarize")
        run.save(visualized_subsample_limit=limit)
    summary = json.loads((OUT/"summary"/"summary_status.json").read_text())
    run.save(status="completed_with_failed_fits" if summary["failed"] else "completed",
             failed_candidate_fits=summary["failed"], reliable_candidate_fits=summary["reliable"],
             finished_utc=now())


if __name__ == "__main__":
    main()
