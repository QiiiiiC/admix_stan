"""Durable runner: compile once, simulate once, fit subsamples in parallel workers, summarise.

Java and the hap-IBD jar default to the copies inside genetics_env; override with
the JAVA / HAPIBD_JAR environment variables.  Workers split the 100 simulated
blocks by --block-slice and the 29 candidate orders by --candidate-slice, so
--jobs 8 runs eight processes at once against disjoint files.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
OUT = HERE / "runs" / "default"
JAVA = os.environ.get("JAVA", "/opt/miniconda3/envs/genetics_env/bin/java")
JAR = os.environ.get("HAPIBD_JAR", "/opt/miniconda3/envs/genetics_env/share/hap-ibd-1.0.rev20May22.818-0/hap-ibd.jar")


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(state):
    state.update(updated_utc=now(), block_caches=len(list(OUT.glob("*/pool_*/block_*.npz"))),
                 candidate_results=len(list(OUT.glob("*/pool_*/fits/rep_*/*/*/*/result.json"))))
    tmp = OUT/"status.tmp.json"; tmp.write_text(json.dumps(state, indent=2)+"\n"); tmp.replace(OUT/"status.json")


def launch(name, extra, log):
    cmd = [sys.executable, "-B", "-u", str(HERE/"run_study.py"), name, "--out", str(OUT), *extra]
    handle = log.open("a", buffering=1)
    handle.write(f"\nStarting {name}: {now()} {list(extra)}\n")
    return subprocess.Popen(cmd, cwd=HERE, stdout=handle, stderr=subprocess.STDOUT), handle


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


def stage(state, name, extra=(), jobs=1):
    state.update(stage=name, stage_started_utc=now())
    if jobs > 1:
        prime_arviz()
    procs = []
    flag = "--block-slice" if name == "simulate" else "--candidate-slice"
    for k in range(jobs):
        args = list(extra) + ([flag, f"{k}/{jobs}"] if jobs > 1 else [])
        procs.append(launch(name, args, OUT/(name+(f"_worker{k}" if jobs > 1 else "")+".log")))
    state["worker_pids"] = [p.pid for p, _ in procs]
    while any(p.poll() is None for p, _ in procs):
        save(state); time.sleep(10)
    codes = [p.returncode for p, _ in procs]
    for _, h in procs: h.close()
    state["last_returncodes"] = codes
    if any(codes):
        state.update(status="failed", finished_utc=now()); save(state)
        raise SystemExit(max(codes))
    save(state)


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
    state = dict(status="running", supervisor_pid=os.getpid(), started_utc=now(), jobs=args.jobs,
                 subsamples_per_scenario=target, planned_subsamples_per_scenario=c["replicates"],
                 expected_blocks=2*c["blocks"]*c["pools"],
                 expected_candidate_fits=len(manifest["candidate_orders"])*2*len(args.sources)*2*c["pools"]*target,
                 fitted_subsample_limit=0, visualized_subsample_limit=0)
    save(state)
    stage(state, "compile")
    stage(state, "simulate", ["--hapibd-command", f"{JAVA} -jar {JAR}", "--sources", *args.sources], jobs=args.jobs)
    for limit in range(1, target+1):
        state["target_subsample_limit"] = limit
        stage(state, "fit", ["--replicate-limit", str(limit), "--sources", *args.sources], jobs=args.jobs)
        state["fitted_subsample_limit"] = limit
        stage(state, "summarize")
        state["visualized_subsample_limit"] = limit
        save(state)
    summary = json.loads((OUT/"summary"/"summary_status.json").read_text())
    state.update(status="completed_with_failed_fits" if summary["failed"] else "completed",
                 failed_candidate_fits=summary["failed"], reliable_candidate_fits=summary["reliable"], finished_utc=now())
    save(state)


if __name__ == "__main__":
    main()
