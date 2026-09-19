"""Durable runner: simulate once, fit successive subsamples, refresh visualizations."""
import datetime
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
OUT = HERE / "runs" / "default"
JAVA = "/opt/anaconda3/pkgs/openjdk-25.0.2-h258754b_0/lib/jvm/bin/java"
JAR = "/Users/qi/software/hap-ibd/hap-ibd.jar"


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(state):
    state.update(updated_utc=now(), block_caches=len(list(OUT.glob("*/pool_*/block_*.npz"))),
                 candidate_results=len(list(OUT.glob("*/pool_*/fits/rep_*/*/*/*/result.json"))))
    tmp=OUT/"status.tmp.json"; tmp.write_text(json.dumps(state,indent=2)+"\n"); tmp.replace(OUT/"status.json")


def stage(state, name, extra=()):
    cmd=[sys.executable,"-B","-u",str(HERE/"run_study.py"),name,"--out",str(OUT),*extra]
    state.update(stage=name,stage_started_utc=now(),command=cmd)
    with (OUT/(name+".log")).open("a",buffering=1) as handle:
        handle.write(f"\nStarting {name}: {now()} {list(extra)}\n")
        child=subprocess.Popen(cmd,cwd=HERE,stdout=handle,stderr=subprocess.STDOUT)
        state["worker_pid"]=child.pid
        while child.poll() is None:
            save(state); time.sleep(10)
        state["last_returncode"]=child.returncode
    if child.returncode:
        state.update(status="failed",finished_utc=now()); save(state)
        raise SystemExit(child.returncode)
    save(state)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subsample-limit", type=int, default=None)
    parser.add_argument("--wait-for-worker", type=int, default=None,
                        help="Wait for an existing fit worker before resuming cached stages.")
    args = parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    # Make/validate the immutable scientific manifest before launching work.
    subprocess.run([sys.executable,"-B",str(HERE/"run_study.py"),"plan","--out",str(OUT)],check=True)
    manifest=json.loads((OUT/"manifest.json").read_text()); c=manifest["config"]
    target = c["replicates"] if args.subsample_limit is None else args.subsample_limit
    if not 1 <= target <= c["replicates"]:
        parser.error("subsample limit must be between 1 and the configured replicates")
    state=dict(status="running",supervisor_pid=os.getpid(),started_utc=now(),
               subsamples_per_scenario=target,planned_subsamples_per_scenario=c["replicates"],diploid_samples_per_population=c["diploid_samples"],
               expected_blocks=2*c["blocks"]*c["pools"],
               expected_candidate_fits=len(manifest["candidate_orders"])*2*2*2*c["pools"]*target,
               fitted_subsample_limit=0,visualized_subsample_limit=0)
    save(state)
    if args.wait_for_worker:
        state.update(stage="waiting_for_existing_worker", worker_pid=args.wait_for_worker)
        while True:
            try:
                os.kill(args.wait_for_worker, 0)
            except ProcessLookupError:
                break
            save(state)
            time.sleep(10)
    stage(state,"simulate",["--hapibd-command",f"{JAVA} -jar {JAR}"])
    # Already-saved results are skipped on every invocation. Both Ne models,
    # IBD sources and scenarios progress before increasing the subsample limit.
    for limit in range(1,target+1):
        state["target_subsample_limit"]=limit
        stage(state,"fit",["--replicate-limit",str(limit)])
        state["fitted_subsample_limit"]=limit
        stage(state,"summarize")
        state["visualized_subsample_limit"]=limit
        save(state)
    summary=json.loads((OUT/"summary"/"summary_status.json").read_text())
    state.update(status="completed_with_failed_fits" if summary["failed"] else "completed",
                 failed_candidate_fits=summary["failed"],reliable_candidate_fits=summary["reliable"],finished_utc=now())
    save(state)


if __name__ == "__main__":
    main()
