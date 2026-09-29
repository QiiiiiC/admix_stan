"""Run one design through a list of stages with parallel workers, then publish its summary.

    python execute_full.py --stages simulate fit summarize
    python execute_full.py --config ../four_leaf_b_admixture_varying_ne/pipeline_config.json \\
        --out ../four_leaf_b_admixture_varying_ne/runs/pipeline \\
        --publish ../four_leaf_b_admixture_varying_ne/summary/pipeline --stages simulate fit summarize

Workers: simulate and fit use --jobs processes; nuts uses --nuts-jobs processes of 4
chains each.  Every stage resumes from saved files.
"""
import argparse
import os
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from methods.parallel import Supervisor, now


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=HERE / "config.json")
    p.add_argument("--out", type=Path, default=HERE / "runs" / "default")
    p.add_argument("--publish", type=Path, default=HERE / "summary")
    p.add_argument("--stages", nargs="+", required=True,
                   choices=["plan", "compile", "simulate", "fit", "nuts", "summarize"])
    p.add_argument("--jobs", type=int, default=11)
    p.add_argument("--nuts-jobs", type=int, default=3)
    a = p.parse_args()
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    progress = lambda: dict(blocks=len(list(out.glob("pool_000/block_*.npz"))),
                            pathfinder=len(list(out.glob("pool_000/pathfinder/cm_*/rep_*/*/*/*/result.json"))),
                            nuts=len(list(out.glob("pool_000/nuts/cm_*/rep_*/*/*/*/result.json"))))
    run = Supervisor(HERE / "run_study.py", out, progress)
    run.save(supervisor_pid=os.getpid(), stages=a.stages)
    extra = ["--config", str(a.config.resolve())]
    for stage in ["plan", "compile"] + [s for s in a.stages if s not in ("plan", "compile")]:
        if stage == "simulate":
            run.stage(stage, extra, jobs=a.jobs, slice_flag="--block-slice")
        elif stage == "fit":
            run.stage(stage, extra, jobs=a.jobs, slice_flag="--task-slice")
        elif stage == "nuts":
            run.stage(stage, extra, jobs=a.nuts_jobs, slice_flag="--task-slice")
        else:
            run.stage(stage, extra)
        if stage == "summarize":
            if a.publish.exists():
                shutil.rmtree(a.publish)
            shutil.copytree(out / "summary", a.publish)
    run.save(status="completed", finished_utc=now())


if __name__ == "__main__":
    main()
