"""Benchmark spicyR's spicyGLM() against spicyglm: fitting time and peak memory.

Usage: python benchmarks/benchmark.py <spicyR> <data_dir> <results.csv>

Synthetic datasets (5,000 cells per image, 8 cell types, 2 images per subject,
condition assigned per subject) are generated into data_dir if missing; a real
dataset is included when data_dir/schurch.csv exists.

Each configuration is one process, run one at a time, repeating the fit inside
the process; the median, min and max fitting times are recorded. Peak memory is
the largest total resident memory of the process and all its children (R's
forked workers when cores > 1), sampled every 0.1 s; shared copy-on-write pages
are counted once per process, so for forked R this overstates real usage. A
run whose total exceeds MEM_LIMIT_GB is stopped and recorded as such. Results
are appended as they finish, and finished configurations are skipped on rerun.
"""

import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

HERE = Path(__file__).resolve().parent
SIZES = [50_000, 100_000, 250_000, 500_000, 1_000_000]  # >= 10 images: at least 2 subjects per condition
CELLS_PER_IMAGE = 5_000
JOBS = [1, 4]
REPEATS = {"spicyglm": 5, "R spicyGLM()": 3}
MEM_LIMIT_GB = 6.0


def synthetic(data_dir, n_cells):
    path = data_dir / f"synthetic_{n_cells}.csv"
    if not path.exists():
        rng = np.random.default_rng(n_cells)
        n_images = n_cells // CELLS_PER_IMAGE
        image = np.repeat(np.arange(n_images), CELLS_PER_IMAGE)
        subject = image // 2
        pd.DataFrame({
            "imageID": [f"img{i:04d}" for i in image],
            "subject": [f"s{s:04d}" for s in subject],
            "condition": np.where(subject < (n_images // 2) // 2, "a", "b"),
            "cellType": rng.choice(list("ABCDEFGH"), n_cells),
            "x": rng.uniform(0, 1000, n_cells),
            "y": rng.uniform(0, 1000, n_cells),
        }).to_csv(path, index=False)
    return path


def run_monitored(cmd):
    """Run cmd; return (per-repeat times, peak total RSS in MB, stopped for memory)."""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    root, peak, stopped = psutil.Process(proc.pid), 0, False
    while proc.poll() is None:
        try:
            tree = [root] + root.children(recursive=True)
            total = sum(p.memory_info().rss for p in tree if p.is_running())
        except psutil.Error:
            total = 0
        peak = max(peak, total)
        if total > MEM_LIMIT_GB * 1024**3:
            for p in root.children(recursive=True) + [root]:
                p.kill()
            stopped = True
            break
        time.sleep(0.1)
    output = proc.communicate()[0]
    if stopped:
        return [], peak / 1024**2, True
    match = re.search(r"elapsed_seconds=([0-9.,]+)", output)
    if proc.returncode != 0 or not match:
        raise RuntimeError(f"benchmark run failed: {' '.join(cmd)}\n{output[-3000:]}")
    return [float(t) for t in match.group(1).split(",")], peak / 1024**2, False


def main():
    spicyR, data_dir, out_csv = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    data_dir.mkdir(parents=True, exist_ok=True)

    runs = []  # (dataset, path, n_cells, family, r_or_k, diagnostics)
    for n in SIZES:
        path = synthetic(data_dir, n)
        runs += [(f"synthetic_{n}", path, n, "poisson", 30, False), (f"synthetic_{n}", path, n, "binomial", 10, False)]
    real = data_dir / "schurch.csv"
    if real.exists():
        n_real = sum(1 for _ in open(real)) - 1
        runs += [("Schurch 2020", real, n_real, "poisson", 50, True), ("Schurch 2020", real, n_real, "binomial", 10, False)]

    rows = pd.read_csv(out_csv).to_dict("records") if out_csv.exists() else []
    done = {(r["dataset"], r["family"], r["implementation"], r["workers"]) for r in rows}
    for dataset, path, n_cells, family, size, diagnostics in runs:
        for impl in ("spicyglm", "R spicyGLM()"):
            for workers in JOBS:
                if (dataset, family, impl, workers) in done:
                    continue
                reps = REPEATS[impl]
                if impl == "spicyglm":
                    cmd = [sys.executable, str(HERE / "run_python.py"), str(path), family, str(size),
                           "1" if diagnostics else "0", str(workers), str(reps)]
                else:
                    cmd = ["Rscript", str(HERE / "run_r_spicyglm.R"), spicyR, str(path), family, str(size),
                           "TRUE" if diagnostics else "FALSE", str(workers), str(reps)]
                times, peak_mb, stopped = run_monitored(cmd)
                rows.append(dict(
                    dataset=dataset, n_cells=n_cells, family=family, diagnostics=diagnostics,
                    implementation=impl, workers=workers, repeats=len(times),
                    seconds_median=np.median(times) if times else np.nan,
                    seconds_min=min(times) if times else np.nan, seconds_max=max(times) if times else np.nan,
                    peak_rss_mb=round(peak_mb), stopped_for_memory=stopped,
                ))
                status = "STOPPED (memory limit)" if stopped else f"median {np.median(times):8.2f} s"
                print(f"{dataset:>18} {family:>8} {impl:>13} workers={workers}: {status}  peak {peak_mb:6.0f} MB",
                      flush=True)
                pd.DataFrame(rows).to_csv(out_csv, index=False)


if __name__ == "__main__":
    main()
