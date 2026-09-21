"""Time spicy_glm on a cells CSV (fitting only; reading the CSV excluded).

Usage: python benchmarks/run_python.py <cells.csv> <family> <r_or_k> <diagnostics 0|1> <n_jobs> <repeats>

Prints one line "elapsed_seconds=<t1>,<t2>,..." with one fitting time per repeat.
"""

import sys
import time

import pandas as pd

from spicyglm import spicy_glm

cells_csv, family, size, diagnostics, n_jobs, repeats = sys.argv[1:7]
cells = pd.read_csv(cells_csv)
kw = dict(r=float(size)) if family == "poisson" else dict(k=int(size))
times = []
for _ in range(int(repeats)):
    start = time.perf_counter()
    spicy_glm(cells, condition="condition", subject="subject", family=family,
              compute_diagnostics=diagnostics == "1", n_jobs=int(n_jobs), **kw)
    times.append(time.perf_counter() - start)
print("elapsed_seconds=" + ",".join(f"{t:.4f}" for t in times))
