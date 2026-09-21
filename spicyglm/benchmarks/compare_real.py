"""Check spicy_glm against R outputs written by benchmarks/run_r.R on a real dataset.

Usage: python benchmarks/compare_real.py <cells.csv> <r_out_dir>

Expects R runs (benchmarks/run_r.R with "tight") named <stem>_poisson (r = 50,
diagnostics), <stem>_binomial (k = 10), and optionally <stem>_naive_poisson and
<stem>_naive_binomial (cr2Method = "naive"), where <stem> is the CSV file name
without extension.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
from r_compare import DIAG_TABLES, compare_diagnostic_table, compare_results  # noqa: E402

from spicyglm import spicy_glm  # noqa: E402

cells_csv, r_dir = Path(sys.argv[1]), Path(sys.argv[2])
cells = pd.read_csv(cells_csv)
stem = cells_csv.stem

# name suffix -> spicy_glm arguments; R runs must use the same settings
runs = {
    "poisson": dict(family="poisson", r=50, compute_diagnostics=True),
    "binomial": dict(family="binomial", k=10),
    "naive_poisson": dict(family="poisson", r=50, cr2_method="naive"),
    "naive_binomial": dict(family="binomial", k=10, cr2_method="naive"),
}
for suffix, kw in runs.items():
    family = kw["family"]
    name = f"{stem}_{suffix}"
    if not (r_dir / f"results_{name}.csv").exists():
        print(f"{name}: no R output, skipped")
        continue
    out = spicy_glm(cells, condition="condition", subject="subject", **kw)
    n_pairs, flipped, diverged = compare_results(out, r_dir, name, family)
    print(f"{name}: {n_pairs - len(diverged)} of {n_pairs} fitted pairs and all {len(out.skipped)} skipped pairs "
          f"match R ({len(flipped)} pairs where R used the second condition as reference)")
    if diverged:
        print(f"  excluded {len(diverged)} pairs where R's brglmFit did not converge: {sorted(diverged)}")
    if kw.get("compute_diagnostics"):
        for table in DIAG_TABLES:
            rows = compare_diagnostic_table(out, r_dir, name, table, flipped)
            print(f"  diagnostics '{table}': {rows} rows match R")

    threaded = spicy_glm(cells, condition="condition", subject="subject", n_jobs=4, **kw)
    pd.testing.assert_frame_equal(out.results, threaded.results)
    pd.testing.assert_frame_equal(out.skipped, threaded.skipped)
    print("  n_jobs=4 output identical to n_jobs=1")
