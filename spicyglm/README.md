# spicyglm

spicyGLM (spatial co-localisation inference with a fast CR2 sandwich variance)
with a C++17 numeric core and a Python front end. It is a port of `spicyGLM()`
from spicyR (`gee` branch). Equation numbers in the source refer to
`spicyClub_Supplementary_Math.pdf`.

Status:

- Poisson family (fixed radius): closed-form MLE and Firth, fast CR2 with
  Satterthwaite degrees of freedom, and the leverage / influence /
  point-estimate shift diagnostics with cross-pair flagging.
- Binomial family (fixed-k nearest neighbours): MLE and Firth by a
  one-dimensional root-find per condition (replacing glm / brglm2), with the
  same CR2 and degrees of freedom. The effect is directional, so every ordered
  pair is fitted (see below). Diagnostics are Poisson-only, as in R.
- `cr2_method="naive"` (spicyR's `cr2Method = "naive"`): model-based variance
  ignoring clustering, with a z-test. `cr2Method = "clubSandwich"` is not ported.

## Build

```sh
python3 -m venv .venv
.venv/bin/pip install scikit-build-core pybind11 numpy pandas scipy pytest
.venv/bin/pip install --no-build-isolation -e .
```

Eigen 3.4 is used from the system if CMake finds it, otherwise downloaded.

## Use

```python
import pandas as pd
from spicyglm import spicy_glm

cells = pd.read_csv("cells.csv")  # one row per cell
out = spicy_glm(cells, condition="condition", subject="subject", r=40)
out.results   # one row per fitted pair
out.skipped   # pairs that could not be fitted, with a reason code

out = spicy_glm(cells, condition="condition", subject="subject", r=40, ref="NR")
# ref: the reference condition (default: first category, otherwise sorted first)

out = spicy_glm(cells, condition="condition", subject="subject", r=40, compute_diagnostics=True)
out.diagnostics["pair"]                   # one row per pair: nu, max influence, max shift
out.diagnostics["patient"]                # leverage, influence, shift and within-pair ranks
out.diagnostics["image"]                  # image-level leverage and influence
out.diagnostics["cross_pair"]["patient"]  # flagging rates with Wilson intervals

out = spicy_glm(cells, condition="condition", subject="subject", family="binomial", k=10)
out.results   # log_odds_ratio / odds_ratio instead of log_rate_ratio / rate_ratio

out = spicy_glm(cells, condition="condition", subject="subject", r=40, cr2_method="naive")
```

`n_jobs` fits pairs on several threads (and builds the neighbour lists in
parallel); output is identical to `n_jobs=1`.

### Which pairs are fitted

The binomial effect is directional. A→B models how many of each A cell's k
nearest neighbours are B, and B→A is a different model:

- nearest-neighbour membership is not symmetric, so the directed edge totals
  differ;
- the logit link keeps the background-proportion offset inside
  `expit(beta + logit(p0))`, so the two score equations differ for `beta != 0`
  even when the edge totals agree;
- the reference cells, and with them the CR2 working variances, differ.

On simulated data the two directions gave log odds ratios of −0.789 and −0.451.
The Poisson effect is direction-invariant. Counts within `r` are symmetric, the
offset total `sum_j A_j B_j pi r^2 / |W_j|` is symmetric, and the log link
reduces the estimate and CR2 to those totals. A→B and B→A give the same log
rate ratio, SE, degrees of freedom and p-value (checked against clubSandwich).

| `from_` / `to` | Poisson | Binomial |
|---|---|---|
| omitted | each unordered pair once, plus self-pairs: n(n+1)/2 | every ordered pair, including self-pairs: n² |
| one type each | that pair | that direction |
| lists | unordered pairs among the union | every pair in `from_` × `to` (an omitted side means all types) |

BH is applied across every fitted pair, so a binomial run adjusts over n² tests.
This matches spicyR's `spicyGLM()` from `gee@129b248`.

## Benchmarks

One laptop (Apple silicon, 10 cores, 24 GB), against spicyR's `spicyGLM()` on
the `gee` branch at `b591b17` (loaded with `devtools::load_all`), default
`cr2Method = "fast"`. The later `gee` commits (`2adb696`, `7ff4f90`) only add
the naive option and untrack files; the timed path is unchanged. Fitting time only; medians
over 5 runs (spicyglm) and 3 runs (R); like-for-like workers. Plots and the raw
`timings.csv` are in `benchmarks/results/`.

| Data | Family | R, 1 core | spicyglm, 1 thread | R, 4 cores | spicyglm, 4 threads |
|---|---|---|---|---|---|
| Synthetic, 250k cells | Poisson | 11.3 s | 0.27 s (42x) | 7.3 s | 0.13 s (54x) |
| Synthetic, 250k cells | Binomial | 21.0 s | 0.35 s (61x) | 7.7 s | 0.18 s (43x) |
| Synthetic, 1M cells | Poisson | 66 s | 1.07 s (62x) | stopped at 6 GB | 0.53 s |
| Synthetic, 1M cells | Binomial | 93 s | 1.36 s (68x) | stopped at 6 GB | 0.58 s |
| Schürch 2020, 258k cells | Poisson + diagnostics | 93 s | 4.3 s (22x) | stopped at 6 GB | 4.0 s |
| Schürch 2020, 258k cells | Binomial | 132 s | 0.47 s (279x) | stopped at 6 GB | 0.22 s |

Caveats:

- The binomial rows were timed when both implementations fitted one direction
  per unordered pair. Both now fit every ordered pair, about twice as many fits,
  so absolute binomial times roughly double for both. The speed-ups should hold,
  but they have not been re-timed.
- R with 4 cores was stopped when its total memory (summed over forked
  workers, which overstates shared pages) passed 6 GB, so the 4-worker
  comparison only covers up to 250k cells.
- Peak memory at 1M cells: R 4.4 GB (1 core, including about 2 GB for loading
  spicyR's dependencies) vs spicyglm 0.43 GB.
- Synthetic cells are uniformly scattered; on real tissue the single-thread
  speed-up ranged from 22x (Poisson with diagnostics, where building the pandas
  tables dominates) to 279x (binomial).

Reproduce:

```sh
# optional real dataset
Rscript -e 'library(SpatialDatasets); library(SpatialExperiment)
  spe <- spe_Schurch_2020(); cd <- colData(spe); xy <- spatialCoords(spe)
  write.csv(data.frame(imageID = cd$imageID, subject = paste0("p", cd$patients),
            condition = paste0("group", cd$group), cellType = cd$cellType,
            x = xy[, "x"], y = xy[, "y"]), "<data_dir>/schurch.csv", row.names = FALSE)'
python benchmarks/benchmark.py /path/to/spicyR <data_dir> benchmarks/results/timings.csv
python benchmarks/plot.py benchmarks/results/timings.csv benchmarks/results
python benchmarks/compare_real.py <data_dir>/schurch.csv <r_out_dir>  # outputs vs R
```

## Differences from the R implementation

Checked on Schürch 2020 (258k cells, 140 images, 35 patients, 29 cell types)
with `benchmarks/compare_real.py`. The first comparison, against `gee@b591b17`,
found the two differences below. Both were R bugs and are fixed in later `gee`
commits.

Rerun against `gee@40260c4` (R references from `benchmarks/run_r.R` with
`tight`), every run matches R to 1e-7, with no pair where R used the other
reference and no pair excluded as non-converged:

| Run | Fitted pairs | Skipped pairs |
|---|---|---|
| Poisson, r = 50, with diagnostics (all five tables match) | 419 / 419 | 16 / 16 |
| Binomial, k = 10 (ordered pairs) | 809 / 809 | 32 / 32 |
| Poisson, r = 50, `cr2_method = "naive"` | 429 / 429 | 6 / 6 |
| Binomial, k = 10, `cr2_method = "naive"` | 829 / 829 | 12 / 12 |

`n_jobs = 4` output is identical to `n_jobs = 1` in every run.

- **Reference condition (fixed in `gee@d9f8b6f`).** At `b591b17`, R's
  `buildGLM()` took the reference level from the first image with data for a
  pair (`dplyr::bind_rows` merges factor levels in order of appearance), so on
  103 of 419 pairs R silently used the second condition and the log ratio's sign
  flipped. `d9f8b6f` factors the condition once on the full dataset, so every
  pair uses the first level as spicyglm does. It also adds a `ref=` argument,
  which spicyglm has too.
- **Non-converged binomial fits (fixed in `gee@08d07d9`).** At `b591b17`, on 4
  separated pairs `brglm2::brglmFit` did not converge (coefficients near -1e15)
  and R reported p-values near 1e-15. `08d07d9` damps brglm2's step
  (`slowit = 0.1`), which converges to the same Jeffreys-penalised maximum that
  spicyglm returns by root-finding. The tight refit in
  `tests/r_reference/run_spicyglm.R` uses the same damping; without it, 10 of
  the 809 binomial pairs diverge in the refit even though `gee`'s own fit
  converges.
- **Ties at the k-th neighbour** are broken exactly as spatstat's `nnwhich()`
  (a port of R's quicksort and spatstat's search), so counts match R even with
  integer coordinates.

The fixtures in `tests/r_reference/` are compared with R's `spicyGLM()` code
path. The binomial fixtures were regenerated against `gee@129b248` (ordered
pairs). Regenerating the Poisson and diagnostic fixtures against that commit
reproduces the committed files to 1e-14 relative, so they were left unchanged.

## Layout

- `cpp/include/spicyglm/core.hpp`, `cpp/src/` — window areas, radius counts, k-nearest neighbours, Poisson and binomial fits, CR2, degrees of freedom, per-pair diagnostics
- `cpp/src/bindings.cpp` — pybind11 module `spicyglm._core`
- `python/spicyglm/api.py` — pair enumeration, skip rules, p-values, FDR
- `python/spicyglm/diagnostics.py` — diagnostic tables, percentile ranks, cross-pair flagging
- `tests/test_core.py` — core against dense implementations from the definitions
- `benchmarks/` — timing and memory against R, plots, and the real-data output comparison
- `tests/test_against_r.py`, `tests/test_diagnostics.py` — end-to-end against R; regenerate
  fixtures with `Rscript tests/r_reference/make_fixtures.R /path/to/spicyR` (needs
  spatstat.geom, dplyr, binom and brglm2)
