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
  same CR2 and degrees of freedom. Diagnostics are Poisson-only, as in R.
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

## R front end

`rpkg/` is an R package that calls the same C++ core through Rcpp, so the two
front ends can never disagree numerically. The core sources are not duplicated:
`rpkg/src/core_*.cpp` are one-line stubs that include `cpp/src/*.cpp`, and
`rpkg/src/Makevars` puts `cpp/include` on the include path.

```sh
Rscript -e 'Rcpp::compileAttributes("rpkg")'
R CMD INSTALL rpkg
```

```r
library(spicyglm)
cells <- read.csv("cells.csv")                       # one row per cell
cells$condition <- factor(cells$condition, levels = c("NR", "R"))   # first level is the reference

out <- spicy_glm(cells, condition = "condition", subject = "subject", r = 40)
out$results   # one row per fitted pair
out$skipped   # pairs that could not be fitted, with a reason code

spicy_glm(cells, condition = "condition", family = "binomial", k = 10)
spicy_glm(cells, condition = "condition", r = 40, cr2_method = "naive")
```

Arguments match the Python front end, with `from_` spelled `from` and `n_jobs`
replaced by `n_threads`, which is used to build the neighbour lists. Pairs are
fitted sequentially, so a run takes a few seconds where the threaded Python
front end takes under one; both are far below the R `spicyGLM()` they replace.

Verified against the Python front end on a 1.22M-cell, 185-image dataset:
Poisson and Binomial, at both 10 and 21 cell types, agree to 6e-17 in the log
effect and 7e-16 in the p-value, on all 55 and all 231 pairs.

`compute_diagnostics = TRUE` returns the same four tables as the Python front
end: `out$diagnostics$pair`, `$patient`, `$image` and `$cross_pair`. On the
dataset above the tables agree column by column with Python to 5e-14 relative
(`S_g`) and 1e-16 absolute (the Wilson bounds).

One cosmetic difference: rows of `cross_pair` whose Wilson lower bound is
mathematically zero can come out in a different order, because that bound is a
difference of two equal terms and the two languages round it to 0 and to
1.2e-17 respectively. Sort by a second key if you need the orders to match.


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

Checked on Schürch 2020 with `benchmarks/compare_real.py`, and reproduced by
calling `spicyGLM()` itself. Everything else matches R to 1e-7 (binomial R fits
refitted with a tight convergence tolerance).

- **Reference condition.** R's `buildGLM()` takes the reference level from the
  first image with data for a pair (`dplyr::bind_rows` merges factor levels in
  order of appearance), so on 103 of 419 pairs R silently uses the second
  condition and the log ratio's sign flips. spicyglm always uses the first level.
- **Non-converged binomial fits.** On 4 separated pairs `brglm2::brglmFit` does
  not converge (coefficients near -1e15) and R reports p-values near 1e-15.
  spicyglm returns the Jeffreys-penalised maximum, which matches R's own
  `optimize()` of the same likelihood.
- **Ties at the k-th neighbour** are broken exactly as spatstat's `nnwhich()`
  (a port of R's quicksort and spatstat's search), so counts match R even with
  integer coordinates.

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
