# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

spicyR is a Bioconductor R package for spatial analysis of in situ cytometry / spatial omics data. It performs
inference on changes in spatial relationships between pairs of cell types across images/samples. The workflow has
three stages:

1. **Summarize** spatial localization between pairs of cell types for each image (`getPairwise()`, based on
   inhomogeneous L-function estimation via `spatstat`).
2. **Model** the variability in those localization summary statistics as a function of cell counts (weighting).
3. **Test** for association between spatial localization and a response variable, via one of several backends:
   mixed-effects models (`lmer`/`coxme`), GLM (logistic-style, `spicyGLM`), or GEE (`spicyGEE`).

This repo (`gee` branch) is being developed as part of a project (University of Sydney, SydneyBioX group) developing spicyGLM. 
In it we primairly now are comparing **spicyGLM** (the novel method developed here) against **spicyR** (the established L-function baseline).
The core methodological contribution is a fast closed-form CR2 sandwich variance estimator: the block structure of
the per-patient design collapses the `N_i x N_i` eigendecomposition to a small `n_i x n_i` reduced problem (`n_i` =
images per patient), giving ~39x speedup over the `clubSandwich` baseline. The reduced problem is
diagonal-plus-rank-1 (DPR1), but the actual speedup on the default path comes from solving it with a **direct
generic `eigen()` call**, not the bespoke `dpr1eig` secular-equation solver: empirically `eigen()` is 6–56x faster
at spicyGLM's scale (`n_i` in the single digits, rarely > 20), where R's interpreted per-eigenpair loop in the
shift-and-invert route dominates. Hence `fastMethod = "direct"` is the default; `dpr1eig` (`fastMethod = "dpr1"`) is
implemented and validated as a reference alternative but is not the fast path.

**Scope note:** this repo is for package development only. All analysis work (simulations, applying the pipeline
to real datasets, generating thesis figures/results) happens separately on the HPC cluster (albona), not here.
Don't assume access to analysis scripts, `.rds` caches, or HPC output from this repo.

## Common commands

Development is done via the R console/RStudio, not a CLI build tool. Typical commands (run from an R session with
`devtools` loaded, working directory = repo root):

```r
devtools::load_all(".")        # load package for interactive development
devtools::document()           # regenerate NAMESPACE and man/*.Rd from roxygen comments
devtools::test()               # run the full testthat suite
devtools::test(filter = "utilities")  # run tests in tests/testthat/test-utilities.R only
testthat::test_file("tests/testthat/test-spicyR.R")  # run a single test file directly
devtools::check()              # full R CMD check (build, install, test, docs)
BiocCheck::BiocCheck(".")      # Bioconductor-specific package checks
```

- Tests use testthat edition 3 (`Config/testthat/edition: 3` in DESCRIPTION).
- After changing any roxygen `#'` doc comments or function signatures/exports, run `devtools::document()` before
  committing — `NAMESPACE` and `man/*.Rd` are generated files and must stay in sync with the source.
- The regression test in `tests/testthat/test-spicyR.R` compares `spicy()` output against a stored snapshot
  (`inst/testdata/original_result.rds`) with `tolerance = 0.01`. Intentional changes to `spicy()`'s numeric output
  require regenerating that snapshot.

## Architecture

### Entry points and dispatch

`spicy()` (`R/spicy.R`) is the main user-facing function. It accepts a `SingleCellExperiment`/`SpatialExperiment`/
data frame, standardizes it via `.format_data()` (`R/utilities.R`), computes pairwise spatial association statistics
via `getPairwise()`, computes per-image weights via `getWeightFunction()`, then dispatches to a modelling backend
based on the type of `condition`:
- A `Surv` object → `spatialSurv()` (survival / mixed-effects Cox model).
- Otherwise → linear mixed-effects model (`lmer`) or, when `subject` has a 1:1 mapping to `imageID`, a plain `lm`.

`spicyGEE()` (`R/spicyGEE.R`) and the GLM path (`R/spicyGLM.R`, `R/spicyGLM_meanModel.R`) are alternative modelling
backends that reuse the same pairwise-association/weighting machinery but fit GEE or GLM-style models per cell-type
pair instead of mixed-effects models. Many arguments (`imageID`, `cellType`, `spatialCoords`, `from`, `to`, `r`,
`window`, `cores`) are shared across `spicy()`, `spicyGEE()`, and `getPairwise()` — keep naming/semantics consistent
when touching any of them.

Several arguments across `spicy()`/`getPairwise()` are deprecated in favor of renamed equivalents (e.g. `nCores` →
`cores`, `BPPARAM` → `cores`, `imageIDCol` → `imageID`, `cellTypeCol` → `cellType`, `spatialCoordCols` →
`spatialCoords`, `Rs` → `r`) — see `NEWS.txt` for the full mapping. New code should use the modern names; deprecated
aliases are kept for backward compatibility via `lifecycle::deprecate_soft()`/`deprecate_warn()` and must continue to
work.

### GLM backend internals (`R/spicyGLM.R`, `R/spicyGLM_meanModel.R`, `R/cr2_fast_multi.R`, `R/dpr1eig.R`)

This is the most numerically involved part of the package:
- `R/spicyGLM_meanModel.R` fits a per-cell-type-pair GLM (`spicyFit` S3 class, with `coef.spicyFit`,
  `fitted.spicyFit`, `residuals.spicyFit` methods), supporting MLE and Firth-corrected estimators via closed-form
  or `brglm2`/base `glm` fallbacks (`fit_pair()` dispatches based on data rank/separation, see
  `check_rank1_design()`).
- `R/cr2_fast_multi.R` implements a fast CR2 (cluster-robust sandwich, bias-reduced) variance estimator and Wald
  test (`vcovCR2_fast_multi()`, `waldTest_CR2_fast()`) as a hand-rolled alternative to calling `clubSandwich::vcovCR`
  directly, for performance on many cell-type pairs. The reduced per-patient eigendecomposition is done by a direct
  `eigen()` call by default (`method = "direct"`, threaded from `fastMethod`); `method = "dpr1"` routes through
  `dpr1eig.R` instead and is a validated reference only — `"direct"` is empirically faster at spicyGLM's scale.
- `R/dpr1eig.R` is a from-scratch implementation of the diagonal-plus-rank-1 secular equation eigensolver
  (`dpr1eig()`). It is a validated *reference* alternative for the reduced eigendecomposition in `cr2_fast_multi.R`
  (`fastMethod = "dpr1"`), **not** the default route: a direct `eigen()` call (`fastMethod = "direct"`) is
  empirically 6–56x faster at spicyGLM's scale. It's low-level numerical linear algebra — treat it as self-contained
  utility code, don't casually "simplify" without understanding the secular equation math (bracket/bisect
  root-finding for tied vs. distinct poles).
- `R/spicyGLM.R` orchestrates fitting across all cell-type pairs (parallelized via `BiocParallel`), handles missing
  cell types per condition group (`diagnoseMissingConditions()`, `computeCellTypePresence()`), and builds the
  `SpicyResults`-like output plus its `print`/`show` method (`.showSpicyGLMResults()`).

### Data handling

`.format_data()` (`R/utilities.R`) is the common ingestion point: it normalizes `SingleCellExperiment`,
`SpatialExperiment`, or plain data frame input into the internal representation spicyR's spatial statistics code
expects (columns for imageID, cellType, and x/y spatial coordinates). Accessor helpers (`getCellType`, `getImageID`,
`getImagePheno`, etc.) are used throughout `spicy.R`/`spicyGEE.R` rather than indexing `colData`/columns directly —
prefer these accessors over ad hoc extraction when working with `cells` objects.

`convPairs()` (`R/convPairs.R`) converts a `colPairs` spatial-graph object (e.g. from `imcRtools::buildSpatialGraph`)
into a per-image cell-type-pair abundance matrix, usable as `alternateResult` input to `spicy()`.

### Results and visualization

- `SpicyResults` (`R/AllClasses.R`) is a thin S4 class wrapping a `list`, returned by `spicy()`.
- `topPairs()` (`R/AllGenerics.R`) is an S4 generic for extracting/ranking significant results from a
  `SpicyResults` object.
- `signifPlot()`, `spicyBoxPlot()`, `plotImage()`, `imageCrossPlot()` provide `ggplot2`-based visualizations of
  results and raw spatial data.

### Data

`data/diabetesData.rda` and `data/spicyTest.rda` are bundled example datasets (see `R/data.R` for docs) used
throughout examples, vignettes, and tests. `inst/extdata/isletCells.txt.gz` and `inst/testdata/original_result.rds`
are additional fixtures used by tests/vignettes — `original_result.rds` is the golden snapshot the main regression
test compares against.

## Documentation source of truth

Function documentation lives as roxygen2 (`#'`) comments directly above each exported function in `R/*.R`.
`man/*.Rd` and `NAMESPACE` are generated — do not hand-edit them; edit the roxygen comments and run
`devtools::document()`. The vignette (`vignettes/spicyR.Rmd`) is the primary end-to-end usage walkthrough.

## Working conventions

- **Plan before code.** Agree on the qualitative approach before writing anything. Don't jump straight to
  implementation.
- **No new `.qmd` files** unless explicitly requested.
- **Minimal code comments.** Placeholder comments in production code are treated as an error.
- **Step-by-step with empirical validation.** Verify each change against known behavior before building on it.
  Verify mathematical claims numerically before accepting them as correct.
- **Be precise.** Don't overstate confidence or make assumptions about data structures without checking.
  Corrections to imprecise notation/framing should be taken on board immediately, not argued with.
- **Keep responses brief.** Avoid verbose or circular explanations.

## Key mathematical results (already established — don't re-derive from scratch)

- Leverage share `ℓ_i` is exactly pre-fit computable (the `exp(β̂)` term cancels).
- Variance-share influence `e_i²/v_L̂` is an exact decomposition of reported variance, not an approximation.
- Closed-form leave-one-out `Δlog RR_i` requires only subtraction of patient totals.
- Satterthwaite `ν` is a structural leverage-concentration penalty, distinct from influence.
- CR2 never dampens residuals — `A_i` eigenvalues are always ≥ 1. Low-leverage/high-residual patients pass through
  nearly untouched — a known blind spot for leverage-only diagnostics.
- Leverage is pre-outcome/structural/controllable (via design weights); influence is data-specific and not
  pre-controllable. High-leverage/low-influence patients can be legitimate high-quality contributors, not
  necessarily a problem.
- Firth closed form: group-sum-of-leverages collapses to 1 regardless of β in the rank-1 design →
  `β̂_Firth = log((Y_g + 0.5)/D_g)`. Jacobian-invariance means the existing fast CR2 machinery applies unmodified.
- `logRateRatio` is mathematically identical for A→B and B→A (exact) — bidirectional pair dedup is principled,
  not a heuristic.
- Under working independence, GEE and GLM+sandwich give identical point estimates and identical CR2 variance.
  No efficiency gain from richer working correlation structure in single-image-per-patient settings.
- CR2 requires ≥2 distinct clusters per condition group. Longitudinal/repeated-measures designs (e.g. pre/post)
  are not handled by the current framework without a paired/mixed-effects extension.
- Binomial (fixed-k NN design): the CR2/Satterthwaite chain is the Poisson chain under the substitution
  `μ_ij → v_ij = k·π̂(1−π̂)` (working variance) everywhere it feeds `Λ_i, f_i, T_i, S_g, H_X,ii, G_i, A_i,
  P_j^diag`, **and also `τ̃_j = 1ᵀ A_j (∂μ_j/∂β) = 1ᵀ A_j v_j`** in the d.o.f. — `∂μ/∂β` is the working
  variance for the logit link, coinciding with the mean only for Poisson (`V = diag(μ)`). The residual stays
  count-scale: `r̂_i = Y_count − k·π̂` (= `k × proportion residual`). `spicyClub_Supplementary_Math.pdf` §12.3
  mis-states `τ̃_j` as using the true mean — that is a doc error; the code uses the working variance and matches
  `clubSandwich` d.o.f. exactly.

## Known methodological discrepancy (documented for thesis, not a bug to fix)

spicyR silently computes results for pairs where a cell type is absent from one condition (assigning near-zero
floor values). spicyGLM explicitly skips and logs these instead. This is intentional and should be preserved,
not "fixed" to match spicyR's behavior.

## Known outstanding issue

`subject = NULL` threading is currently broken somewhere in `getPairwiseAssoc()` / the CR2 machinery — needs
resolving.

## repomix workflow (for generating a full-repo context dump)

```bash
npx repomix --remote SydneyBioX/spicyR --remote-branch gee --style xml
```
Output: `repomix-output.xml` in the current working directory.

## Key references

- Pustejovsky & Tipton (2018) — clubSandwich
- Bell & McCaffrey (2002) — BRL
- Jakovčević Stor et al. (2015) — DPR1 eigensolver
- Thesis mathematical derivation doc: `spicyClub_maths.pdf`
