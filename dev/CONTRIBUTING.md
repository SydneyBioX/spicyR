# Contributing to spicyR

spicyR’s cell-level test runs in a C++17 core (`src/core`, with Eigen)
that is shared verbatim with the Python package
[spicyr](https://github.com/SydneyBioX/spicyr-py). The files in
`src/*.cpp` are the R bindings only, and the R code in `R/cell.R` and
`R/spicy_main.R` prepares inputs and assembles results without doing any
statistics.

- Change the core here, then copy it to spicyr with its
  `tools/sync_core.sh`. spicyr’s tests check the copy.
- `tests/testthat/test-cell.R` checks the core against a slow plain-R
  implementation written from the definitions
  (`tests/testthat/helper-reference.R`). Keep them in step.
- spicyr’s `tests/shared_cases/make_cases.R` writes cases from this
  package that spicyr must reproduce; rerun it after any change to the
  core or to `R/cell.R` or `R/spicy_main.R`.
- `R CMD check` and `BiocCheck` should stay clean.
