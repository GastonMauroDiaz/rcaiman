# Working with the tests

The suite protects scientific calculations and the workflows that connect them.
The number printed by testthat counts individual expectations, including those
inside loops; it is not the number of independent scenarios.

## During development

Run the affected files while changing a function. `filter` matches test filenames
and accepts regular expressions, so a small workflow can be tested together:

```r
testthat::test_local(filter = "is-outlier|rem-outliers|validate-cie-model")
testthat::test_local(filter = "direction|reproject|crop")
testthat::test_local(filter = "thr-|binarize-by-region")
testthat::test_local(filter = "cie|sun")
```

Run `testthat::test_local()` after a change spanning several workflows, and use
`R CMD check` on the built package before release. Documentation-only edits
usually need their examples checked, rather than the entire numerical suite.

Parallel tests also load rcaiman in separate R processes. When testing source
code, those processes must see the same package version: a stale installed
package can produce differences unrelated to the change. Testing a freshly
built package in an isolated library avoids that mismatch. Tests marked with
`skip_on_cran()` run when `NOT_CRAN` is `true`.

## Choosing what to test

- Keep known scientific results, independent reference calculations, and
  regressions for bugs observed with real data. Keep end-to-end examples that
  exercise units, orientation, weighting, missing coverage, and table/raster
  correspondence.
- Use the smallest fixture that still expresses the problem. Use real files
  for file-format behavior; use small synthetic rasters for raster arithmetic.
- Test a shared numerical rule or validator where it lives. Its wrappers need
  a representative integration check, not every combination already tested
  below them.
- Retain meaningful boundaries, such as the horizon, poles, image edges,
  singular fits, and zero area. Prefer representative combinations over a
  Cartesian product of otherwise independent options.
- Check returned results and useful public contracts. Avoid freezing private
  intermediate objects, repeating a function call just to compare it with
  itself, or permanently testing that retired functions remain absent.
- A deliberate redesign can change an expected result. Review the scientific
  reason for the test before updating or removing it; a failing test does not
  by itself make the old design preferable.

There is no target number of tests. Add a test when it protects a distinct
behavior or catches a plausible regression, not simply because code changed.
