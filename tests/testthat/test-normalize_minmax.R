test_that("normalize_min_max() scales image values", {
  path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
  normalized <- normalize_min_max(read_caim(path), 0, 255)

  expect_equal(min(normalized[]), 0)
  expect_equal(max(normalized[]), 1)
})
