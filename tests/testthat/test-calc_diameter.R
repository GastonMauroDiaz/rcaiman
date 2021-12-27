test_that("calc_diameter works", {
  expect_equal(calc_diameter(lens("Nikkor_10.5_mm"), 1202, 53), 3754)
})
