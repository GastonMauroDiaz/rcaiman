test_that("angle unit conversion works", {
  expect_equal(radian2degree(pi/2), 90)
  expect_equal(degree2radian(90), pi/2)
})
