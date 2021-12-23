test_that("normalize() works", {
  expect_equal(min(normalize(read_caim(), 0, 255)[]),
               0)
  expect_equal(max(normalize(read_caim(), 0, 255)[]),
               1)
})
