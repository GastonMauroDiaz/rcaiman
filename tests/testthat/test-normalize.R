test_that("normalize() works", {
  path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
  expect_equal(min(normalize(read_caim(path), 0, 255)[]),
               0)
  expect_equal(max(normalize(read_caim(path), 0, 255)[]),
               1)
})
