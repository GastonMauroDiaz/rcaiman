test_that("a call to read_caim() returs a RasterBrick", {
  expect_is(r <- read_caim(), "RasterBrick")
  expect_equal(nlayers(r), 3)
})

test_that("read_caim() assigns good layer names", {
  expect_setequal(names(read_caim()), c("Red", "Green", "Blue"))
})
