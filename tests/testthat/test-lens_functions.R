test_that("lens return numeric vectors", {
  local_edition(3)
  expect_type(lens("Nikon_FCE9"), "double")
})

