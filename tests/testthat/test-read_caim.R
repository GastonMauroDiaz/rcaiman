test_that("read_caim() assigns RGB layer names", {
  path <- system.file("external/DSCN4500.JPG", package = "rcaiman")

  expect_setequal(
    names(read_caim(path)),
    c("Red", "Green", "Blue")
  )
})
