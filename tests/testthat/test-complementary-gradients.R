complementary_test_cfa <- function(values, type = "bayer") {
  r <- set_rcaiman_geometry(terra::rast(nrows = 1, ncols = nrow(values), nlyrs = 4))
  terra::values(r) <- values
  names(r) <- if (type == "bayer") c("Red", "Green1", "Green2", "Blue") else
    c("Cyan", "Yellow", "Green", "Magenta")
  r
}

test_that("Bayer retains the original gradients and averages both greens", {
  r <- complementary_test_cfa(rbind(c(1, 0, 0, 0), c(0, 1, 1, 0),
                                    c(0, 0, 0, 1), c(1, 2, 6, 9)))
  out <- complementary_gradients(r)
  expect_equal(unname(out[]), rbind(c(1, -1, -1), c(-1, -1, 1),
                                    c(1, 1, 1), c(6, 4, 12) / 14))
  expect_identical(names(out), c("green_magenta", "yellow_blue", "red_cyan"))
  expect_true(terra::compareGeom(out, r))
  swapped <- r
  swapped$Green1 <- r$Green2
  swapped$Green2 <- r$Green1
  expect_equal(complementary_gradients(swapped)[], out[])
})

test_that("CYGM uses direct contrasts with the agreed polarity", {
  r <- complementary_test_cfa(rbind(c(8, 2, 2, 8), c(2, 8, 8, 2),
                                    c(8, 2, 2, 2)), "cygm")
  out <- complementary_gradients(r, "cygm")
  expect_equal(unname(out[]), rbind(c(.6, 7/9, -1/9),
                                    c(-.6, -1/3, -2/3), c(0, 2/3, 1/3)))
  expect_identical(names(out), c("green_magenta", "yellow_blue", "red_cyan"))
  expect_true(terra::compareGeom(out, r))
})

test_that("equal positive bands have the documented nonzero responses", {
  for (type in c("bayer", "cygm")) {
    r <- complementary_test_cfa(matrix(10, nrow = 4, ncol = 4), type)
    expected <- if (type == "bayer") c(1, -1, 1)/3 else c(0, 1, -1)/3
    expect_equal(unname(complementary_gradients(r, type)[]),
                 matrix(rep(expected, each = 4), ncol = 3))
  }
})

test_that("band order and common scaling do not change the gradients", {
  values <- rbind(c(1, 2, 4, 3), c(4, 3, 1, 8), c(6, 2, 5, 4))
  for (type in c("bayer", "cygm")) {
    r <- complementary_test_cfa(values, type)
    expected <- complementary_gradients(r, type)[]
    expect_equal(complementary_gradients(r[[c(4, 2, 1, 3)]], type)[], expected)
    expect_equal(complementary_gradients(r * 100, type)[], expected)
    padded <- complementary_test_cfa(rbind(matrix(0, 40, 4), values,
                                           matrix(NA_real_, 4, 4)), type)
    expect_equal(complementary_gradients(padded, type)[][41:43, ], expected)
  }
})

test_that("missing values and zero denominators affect the relevant gradients", {
  r <- complementary_test_cfa(rbind(c(0, 0, 0, 0), c(NA, 2, 3, 4),
    c(1, NA, 3, 4), c(1, 2, NA, 4), c(1, 2, 3, NA), c(1, 2, 3, 4)))
  out <- complementary_gradients(r)[]
  expect_true(all(is.na(out[1:5, ])))
  expect_true(all(is.finite(out[6, ])))
  r <- complementary_test_cfa(rbind(c(NA, 2, 3, 4), c(1, NA, 3, 4),
    c(1, 2, NA, 4), c(1, 2, 3, NA), c(0, 0, 3, 0), c(2, 1, 0, 0)), "cygm")
  out <- unname(complementary_gradients(r, "cygm")[])
  expect_equal(out, rbind(c(1/7, NA, NA), c(1/7, NA, NA),
    c(NA, 3/7, -5/7), c(NA, NA, NA), c(-1, NA, NA), c(NA, 1/3, 1/3)))
})

test_that("CFA band names are explicit contracts", {
  r <- complementary_test_cfa(matrix(1, 1, 4))
  expect_error(complementary_gradients(r, "cygm"), "Cyan, Yellow, Green, Magenta")
  expect_error(complementary_gradients(r, "rgb"), "arg")
  expect_error(complementary_gradients(r[[1:3]]), "four uniquely named")
  names(r) <- c("Red", "Green", "Green", "Blue")
  expect_error(complementary_gradients(r), "Green1, Green2")
})

test_that("split_cfa supplies the bands directly for either family", {
  for (type in c("bayer", "cygm")) {
    if (type == "bayer") {
      pattern <- matrix(c("Red", "Green1", "Green2", "Blue"), 2, byrow = TRUE)
      block <- c(1, 2, 6, 9)
      expected <- c(6, 4, 12)/14
    } else {
      pattern <- matrix(c("Cyan", "Yellow", "Green", "Magenta"), 2, byrow = TRUE)
      block <- c(8, 2, 2, 8)
      expected <- c(.6, 7/9, -1/9)
    }
    raw <- set_rcaiman_geometry(terra::rast(nrows = 4, ncols = 4))
    terra::values(raw) <- rep(c(rep(block[1:2], 2), rep(block[3:4], 2)), 2)
    caim <- split_cfa(raw, pattern)
    out <- complementary_gradients(caim, type)
    expect_equal(unname(out[]), matrix(rep(expected, each = 4), ncol = 3))
    expect_true(terra::compareGeom(out, caim))
  }
})
