centroid_features <- function(z, a, w, values = rep(1, length(z))) {
  layer <- function(x) {
    r <- terra::rast(nrows = 1, ncols = length(z))
    terra::values(r) <- x
    r
  }
  d <- c(layer(z), layer(a))
  names(d) <- c("zenith", "azimuth")
  compute_features(layer(values), layer(rep(1, length(z))), d, layer(w))
}
centroid_angles <- function(x) {
  unlist(x[c("spherical_centroid_zenith", "spherical_centroid_azimuth")],
         use.names = FALSE)
}

test_that("spherical centroid follows represented area rather than pixel density", {
  balanced <- centroid_features(c(60, 60), c(0, 90), c(1, 1))
  # Subdivide only the eastward contribution without changing its total area.
  subdivided <- centroid_features(rep(60, 5), c(0, rep(90, 4)), c(1, rep(.25, 4)))
  expect_equal(centroid_angles(subdivided), centroid_angles(balanced))
  unweighted <- centroid_features(rep(60, 5), c(0, rep(90, 4)), rep(1, 5))
  expect_gt(unweighted$spherical_centroid_azimuth, balanced$spherical_centroid_azimuth)
  scaled <- centroid_features(rep(60, 5), c(0, rep(90, 4)), 13 * c(1, rep(.25, 4)))
  expect_equal(centroid_angles(scaled), centroid_angles(balanced))
  expect_equal(subdivided$planar_centroid_col, 3)
  expect_equal(subdivided$planar_centroid_row, 1)
})

test_that("area weighted directions cross the azimuth seam and ignore radiometry", {
  f <- centroid_features(c(60, 60), c(350, 10), c(1, 3), c(NA, 100))
  # Both directions have the same vertical component; their weighted horizontal
  # components yield this independent azimuth calculation.
  expected <- atan2(2 * sin(10*pi/180), 4 * cos(10*pi/180))
  expect_equal(f$spherical_centroid_azimuth, expected)
  expect_equal(centroid_angles(f),
               centroid_angles(centroid_features(c(60,60), c(350,10), c(1,3))))
})

test_that("zero area and cancelling directions have no spherical centroid", {
  f <- centroid_features(c(45, NA, NA), c(25, NA, NA), c(1, 0, NA))
  expect_equal(centroid_angles(f), c(45, 25) * pi/180)
  expect_true(all(is.na(centroid_angles(centroid_features(c(45,60), c(0,90), c(0,NA))))))
  expect_true(all(is.na(centroid_angles(centroid_features(c(90,90), c(0,180), c(1,1))))))
  expect_true(all(is.na(centroid_angles(centroid_features(c(45,NA), c(0,NA), c(1,1))))))
})
