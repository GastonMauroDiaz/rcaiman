landmark_g <- function(polar = c(12, 24, 38, 52, 67, 81),
                        azimuth = c(15, 85, 160, 220, 275, 330)) {
  g <- terra::rast(nrows = 1, ncols = length(polar), nlyrs = 2)
  names(g) <- c("off_axis", "internal_azimuth")
  terra::values(g) <- cbind(polar, azimuth)
  g
}

landmark_tables <- function(g, orientation) {
  indices <- data.frame(col = seq_len(terra::ncell(g)), row = 1)
  d <- do.call(direction_image, c(list(g), as.list(orientation)))
  list(internal = add_directions(indices, g = g),
       reference = add_directions(indices, d = d))
}

test_that("landmarks recover arbitrary orientations through the public workflow", {
  g <- landmark_g()
  orientations <- list(c(242, 240, 83), c(90, 90, 32.5), c(0, 0, 0),
                       c(17, 231, -24), c(270, 32, 89.99999),
                       c(37, 210, -89.99999))
  for (orientation in orientations) {
    tables <- landmark_tables(g, orientation)
    before <- tables
    fit <- fit_image_orientation(tables$internal, tables$reference)
    expect_equal(unname(unlist(fit$orientation)), orientation, tolerance = 1e-6)
    expect_lt(max(fit$residuals$angular_error), 1e-12)
    d <- do.call(direction_image, c(list(g), fit$orientation))
    expect_equal(terra::values(d) * pi / 180,
                 as.matrix(tables$reference[c("zenith", "azimuth")]),
                 tolerance = 1e-12)
    expect_identical(tables, before)
  }
})

test_that("vertical optical axes have a deterministic equivalent orientation", {
  g <- landmark_g()
  for (elevation in c(-90, 90)) {
    tables <- landmark_tables(g, c(37, 123, elevation))
    fit <- fit_image_orientation(tables$internal, tables$reference)
    expect_equal(fit$orientation$sensor_rotation, 0)
    expect_equal(fit$orientation$elevation, elevation)
    expected_yaw <- if (elevation > 0) 123 - 37 else 123 + 37
    expect_equal(fit$orientation$azimuth_rotation, expected_yaw)
    d <- do.call(direction_image, c(list(g), fit$orientation))
    expect_equal(terra::values(d) * pi / 180,
                 as.matrix(tables$reference[c("zenith", "azimuth")]),
                 tolerance = 1e-12)
  }
})

test_that("two nonparallel landmarks suffice and poles need no azimuth", {
  internal <- data.frame(off_axis = c(0, pi / 2), internal_azimuth = c(NA, 0))
  reference <- data.frame(zenith = c(pi / 2, 0), azimuth = c(pi, NA))
  fit <- fit_image_orientation(internal, reference)
  expect_equal(unname(unlist(fit$orientation)), c(0, 0, 0), tolerance = 1e-12)
  expect_equal(fit$residuals$zenith, reference$zenith)
  expect_true(is.na(fit$residuals$azimuth[2]))
  expect_lt(fit$summary[["max"]], 1e-12)
})

test_that("a reflection is not accepted as a rotation", {
  # Reflect east-west, with unequal numbers of rays along three orthogonal
  # axes. A proper rotation must leave an error; the unique best fit reverses
  # the least represented axis instead of reflecting all rays exactly.
  internal <- data.frame(off_axis = c(rep(pi / 2, 5), 0),
                         internal_azimuth = c(rep(pi / 2, 3), 0, 0, 0))
  reference <- data.frame(zenith = internal$off_axis,
                          azimuth = internal$internal_azimuth)
  fit <- fit_image_orientation(internal, reference)
  expect_equal(fit$residuals$angular_error, c(rep(0, 5), pi), tolerance = 1e-12)
  expect_equal(fit$summary[["rmse"]], pi / sqrt(6), tolerance = 1e-12)
})

test_that("noisy correspondences retain all pairs and expose angular errors", {
  g <- landmark_g()
  tables <- landmark_tables(g, c(24, 59, 81))
  tables$reference$zenith <- tables$reference$zenith +
    c(0.004, -0.005, 0.009, -0.003, 0.006, -0.007)
  tables$reference$id <- rev(seq_len(nrow(tables$reference)))
  fit <- fit_image_orientation(tables$internal, tables$reference)
  expect_equal(nrow(fit$residuals), nrow(tables$reference))
  error <- calc_spherical_distance(fit$residuals$zenith, fit$residuals$azimuth,
                                   tables$reference$zenith, tables$reference$azimuth)
  expect_equal(fit$residuals$angular_error, error, tolerance = 1e-10)
  expect_equal(fit$summary, c(mean = mean(error), rmse = sqrt(mean(error^2)),
                             max = max(error)), tolerance = 1e-10)
  # An independent optimization of the chord-distance objective should not
  # improve on the SVD solution, despite a perturbed starting orientation.
  objective <- function(angles) {
    d <- do.call(external_directions_from_internal,
                 c(list(tables$internal), as.list(angles)))
    v <- as.matrix(d[c("zenith", "azimuth")])
    angles_to_vectors <- function(z, a) cbind(sin(z) * sin(a), sin(z) * cos(a), cos(z))
    sum((angles_to_vectors(v[, 1], v[, 2]) -
           angles_to_vectors(tables$reference$zenith, tables$reference$azimuth))^2)
  }
  initial <- unlist(fit$orientation)
  optimized <- stats::optim(initial + c(1, -1, 0.5), objective,
                            control = list(reltol = 1e-12))
  expect_equal(objective(initial), optimized$value, tolerance = 1e-8)
})

test_that("incomplete and degenerate landmark sets fail explicitly", {
  tables <- landmark_tables(landmark_g(), c(0, 0, 90))
  internal <- tables$internal
  reference <- tables$reference
  expect_error(fit_image_orientation(internal[-1, ], reference), "same number")
  expect_error(fit_image_orientation(internal[1, ], reference[1, ]), "At least two")
  expect_error(fit_image_orientation(internal[c(1, 1), ], reference[c(1, 1), ]),
                "unique rotation")
  antipodal <- data.frame(off_axis = c(0, pi), internal_azimuth = c(NA, NA))
  expect_error(fit_image_orientation(antipodal,
                                     data.frame(zenith = c(0, pi), azimuth = NA)),
                "unique rotation")
  bad <- internal
  bad$off_axis[1] <- NA
  expect_error(fit_image_orientation(bad, reference), "NA")
  bad <- internal
  bad$internal_azimuth[1] <- NA
  expect_error(fit_image_orientation(bad, reference), "Azimuth may be NA")
  bad <- reference
  bad$zenith[1] <- Inf
  expect_error(fit_image_orientation(internal, bad), "Inf")
  bad <- reference
  bad$azimuth[1] <- 2 * pi
  expect_error(fit_image_orientation(internal, bad), "radians")
  expect_error(fit_image_orientation(internal * 180 / pi, reference), "radians")
  expect_error(fit_image_orientation(reference, internal), "off_axis")
})
