fov_fixture <- function(mask) {
  r <- terra::rast(nrows = nrow(mask), ncols = ncol(mask), nlyrs = 2,
                   xmin = 0, xmax = ncol(mask), ymin = 0, ymax = nrow(mask), crs = "")
  id <- seq_len(terra::ncell(r))
  angles <- cbind(off_axis = id / 100, internal_azimuth = 100 + id / 100)
  angles[!as.vector(t(mask)), ] <- NA_real_
  terra::values(r) <- angles
  names(r) <- colnames(angles)
  r
}

test_that("a rectangular field follows border centers clockwise and closes", {
  g <- fov_fixture(matrix(TRUE, 3, 4))
  before <- terra::values(g)
  ids <- c(1,2,3,4,8,12,11,10,9,5,1)
  fov <- extract_fov(g)
  expect_identical(fov, as.data.frame(before[ids, ] * pi / 180))
  expect_identical(terra::values(g), before)
  expect_identical(class(fov), "data.frame")
  expect_identical(names(fov), c("off_axis", "internal_azimuth"))
  expect_identical(extract_fov(terra::shift(g, 17, -25)), fov)
})

test_that("a concave window retains its indentation without a hull", {
  mask <- matrix(TRUE, 4, 4)
  mask[3:4, 3:4] <- FALSE
  g <- fov_fixture(mask)
  ids <- c(1,2,3,4,8,7,10,14,13,9,5,1)
  expect_identical(extract_fov(g),
                   as.data.frame(terra::values(g)[ids, ] * pi / 180))
})

test_that("separate windows, corner contacts, and holes are explicit errors", {
  mask <- matrix(TRUE, 5, 5)
  mask[3,3] <- FALSE
  expect_error(extract_fov(fov_fixture(mask)), "more than one contour")
  mask <- matrix(TRUE, 5, 5)
  mask[,3] <- FALSE
  expect_error(extract_fov(fov_fixture(mask)), "more than one contour")
  expect_error(extract_fov(fov_fixture(diag(TRUE, 2))),
               "more than one contour")
})

test_that("thin fields preserve the walk and empty fields keep the schema", {
  for (mask in list(matrix(TRUE, 1, 4), matrix(TRUE, 4, 1))) {
    g <- fov_fixture(mask)
    ids <- c(1,2,3,4,3,2,1)
    expect_identical(extract_fov(g),
                     as.data.frame(terra::values(g)[ids, ] * pi / 180))
  }
  g <- fov_fixture(matrix(TRUE,1,1))
  expect_identical(extract_fov(g), as.data.frame(terra::values(g)*pi/180))
  g <- fov_fixture(matrix(FALSE,3,4))
  expect_identical(extract_fov(g),
                   data.frame(off_axis=numeric(), internal_azimuth=numeric()))
})

test_that("either missing angle excludes a cell and angles do not order the walk", {
  g <- fov_fixture(matrix(TRUE,3,4))
  values <- terra::values(g)
  values[1,1] <- NA_real_
  values[12,2] <- NA_real_
  values[2,2] <- 359.9
  values[3,2] <- .1
  terra::values(g) <- values
  ids <- c(2,3,4,8,11,10,9,5,2)
  expect_identical(extract_fov(g), as.data.frame(values[ids, ]*pi/180))
  expect_error(extract_fov(g[[1]]), "layers")
  values[2,1] <- Inf
  terra::values(g) <- values
  expect_error(extract_fov(g), "Inf")
})

test_that("radial and anisotropic geometries yield original border directions", {
  radial <- set_internal_geometry(c(40,30), list(model="radial_projection",
    parameters=equidistant_coef, center_of_vision=c(20,15),
    normalization_radius=18, normalization_angle=90))
  anisotropic <- set_internal_geometry(c(20,26), list(model="anisotropic_projection",
    center_of_vision=c(10,13), parameters=list(off_axis_coef=c(1.1,.2),
      sin_internal_azimuth_coef=c(0,1), cos_internal_azimuth_coef=c(-1,0), atan_scale=.7)))
  for (g in list(radial, anisotropic)) {
    fov <- extract_fov(g)
    expect_gt(nrow(fov), 3)
    expect_equal(unname(as.numeric(fov[1, ])), unname(as.numeric(fov[nrow(fov), ])))
    # Independently identify the border from four-neighbor erosion.
    mask <- !is.na(g$off_axis) & !is.na(g$internal_azimuth)
    neighbors <- terra::focal(mask, w=matrix(c(NA,1,NA,1,1,1,NA,1,NA),3),
                              fun="sum", fillvalue=0)
    border <- terra::values(mask,mat=FALSE) & terra::values(neighbors,mat=FALSE) < 5
    expected <- as.data.frame(terra::values(g)[border, , drop=FALSE]*pi/180)
    key <- function(x) paste(x$off_axis, x$internal_azimuth)
    expect_setequal(key(fov),key(expected))
  }
})
