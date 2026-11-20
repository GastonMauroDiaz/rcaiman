openness_test_raster <- function(values) {
  r <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2)
  terra::values(r) <- values
  r
}


openness_test_d <- function(zenith = rep(45, 6)) {
  d <- c(openness_test_raster(zenith), openness_test_raster(rep(NA_real_, 6)))
  names(d) <- c("zenith", "azimuth")
  d
}

test_that("openness integrates weighted gap area and normalizes every selection", {
  bin <- openness_test_raster(c(TRUE,FALSE,TRUE,FALSE,TRUE,FALSE))
  d <- openness_test_d()
  w <- openness_test_raster(c(1,2,3,4,5,6))
  expect_equal(compute_canopy_openness(bin,d,w), 9/21)
  all <- bin | !bin
  expect_equal(compute_canopy_openness(bin,d,w,all), compute_canopy_openness(bin,d,w))
  m <- openness_test_raster(c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE))
  expect_equal(compute_canopy_openness(bin,d,w,m), 1/3)
  expect_equal(compute_canopy_openness(bin,d,w*1e100,m), 1/3)
  expect_equal(compute_canopy_openness(bin,d,w*0+1), .5)
  expect_equal(compute_canopy_openness(all,d,w,m), 1)
  expect_equal(compute_canopy_openness(!all,d,w,m), 0)
  expect_equal(compute_canopy_openness(bin,d,w*0+1e308), .5)
})

test_that("missing observations and zero weights do not act as vegetation", {
  bin <- openness_test_raster(c(TRUE,FALSE,NA,TRUE,FALSE,TRUE))
  d <- openness_test_d()
  w <- openness_test_raster(c(1,2,100,NA,0,3))
  expect_equal(compute_canopy_openness(bin,d,w), 4/6)
  m <- openness_test_raster(c(FALSE,FALSE,TRUE,TRUE,TRUE,FALSE))
  expect_error(compute_canopy_openness(bin,d,w,m), "positive weights")
  expect_error(compute_canopy_openness(bin,d,w*0), "positive weights")
  bin[] <- NA
  expect_error(compute_canopy_openness(bin,d,w), "positive weights")
})

test_that("openness validates grids, binary values and weights", {
  bin <- openness_test_raster(rep(c(TRUE,FALSE),3))
  d <- openness_test_d()
  w <- openness_test_raster(rep(1,6))
  expect_error(compute_canopy_openness(w,d,w), "logical")
  bad <- terra::deepcopy(w)
  bad[1] <- -1
  expect_error(compute_canopy_openness(bin,d,bad), "non-negative")
  bad[1] <- Inf
  expect_error(compute_canopy_openness(bin,d,bad), "Inf")
  m <- terra::deepcopy(bin)
  m[1] <- NA
  expect_error(compute_canopy_openness(bin,d,w,m), "NA-free")
  expect_error(compute_canopy_openness(bin,d,terra::shift(w,dx=1)), "same geometry")
  expect_error(compute_canopy_openness(bin,d,w,terra::shift(bin,dy=1)), "same geometry")
  expect_error(compute_canopy_openness(bin,terra::shift(d,dx=1),w), "same geometry")
  expect_error(compute_canopy_openness(bin,d[[1]],w), "two layers")
  invalid <- d
  names(invalid) <- c("off_axis", "internal_azimuth")
  expect_error(compute_canopy_openness(bin,invalid,w), "zenith and azimuth")
  expect_error(compute_canopy_openness(bin,openness_test_d(c(-1,rep(45,5))),w),
               "degrees")
})

test_that("the external horizon excludes ground and unknown directions from both sums", {
  bin <- openness_test_raster(c(TRUE,FALSE,TRUE,TRUE,FALSE,TRUE))
  w <- openness_test_raster(c(1,2,3,100,100,100))
  d <- openness_test_d(c(0,45,90,90.01,180,NA))
  expect_equal(compute_canopy_openness(bin,d,w), 4/6)
  all <- bin | !bin
  expect_equal(compute_canopy_openness(bin,d,w,all), 4/6)
  known <- !is.na(d$zenith)
  expect_equal(compute_canopy_openness(bin,d,w,known), 4/6)
  ground <- !is.na(d$zenith) & d$zenith > 90
  expect_error(compute_canopy_openness(bin,d,w,ground), "positive weights")
  expect_error(compute_canopy_openness(bin,openness_test_d(rep(NA_real_,6)),w),
               "positive weights")
})

test_that("inclined radial coverage composes with external selection and renormalization", {
  geometry <- list(model="radial_projection",parameters=equidistant_coef,
    center_of_vision=c(40,30),normalization_radius=30,normalization_angle=90)
  g <- set_internal_geometry(c(80,60),geometry)
  d <- direction_image(g,elevation=32.5)
  w <- area_weight_image(g,geometry)
  bin <- g$internal_azimuth < 180
  upper <- !is.na(d$zenith) & d$zenith <= 90
  expect_true(any(d$zenith[] > 90, na.rm = TRUE))
  expected <- sum((bin * w)[upper], na.rm = TRUE) / sum(w[upper], na.rm = TRUE)
  expect_equal(compute_canopy_openness(bin,d,w), expected)
  expect_equal(compute_canopy_openness(bin,d,w,upper), expected)
  m <- !is.na(d$zenith) & d$zenith <= 70
  result <- compute_canopy_openness(bin,d,w,m)
  expect_true(is.finite(result) && result > 0 && result < 1)
  expect_equal(compute_canopy_openness(!is.na(g$off_axis),d,w,m),1)
  geometry$normalization_radius <- 20
  geometry$normalization_angle <- 60
  geometry$parameters <- normalize_projection_coef(equidistant_coef,60)
  expect_equal(compute_canopy_openness(bin,d,area_weight_image(g,geometry),m),result)
  selected_weights <- terra::values(w,mat=FALSE)
  selected_weights[!terra::values(m,mat=FALSE)] <- NA
  masked_weights <- terra::rast(w)
  terra::values(masked_weights) <- selected_weights
  expect_equal(compute_canopy_openness(bin,d,masked_weights),result)
})
