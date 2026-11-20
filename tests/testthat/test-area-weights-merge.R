weight_test_geometry <- function() list(model = "radial_projection",
  parameters = equidistant_coef, center_of_vision = c(40,30),
  normalization_radius = 30, normalization_angle = 90)

test_that("radial weights use off-axis angle and ignore acquisition orientation", {
  geometry <- weight_test_geometry()
  g <- set_internal_geometry(c(80, 60), geometry)
  w <- area_weight_image(g, geometry)
  theta <- terra::values(g$off_axis, mat = FALSE) * pi/180
  expected <- sin(theta) / theta * (pi/2)^2
  expect_equal(terra::values(w, mat = FALSE), expected)
  d <- direction_image(g, elevation = 32.5)
  expect_true(any(terra::values(d$zenith) > 90, na.rm = TRUE))
  expect_true(all(is.finite(terra::values(w)[!is.na(theta)])))
  expect_null(attr(w, "projection_coef"))
  copy <- terra::deepcopy(g)
  copy$internal_azimuth <- NA_real_
  expect_equal(terra::values(area_weight_image(copy, geometry)), terra::values(w))
  copy$off_axis <- 0
  expect_equal(unique(as.numeric(terra::values(area_weight_image(copy, geometry)))), (pi/2)^2)
  extent <- terra::ext(20, 60, 30, 50)
  expect_equal(terra::values(area_weight_image(terra::crop(g, extent), geometry)),
               terra::values(terra::crop(w, extent)))
})

test_that("normalization changes relative weights but preserves physical solid angles", {
  geometry <- weight_test_geometry()
  g <- set_internal_geometry(c(80, 60), geometry)
  first <- area_weight_image(g, geometry) / geometry$normalization_radius^2
  updated <- geometry
  updated$parameters <- normalize_projection_coef(geometry$parameters, 60)
  updated$normalization_angle <- 60
  updated$normalization_radius <- geometry$normalization_radius *
    calc_relative_radius(60, geometry$parameters)
  second <- area_weight_image(g, updated) / updated$normalization_radius^2
  expect_equal(terra::values(first), terra::values(second))
  geometry$center_of_vision <- c(100,100)
  geometry$normalization_radius <- 90
  g <- set_internal_geometry(c(200,200), geometry)
  solid_angle <- area_weight_image(g, geometry) / geometry$normalization_radius^2
  expect_lt(abs(sum(terra::values(solid_angle), na.rm=TRUE) - 2*pi), .01)
  anisotropic <- list(model="anisotropic_projection", center_of_vision=c(40,30),
    parameters=list(off_axis_coef=1, sin_internal_azimuth_coef=c(0,1),
      cos_internal_azimuth_coef=c(-1,0), atan_scale=.6))
  expect_error(area_weight_image(g, anisotropic), "radial")
})

test_that("merging handles partial inclined grids and sparse IDs without attributes", {
  geometry <- weight_test_geometry()
  g <- set_internal_geometry(c(80,60), geometry)
  d <- direction_image(g, elevation=32.5)
  seg <- healpix_segmentation(d, 4)
  w <- area_weight_image(g, geometry)
  r <- g$off_axis
  before <- terra::values(seg)
  merged <- merge_segments(r, w, seg, scale_parameter=2)
  expect_true(terra::compareGeom(merged, seg))
  expect_equal(terra::values(seg), before)
  expect_true(all(terra::values(merged)[before == 0] == 0))
  expect_true(all(terra::values(merged) %in% before))
  expect_lte(length(unique(terra::values(merged))), length(unique(before)))

  seg <- terra::rast(nrows=2,ncols=3)
  terra::values(seg) <- c(0, 800000001, 800000002, NA, 800000001, 800000002)
  r <- seg * 0 + 1
  w <- seg * 0 + 1
  r[1] <- NA
  w[1] <- NA
  frozen <- merge_segments(r,w,seg,1,max_iterations=0)
  expect_equal(as.numeric(terra::values(frozen)), c(0,800000001,800000002,0,800000001,800000002))
  merged <- merge_segments(r,w,seg,2)
  expect_equal(length(unique(terra::values(merged)[terra::values(merged)>0])),1L)
  expect_true(all(terra::values(merged) %in% c(0,800000001,800000002)))
  r[2] <- NA
  expect_error(merge_segments(r,w,seg,1), "included segments")
})

test_that("merging empty or isolated segments returns their labels safely", {
  seg <- terra::rast(nrows=2,ncols=3)
  terra::values(seg) <- c(1,0,2,1,0,2)
  r <- w <- seg * 0 + 1
  out <- merge_segments(r,w,seg,2)
  expect_equal(terra::values(out), terra::values(seg))
  seg[] <- 0
  expect_equal(as.numeric(terra::values(merge_segments(r,w,seg,2))), rep(0,6))
  seg[] <- c(0,1,1,0,1,1)
  out <- merge_segments(r,w,seg,2,ignore_label_0=FALSE)
  expect_equal(length(unique(terra::values(out))),1L)
  expect_error(merge_segments(r,w,seg,2,max_iterations=.5), "integer")
})
