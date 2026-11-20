test_that("a sky mask can be used directly as a segmentation", {
  geometry <- list(model="radial_projection",parameters=equidistant_coef,
    center_of_vision=c(10,8),normalization_radius=8,normalization_angle=90)
  g <- set_internal_geometry(c(20,16),geometry)
  d <- direction_image(g)
  m <- !is.na(d$zenith) & d$zenith >= 30 & d$zenith <= 60
  seg <- as.numeric(m)
  r <- d$zenith
  r[is.na(r)] <- 0
  w <- area_weight_image(g,geometry)
  expect_equal(compute_mean(r,m,return="vector"),compute_mean(r,seg,return="vector"))
  expect_equal(compute_mean(r,m,return="vector",ignore_label_0=FALSE),
               compute_mean(r,seg,return="vector",ignore_label_0=FALSE))
  f <- compute_features(r,m,d,w)
  expect_equal(f,compute_features(r,seg,d,w))
  expect_equal(terra::values(feature_image(m,f,"mean")),
               terra::values(feature_image(seg,f,"mean")))
  samples <- data.frame(x=c(5.5,8.5),y_down=c(5.5,8.5),dn=c(2,4))
  expect_equal(aggregate_by_segment(samples,m,return="vector"),
               aggregate_by_segment(samples,seg,return="vector"))
  expect_equal(terra::values(merge_segments(r,w,m,scale_parameter=0)),
               terra::values(merge_segments(r,w,seg,scale_parameter=0)))
  expect_equal(terra::values(binarize_by_region(r,m,"thr_isodata")),
               terra::values(binarize_by_region(r,seg,"thr_isodata")))
  expect_equal(terra::values(defuzzify(r/180,m)),terra::values(defuzzify(r/180,seg)))
  expect_error(.check_vector(TRUE,"integerish"),"integer")
  for (bad in c(-1,.5,Inf)) {
    invalid <- seg
    invalid[1] <- bad
    expect_error(compute_mean(r,invalid))
  }
})
