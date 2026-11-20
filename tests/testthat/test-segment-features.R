feature_test_layers <- function() {
  r <- set_rcaiman_geometry(terra::rast(nrows=2,ncols=3))
  terra::values(r) <- c(2,4,NA,8,NA,99)
  seg <- terra::rast(r)
  terra::values(seg) <- c(9000000000,9000000000,20,0,20,NA)
  d <- c(r,r)
  terra::values(d) <- cbind(c(120,130,40,50,60,70),c(10,20,30,40,50,60))
  names(d) <- c("zenith","azimuth")
  w <- terra::rast(r)
  terra::values(w) <- c(1,3,1,1,1,NA)
  list(r=r,seg=seg,d=d,w=w)
}

test_that("features use external directions and preserve sparse segment IDs", {
  x <- feature_test_layers()
  f <- compute_features(x$r,x$seg,x$d,x$w,ignore_label_0=FALSE)
  expect_equal(f$id,c(0,20,9000000000))
  expect_equal(f$mean,c(8,NA,3.5))
  expect_equal(f$pixel_count,c(1,2,2))
  expect_equal(f$valid_count,c(1,0,2))
  expect_equal(f$zenith[3],130*pi/180)
  expect_equal(f$azimuth[3],20*pi/180)
  expect_false(any(c("zen","azi") %in% names(f)))
  image <- feature_image(x$seg,f,c("mean","area"),ignore_label_0=FALSE)
  expect_equal(terra::values(image)[,1],c(3.5,3.5,NA,8,NA,NA))
  expect_equal(names(image),c("mean","area"))
  expect_true(terra::compareGeom(image[[1]],x$r))
})

test_that("means include empty measurements and consistently handle label zero", {
  x <- feature_test_layers()
  expect_equal(compute_mean(x$r,x$seg,return="vector"),
               c(`20`=NA_real_,`9000000000`=3))
  expect_equal(compute_mean(x$r,x$seg,return="vector",ignore_label_0=FALSE),
               c(`0`=8,`20`=NA_real_,`9000000000`=3))
  expect_equal(terra::values(compute_mean(x$r,x$seg,ignore_label_0=FALSE),mat=FALSE),
               c(3,3,NA,8,NA,NA))
  x$seg[] <- NA
  expect_length(compute_mean(x$r,x$seg,return="vector"),0)
  expect_true(all(is.na(terra::values(compute_mean(x$r,x$seg)))))
  expect_equal(nrow(compute_features(x$r,x$seg,x$d,x$w)),0)
})

test_that("feature tables are matched by unique ID rather than row position", {
  x <- feature_test_layers()
  f <- data.frame(id=c(9000000000,0,100),value=c(7,8,9),selected=c(TRUE,FALSE,TRUE))
  out <- feature_image(x$seg,f,c("value","selected"),ignore_label_0=FALSE)
  expect_equal(terra::values(out)[,1],c(7,7,NA,8,NA,NA))
  expect_equal(terra::values(out)[,2],c(1,1,NA,0,NA,NA))
  expect_true(is.na(terra::values(feature_image(x$seg,f,"value"))[4]))
  expect_error(feature_image(x$seg,rbind(f,f[1,]),"value"),"unique")
  x$seg[1] <- 1.5
  expect_error(compute_mean(x$r,x$seg),"integer")
})

test_that("missing weights affect weighted statistics but not observed extrema", {
  x <- feature_test_layers()
  x$w[] <- NA
  f <- compute_features(x$r,x$seg,x$d,x$w)
  expect_true(all(is.na(f$mean)))
  expect_equal(f$min,c(NA,2))
  expect_equal(f$max,c(NA,4))
  expect_true(all(is.na(f$shape_complexity)))
})
