aggregate_test_seg <- function() {
  seg <- set_rcaiman_geometry(terra::rast(nrows=2,ncols=3))
  terra::values(seg) <- c(1,1,2,1,2,3)
  seg
}

test_that("aggregation summarizes observations rather than occupied cell means", {
  seg <- aggregate_test_seg()
  samples <- data.frame(x=c(.2,1.4,2.6,1.8,.2),
    y_down=c(.3,.7,.4,1.2,.3),dn=c(2,6,10,20,100))
  result <- aggregate_by_segment(samples,seg,return="vector")
  expect_equal(result,c(`1`=36,`2`=15,`3`=NA_real_))
  raster <- aggregate_by_segment(samples,seg)
  expect_equal(as.numeric(terra::values(raster)),c(36,36,15,36,15,NA))
  expect_true(terra::compareGeom(raster,seg))
  expect_identical(names(raster),"dn")
  expect_equal(aggregate_by_segment(samples,seg,fun=length,return="vector"),
    c(`1`=3,`2`=2,`3`=NA_real_))
  expect_equal(aggregate_by_segment(samples,seg,col_id=3,return="vector"),result)
  expect_equal(aggregate_by_segment(samples,seg,fun=quantile,probs=.5,return="vector"),
    c(`1`=6,`2`=15,`3`=NA_real_))
})

test_that("missing values belong to the chosen summary function", {
  seg <- aggregate_test_seg()
  samples <- data.frame(x=c(.5,1.5,2.5),y_down=.5,dn=c(2,NA,NA))
  expect_true(all(is.na(aggregate_by_segment(samples,seg,return="vector"))))
  result <- aggregate_by_segment(samples,seg,na.rm=TRUE,return="vector")
  expect_equal(result[1],c(`1`=2))
  expect_true(is.nan(result[2]))
  expect_true(is.na(result[3]))
  expect_equal(aggregate_by_segment(samples,seg,fun=length,return="vector"),
    c(`1`=2,`2`=1,`3`=NA_real_))
  expect_error(aggregate_by_segment(samples,seg,fun=range),"one numeric or logical")
  samples$dn <- c(TRUE,FALSE,TRUE)
  expect_equal(aggregate_by_segment(samples,seg,return="vector"),
    c(`1`=.5,`2`=1,`3`=NA_real_))
})

test_that("containing cells, background and invalid coordinates determine membership", {
  seg <- aggregate_test_seg()
  terra::values(seg) <- c(10,20,0,10,NA,9000000000)
  # x=1 selects the second cell; the right and bottom outer edges are outside.
  samples <- data.frame(x=c(0,1,3,0,-.1,NA,Inf,2.5,1.5,2.5),
    y_down=c(0,0,0,2,0,0,0,.5,1.5,1.5),dn=seq_len(10))
  expected <- c(`10`=1,`20`=2,`9000000000`=10)
  expect_equal(aggregate_by_segment(samples,seg,return="vector"),expected)
  expect_equal(aggregate_by_segment(samples,seg,return="vector",ignore_label_0=FALSE),
    c(`0`=8,expected))
  terra::ext(seg) <- terra::ext(100,106,200,204)
  expect_equal(aggregate_by_segment(samples,seg,return="vector"),expected)
  out <- aggregate_by_segment(samples,seg)
  expect_true(all(is.na(terra::values(out)[c(3,5)])))
})

test_that("empty segments are represented without invoking the summary", {
  seg <- aggregate_test_seg()
  samples <- data.frame(x=numeric(),y_down=numeric(),dn=numeric())
  fun <- function(x) stop("must not be called")
  expect_equal(aggregate_by_segment(samples,seg,fun=fun,return="vector"),
    c(`1`=NA_real_,`2`=NA_real_,`3`=NA_real_))
  expect_true(all(is.na(terra::values(aggregate_by_segment(samples,seg,fun=fun)))))
  seg[] <- NA
  expect_length(aggregate_by_segment(samples,seg,fun=fun,return="vector"),0)
  expect_true(all(is.na(terra::values(aggregate_by_segment(samples,seg,fun=fun)))))
})

test_that("summaries see original row order and malformed inputs are rejected", {
  seg <- aggregate_test_seg()
  samples <- data.frame(x=c(.2,1.5,.7),y_down=.5,dn=c(10,20,30))
  expect_equal(aggregate_by_segment(samples,seg,fun=function(x) x[1],return="vector"),
    c(`1`=10,`2`=NA_real_,`3`=NA_real_))
  expect_error(aggregate_by_segment(samples,seg,col_id=2.5),"integer")
  expect_error(aggregate_by_segment(samples,seg,fun=function(x) "text"),"numeric or logical")
  expect_error(aggregate_by_segment(samples,seg,col_id="absent"),"required columns")
  seg[1] <- 1.5
  expect_error(aggregate_by_segment(samples,seg),"integer")
})

test_that("pixel-center observations reproduce raster segment means", {
  seg <- aggregate_test_seg()
  r <- terra::rast(seg)
  terra::values(r) <- c(2,4,10,6,20,30)
  indices <- terra::rowColFromCell(r,seq_len(terra::ncell(r)))
  samples <- coord_from_indices(indices[,2],indices[,1])
  samples$dn <- terra::values(r,mat=FALSE)
  expect_equal(terra::values(aggregate_by_segment(samples,seg)),
    terra::values(compute_mean(r,seg)),ignore_attr=TRUE)
})

test_that("panoramic cells aggregate all their observations with no interpolation", {
  directions <- data.frame(zenith=c(5,6,15,15)*pi/180,
    azimuth=c(5,6,15,16)*pi/180,dn=c(2,8,20,40))
  pano <- to_panoramic(directions,angle_width=10,padding_deg=0)
  averaged <- aggregate_by_segment(pano$sampling_coord,pano$raster)
  values <- terra::values(averaged,mat=FALSE)
  expect_equal(values[c(1,38)],c(5,30))
  expect_equal(sum(!is.na(values)),2)
  peaks <- aggregate_by_segment(pano$sampling_coord,pano$raster,fun=quantile,probs=.99)
  expect_equal(terra::values(peaks,mat=FALSE)[c(1,38)],c(7.94,39.8))
  expect_true("aggregate_by_segment" %in% getNamespaceExports("rcaiman"))
})
