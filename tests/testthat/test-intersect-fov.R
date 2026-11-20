coverage_plane <- function(x, y, rotation = diag(3)) {
  grid <- expand.grid(x = x, y = rev(y))
  xyz <- cbind(grid$x, grid$y, 1)
  xyz <- xyz / sqrt(rowSums(xyz^2))
  xyz <- xyz %*% t(rotation)
  angles <- cbind(acos(pmax(-1, pmin(1, xyz[, 3]))) * 180 / pi,
                  (atan2(xyz[, 1], xyz[, 2]) * 180 / pi) %% 360)
  r <- terra::rast(nrows = length(y), ncols = length(x), nlyrs = 2,
                   xmin = 0, xmax = length(x), ymin = 0, ymax = length(y), crs = "")
  terra::values(r) <- angles
  names(r) <- c("off_axis", "internal_azimuth")
  r
}

plane_fov <- function(xy) {
  # xy vertices run clockwise in the image plane (y points upwards).
  r <- sqrt(rowSums(xy^2))
  data.frame(off_axis=atan(r), internal_azimuth=atan2(xy[,1],xy[,2]) %% (2*pi))
}

rectangle_fov <- function(xmin, xmax, ymin, ymax) {
  plane_fov(rbind(c(xmin,ymax),c(xmax,ymax),c(xmax,ymin),c(xmin,ymin),c(xmin,ymax)))
}

test_that("contour intersection and raster representation match known windows", {
  x <- seq(-1,1,length.out=21)
  y <- seq(-.8,.8,length.out=17)
  g <- coverage_plane(x,y)
  before <- terra::values(g)
  first <- extract_fov(g)
  second <- rectangle_fov(-.6,.8,-.4,.5)
  common <- intersect_fov(first,second)
  result <- clip_fov(g,common)
  xy <- expand.grid(x=x,y=rev(y))
  selected <- with(xy,x >= -.6-1e-12 & x <= .8+1e-12 & y >= -.4-1e-12 & y <= .5+1e-12)
  expected <- before
  expected[!selected,] <- NA_real_
  expect_identical(terra::values(result),expected)
  expect_identical(terra::values(g),before)
  expect_true(terra::compareGeom(result,g))
  expect_identical(names(result),names(g))
  expect_identical(names(common),c("off_axis","internal_azimuth"))
  expect_equal(as.numeric(common[1,]),as.numeric(common[nrow(common),]))
  shifted <- terra::shift(g,dx=35,dy=-9)
  expect_identical(terra::values(clip_fov(shifted,common)),expected)
  queries <- as.data.frame(before*pi/180)
  expect_identical(select_directions_in_fov(queries,common),queries[selected,])
})

test_that("crossing arcs add vertices and either argument order gives the same area", {
  a <- rectangle_fov(-1,1,-.4,.4)
  b <- rectangle_fov(-.4,.4,-1,1)
  common <- intersect_fov(a,b)
  expect_equal(nrow(common),5L)
  xy <- cbind(tan(common$off_axis)*sin(common$internal_azimuth),
              tan(common$off_axis)*cos(common$internal_azimuth))
  expect_equal(abs(xy),matrix(.4,5,2),tolerance=1e-12)
  q <- fibonacci_directions(5)
  q$id <- seq_len(nrow(q))
  sequential <- select_directions_in_fov(select_directions_in_fov(q,a),b)
  expect_identical(select_directions_in_fov(q,common),sequential)
  expect_identical(select_directions_in_fov(q,intersect_fov(b,a)),sequential)
  expect_identical(select_directions_in_fov(q,intersect_fov(a,a)),
                   select_directions_in_fov(q,a))
})

test_that("concavity is retained without filling an indentation", {
  concave <- plane_fov(rbind(c(-1,1),c(0,1),c(0,0),c(1,0),c(1,-1),c(-1,-1),c(-1,1)))
  common <- intersect_fov(concave,rectangle_fov(-2,2,-2,2))
  queries <- plane_fov(rbind(c(-.5,.5),c(.5,.5),c(.5,-.5),c(0,0)))
  expect_identical(select_directions_in_fov(queries,common),queries[c(1,3,4),])
})

test_that("azimuth wrapping and fields larger than a hemisphere retain their interior", {
  cap <- function(angle) data.frame(off_axis=rep(angle*pi/180,5),
                                    internal_azimuth=c(0,pi/2,pi,3*pi/2,0))
  large <- cap(110)
  common <- intersect_fov(large,cap(100))
  q <- data.frame(off_axis=c(0,95,105,180)*pi/180, internal_azimuth=0)
  expect_identical(select_directions_in_fov(q,common),q[1:2,])
  small <- cap(60)
  seam <- data.frame(off_axis=rep(30*pi/180,3),internal_azimuth=c(1e-5,2*pi-1e-5,pi))
  expect_identical(select_directions_in_fov(seam,small),seam)
  # The order of a spherical loop distinguishes its interior from its complement.
  poles <- data.frame(off_axis=c(0,pi),internal_azimuth=0)
  expect_identical(select_directions_in_fov(poles,small),poles[1,,drop=FALSE])
  expect_identical(select_directions_in_fov(poles,small[5:1,]),poles[2,,drop=FALSE])
})

test_that("empty and zero-area intersections stay empty and multiple contours fail", {
  a <- rectangle_fov(-1,0,-1,1)
  b <- rectangle_fov(.1,1,-1,1)
  empty <- data.frame(off_axis=numeric(),internal_azimuth=numeric())
  expect_identical(intersect_fov(a,b),empty)
  expect_identical(intersect_fov(a,rectangle_fov(0,1,-1,1)),empty)
  expect_identical(intersect_fov(a,empty),empty)
  expect_identical(intersect_fov(empty,a),empty)
  # A U-shaped region cut across its arms yields two disconnected areas.
  u <- plane_fov(rbind(c(-2,2),c(-1,2),c(-1,-1),c(1,-1),c(1,2),c(2,2),
                        c(2,-2),c(-2,-2),c(-2,2)))
  expect_error(intersect_fov(u,rectangle_fov(-3,3,0,1)),"more than one contour")
  g <- coverage_plane(c(-1,0,1),c(-1,0,1))
  expect_true(all(is.na(terra::values(clip_fov(g,empty)))))
})

test_that("direction filtering preserves metadata, missing-value behavior and boundaries", {
  fov <- rectangle_fov(-1,1,-1,1)
  q <- data.frame(off_axis=c(0,atan(1),atan(1.01),0,NA,.2),
    internal_azimuth=c(0,pi/2,pi/2,0,0,NA),col=999,row=-.5,
    note=factor(c("a","b","c","a","b","c")),other=NA_real_,row.names=letters[1:6])
  before <- q
  expect_identical(select_directions_in_fov(q,fov),q[c(1,2,4),])
  expect_identical(q,before)
  expect_identical(select_directions_in_fov(q[FALSE,],fov),q[FALSE,])
  expect_identical(select_directions_in_fov(q,fov[FALSE,]),q[FALSE,])
  expect_identical(select_directions_in_fov(fov,fov),fov)
  g <- coverage_plane(c(-.5,0,.5),c(-.5,0,.5))
  values <- terra::values(g)
  values[1,1] <- NA_real_
  values[2,2] <- NA_real_
  terra::values(g) <- values
  result <- terra::values(clip_fov(g,fov))
  expect_true(all(is.na(result[1:2,])))
  expect_identical(result[3:9,],values[3:9,])
  expect_identical(terra::values(g),values)
})

test_that("extracted radial and anisotropic contours retain their own pixel centers", {
  radial <- set_internal_geometry(c(40,30),list(model="radial_projection",
    parameters=equidistant_coef,center_of_vision=c(20,15),normalization_radius=25,
    normalization_angle=90))
  anisotropic <- set_internal_geometry(c(20,26),list(model="anisotropic_projection",
    center_of_vision=c(10,13),parameters=list(off_axis_coef=c(1.1,.2),
      sin_internal_azimuth_coef=c(0,1),cos_internal_azimuth_coef=c(-1,0),atan_scale=.7)))
  for (g in list(radial,anisotropic)) {
    fov <- extract_fov(g)
    expect_identical(terra::values(clip_fov(g,fov)),terra::values(g))
  }
  geometry <- list(model="radial_projection",parameters=equidistant_coef,
            center_of_vision=c(15,20),normalization_radius=25,normalization_angle=90)
  other <- set_internal_geometry(c(30,40),geometry)
  common <- intersect_fov(extract_fov(radial),extract_fov(other))
  expect_equal(sum(complete.cases(terra::values(clip_fov(radial,common)))),900)
  expect_equal(sum(complete.cases(terra::values(clip_fov(other,common)))),900)
})

test_that("angular units and invalid contours are checked without silent repair", {
  fov <- rectangle_fov(-1,1,-1,1)
  bad <- fov
  bad$off_axis[1] <- 90
  expect_error(intersect_fov(bad,fov),"radians")
  bad$off_axis[1] <- NA_real_
  expect_error(intersect_fov(bad,fov),"NA")
  bad$off_axis[1] <- Inf
  expect_error(select_directions_in_fov(bad,fov),"Inf")
  expect_error(intersect_fov(fov[1,,drop=FALSE],fov),"three vertices")
  bow <- fov[c(1,3,2,4,1),]
  expect_error(intersect_fov(bow,fov),"Invalid contour")
  g <- coverage_plane(c(-1,0,1),c(-1,0,1))
  expect_error(clip_fov(g[[1]],fov),"layers")
  expect_error(select_directions_in_fov(fov["off_axis"],fov),"internal_azimuth")
})

