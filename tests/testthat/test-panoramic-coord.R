panoramic_test_d <- function() {
  geometry <- list(model="radial_projection",parameters=equidistant_coef,
    center_of_vision=c(20,15),normalization_radius=15,normalization_angle=90)
  direction_image(set_internal_geometry(c(40,30),geometry),elevation=32.5)
}

test_that("panoramic coordinates preserve observations and continue across zenith", {
  samples <- data.frame(zenith=c(5,30,60,85)*pi/180,
    azimuth=c(5,355,150,200)*pi/180,dn=1:4,x=99,y_down=99,
    col=100,row=100,y=99,relative_radius=.5,image_angle=0)
  pano <- to_panoramic(samples,10,10)
  coord <- pano$sampling_coord
  expect_named(pano,c("sampling_coord","raster"))
  expect_named(coord,c("x","y_down","zenith","azimuth","dn"))
  expect_equal(coord$x,c(1.5,36.5,16,21,19.5,37.5,.5))
  expect_equal(coord$y_down,c(1.5,4,7,9.5,.5,1.5,4))
  source <- c(1,2,3,4,1,1,2)
  expect_equal(coord$dn,samples$dn[source])
  expect_equal(coord$zenith,samples$zenith[source])
  expect_equal(coord$azimuth,samples$azimuth[source])
  # Padding coordinates and the source angles identify the same unit vector.
  theta <- (coord$y_down*10-10)*pi/180
  azimuth <- (coord$x*10-10)*pi/180
  xyz <- function(z,a) cbind(sin(z)*sin(a),sin(z)*cos(a),cos(z))
  expect_equal(xyz(theta,azimuth),xyz(coord$zenith,coord$azimuth))
  expect_equal(dim(pano$raster),c(10,38,1))
  expect_equal(unname(as.vector(terra::ext(pano$raster))),c(0,38,0,10))
})

test_that("zenith padding also generates seam copies at the upper corners", {
  samples <- data.frame(zenith=5*pi/180,azimuth=185*pi/180,dn=7)
  coord <- to_panoramic(samples,5,10)$sampling_coord
  expect_equal(coord$x,c(39,3,75))
  expect_equal(coord$y_down,c(3,1,1))
  expect_equal(coord$dn,rep(7,3))
})

test_that("unpadded coordinates are continuous and missing rays stay missing", {
  samples <- data.frame(zenith=c(0,22.5,90,100,NA,0)*pi/180,
    azimuth=c(0,12.5,359.9,20,30,NA)*pi/180,id=1:6)
  coord <- to_panoramic(samples,5,0)$sampling_coord
  expect_equal(coord$x[1:3],c(0,2.5,71.98))
  expect_equal(coord$y_down[1:3],c(0,4.5,18))
  expect_true(all(is.na(coord$x[4:6]) & is.na(coord$y_down[4:6])))
  expect_identical(coord$id,samples$id)
  empty <- to_panoramic(samples[FALSE,],5,0)
  expect_equal(nrow(empty$sampling_coord),0)
  expect_named(empty$sampling_coord,c("x","y_down","zenith","azimuth","id"))
  expect_equal(nrow(to_panoramic(samples[4:6,],5,10)$sampling_coord),3)
})

test_that("padding removal preserves the central grid and multiple bands", {
  pano <- to_panoramic(data.frame(zenith=numeric(),azimuth=numeric()),5,10)$raster
  pano <- c(pano,pano+1000)
  names(pano) <- c("Blue","Red")
  out <- rem_panoramic_padding(pano,5,10)
  expected <- terra::values(pano,mat=TRUE)[terra::cellFromRowColCombine(pano,3:20,3:74),,drop=FALSE]
  expect_equal(terra::values(out),expected)
  expect_equal(dim(out),c(18,72,2))
  expect_identical(names(out),c("Blue","Red"))
  expect_equal(unname(as.vector(terra::ext(out))),c(0,72,0,18))
  expect_equal(terra::values(rem_panoramic_padding(out,5,0)),terra::values(out))
})

test_that("panoramic return selects containing cells on an inclined rectangular grid", {
  d <- panoramic_test_d()
  for (width in c(5,10)) {
    pano <- to_panoramic(data.frame(zenith=0,azimuth=0),width,0)$raster
    pano <- c(pano,pano+1000)
    names(pano) <- c("Blue","Red")
    out <- from_panoramic(pano,d,width)
    zen <- terra::values(d$zenith,mat=FALSE)
    azi <- terra::values(d$azimuth,mat=FALSE)
    valid <- !is.na(zen) & !is.na(azi) & zen <= 90
    xy <- cbind(azi[valid]/width,nrow(pano)-zen[valid]/width)
    cells <- terra::cellFromXY(pano,xy)
    expect_equal(terra::values(out)[valid,],cbind(Blue=cells,Red=cells+1000))
    expect_true(all(is.na(terra::values(out)[!valid,])))
    expect_true(terra::compareGeom(out,d))
    expect_identical(names(out),names(pano))
  }
})

test_that("return includes the horizon but does not clamp missing or lower-hemisphere rays", {
  d <- terra::rast(nrows=2,ncols=3,nlyrs=2)
  terra::values(d) <- cbind(c(0,90,90.1,0,NA,45),c(0,359.9,20,NA,10,90))
  names(d) <- c("zenith","azimuth")
  pano <- to_panoramic(data.frame(zenith=0,azimuth=0),10,0)$raster
  out <- from_panoramic(pano,d,10)
  expect_equal(as.numeric(terra::values(out)),c(1,324,NA,NA,NA,154))
  d[] <- NA
  expect_true(all(is.na(terra::values(from_panoramic(pano,d,10)))))
})

test_that("the public panoramic workflow composes with add_directions and interpolation", {
  d <- panoramic_test_d()
  indices <- expand.grid(col=seq(1,ncol(d),2),row=seq(1,nrow(d),2))
  m <- !is.na(d$zenith) & !is.na(d$azimuth) & d$zenith <= 90
  indices <- rem_unselected(indices,m)
  samples <- add_directions(indices,d=d)
  samples$dn <- 7
  pano <- to_panoramic(samples,5,10)
  fitted <- interpolate_planar(pano$sampling_coord,r=pano$raster,
    search_radius=100,parallel=FALSE)
  out <- from_panoramic(rem_panoramic_padding(fitted,5,10),d,5)
  selected <- terra::values(m,mat=FALSE)
  expect_equal(as.numeric(terra::values(out))[selected],rep(7,sum(selected)))
  expect_true(all(is.na(terra::values(out)[!selected])))
  expect_true(all(c("to_panoramic","rem_panoramic_padding","from_panoramic") %in%
    getNamespaceExports("rcaiman")))
  expect_false("fisheye_to_pano" %in% getNamespaceExports("rcaiman"))
  expect_false(exists("fisheye_to_pano",envir=asNamespace("rcaiman"),inherits=FALSE))
})

test_that("panoramic grid arguments and angular units are checked", {
  samples <- data.frame(zenith=pi/4,azimuth=pi)
  expect_error(to_panoramic(samples,7,0),"whole number")
  expect_error(to_panoramic(samples,5,12),"whole number")
  expect_error(to_panoramic(samples,1,91),"90")
  expect_error(to_panoramic(data.frame(zenith=40,azimuth=20)),"radians")
  pano <- to_panoramic(samples,5,10)$raster
  expect_error(from_panoramic(pano,panoramic_test_d(),5),"dimensions")
  expect_error(rem_panoramic_padding(pano,5,5),"dimensions")
  blank <- terra::rast(to_panoramic(samples,5,0)$raster)
  expect_error(from_panoramic(blank,panoramic_test_d(),5),"cell values")
})
