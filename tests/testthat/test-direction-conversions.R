direction_vectors <- function(polar, azimuth) {
  axial <- !is.na(polar) & polar %in% c(0, pi)
  azimuth[axial] <- 0
  cbind(sin(polar)*sin(azimuth), sin(polar)*cos(azimuth), cos(polar))
}

test_that("default orientation relates the two azimuth conventions", {
  q <- data.frame(off_axis=rep(pi/4,4),internal_azimuth=(0:3)*pi/2)
  external <- external_directions_from_internal(q)
  expect_equal(external$zenith,q$off_axis)
  expect_equal(external$azimuth,c(0,3*pi/2,pi,pi/2))
  back <- internal_directions_from_external(external)
  expect_equal(back[c("off_axis","internal_azimuth")],q)
})

test_that("arbitrary rotations round trip and preserve angular separations", {
  q <- data.frame(off_axis=c(0,.01,.5,1.3,2.7,pi),
                  internal_azimuth=c(NA,.2,2,4,6,NA))
  xyz <- direction_vectors(q$off_axis,q$internal_azimuth)
  for (orientation in list(c(0,0,90),c(90,60,32.5),c(-37,721,-23),c(270,180,-90))) {
    external <- external_directions_from_internal(q,orientation[1],orientation[2],orientation[3])
    back <- internal_directions_from_external(external,orientation[1],orientation[2],orientation[3])
    expect_equal(back[c("off_axis","internal_azimuth")],q,tolerance=1e-12)
    e <- direction_vectors(external$zenith,external$azimuth)
    expect_equal(tcrossprod(e),tcrossprod(xyz),tolerance=1e-12)
    external_again <- external_directions_from_internal(back,orientation[1],orientation[2],orientation[3])
    expect_equal(direction_vectors(external_again$zenith,external_again$azimuth),e,tolerance=1e-12)
  }
})

test_that("geographic poles, axial rays and incomplete directions are explicit", {
  external <- data.frame(zenith=c(0,pi,NA,.4),azimuth=c(NA,NA,0,NA))
  internal <- internal_directions_from_external(external,elevation=0)
  expect_equal(internal$off_axis,c(pi/2,pi/2,NA,NA))
  expect_equal(internal$internal_azimuth,c(0,pi,NA,NA))
  back <- external_directions_from_internal(internal,elevation=0)
  expect_equal(back$zenith,c(0,pi,NA,NA))
  expect_true(all(is.na(back$azimuth)))
  q <- data.frame(off_axis=c(0,pi),internal_azimuth=c(NA,NA))
  turned <- external_directions_from_internal(q,73,90,32.5)
  expect_equal(turned$zenith,c(57.5,122.5)*pi/180)
  expect_equal(turned$azimuth,c(270,90)*pi/180)
})

test_that("conversions preserve table metadata and overwrite only destination angles", {
  q <- data.frame(off_axis=c(.4,.4,.8),internal_azimuth=c(2,2,4),
    zenith=-1,azimuth=-1,x=c(1,1,2),row=5,id=factor(c('a','a','b')),
    row.names=c('first','repeat','last'))
  original <- q
  external <- external_directions_from_internal(q,azimuth_rotation=20)
  expect_identical(external[setdiff(names(q),c('zenith','azimuth'))],
                   q[setdiff(names(q),c('zenith','azimuth'))])
  expect_identical(q,original)
  external$off_axis <- external$internal_azimuth <- -999
  internal <- internal_directions_from_external(external,azimuth_rotation=20)
  expect_equal(internal[c('off_axis','internal_azimuth')],q[c('off_axis','internal_azimuth')])
  expect_identical(internal[setdiff(names(external),c('off_axis','internal_azimuth'))],
                   external[setdiff(names(external),c('off_axis','internal_azimuth'))])
  expect_equal(nrow(external_directions_from_internal(q[FALSE,])),0)
  expect_equal(nrow(internal_directions_from_external(external[FALSE,])),0)
})

test_that("tables agree with direction images and reach continuous image coordinates", {
  geometry <- list(model='radial_projection',parameters=2/pi,
    center_of_vision=c(20,15),normalization_radius=25,normalization_angle=90)
  g <- set_internal_geometry(c(40,30),geometry)
  d <- direction_image(g,37,123,32.5)
  indices <- data.frame(col=c(10,15,20,30),row=c(12,18,22,15))
  q <- add_directions(indices,g=g)
  external <- external_directions_from_internal(q,37,123,32.5)
  extracted <- add_directions(indices,d=d)
  expect_equal(external[c('zenith','azimuth')],extracted[c('zenith','azimuth')],tolerance=1e-12)
  internal <- internal_directions_from_external(extracted,37,123,32.5)
  coord <- coord_from_internal_directions(internal,geometry)
  expected <- coord_from_indices(indices$col,indices$row)
  expect_equal(coord[c('x','y_down')],expected,tolerance=1e-10)
  direct <- reproject_points(extracted,geometry,37,123,32.5)
  expect_equal(coord[c('x','y_down')],direct[c('x','y_down')])
  # Reprojection retains source internal columns rather than destination angles.
  external$off_axis <- -99
  expect_identical(reproject_points(external,geometry,37,123,32.5)$off_axis,external$off_axis)
})

test_that("direction conversion validates units and keeps the full sphere", {
  expect_error(internal_directions_from_external(data.frame(zenith=45,azimuth=0)),'radians')
  expect_error(external_directions_from_internal(data.frame(off_axis=.4,internal_azimuth=2*pi)),'radians')
  expect_error(internal_directions_from_external(data.frame(zenith=0,azimuth=0),elevation=91),'elevation')
  expect_error(external_directions_from_internal(data.frame(off_axis=.4,internal_azimuth=0),sensor_rotation=Inf),'Inf')
  q <- data.frame(zenith=2.5,azimuth=.3)
  expect_equal(internal_directions_from_external(q)$off_axis,2.5)
})
