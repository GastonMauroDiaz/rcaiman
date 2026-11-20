radial_coord_geometry <- function() {
  list(model="radial_projection",parameters=2/pi,center_of_vision=c(20,15),
       normalization_radius=30,normalization_angle=90)
}

test_that("clock azimuth places directions on the four image axes", {
  geometry <- radial_coord_geometry()
  directions <- data.frame(off_axis=rep(pi/4,4),internal_azimuth=(0:3)*pi/2,
                            id=letters[1:4],row.names=letters[4:1])
  before <- directions
  result <- coord_from_internal_directions(directions,geometry)
  expect_equal(result$x,c(20,35,20,5),tolerance=1e-12)
  expect_equal(result$y_down,c(0,15,30,15),tolerance=1e-12)
  expect_identical(result[names(directions)],directions)
  expect_identical(directions,before)
  back <- internal_directions_from_coord(result,geometry)
  expect_equal(back$off_axis,directions$off_axis,tolerance=1e-11)
  expect_equal(back$internal_azimuth,directions$internal_azimuth,tolerance=1e-12)
})

test_that("subpixel model evaluation agrees with internal raster geometry", {
  geometry <- list(model="radial_projection",parameters=c(0,1/(pi/3)^2),
            center_of_vision=c(3,2),normalization_radius=3.5,
            normalization_angle=60,max_off_axis_angle=140)
  g <- set_internal_geometry(c(8,6),geometry)
  indices <- expand.grid(col=1:8,row=1:6)
  indices <- rem_unselected(indices,!is.na(g$off_axis))
  coord <- coord_from_indices(indices$col,indices$row)
  calculated <- internal_directions_from_coord(coord,geometry)
  extracted <- add_directions(indices,g = g)
  expect_equal(calculated[c("off_axis","internal_azimuth")],
                extracted[c("off_axis","internal_azimuth")],ignore_attr=TRUE,tolerance=1e-12)
  back <- coord_from_internal_directions(calculated,geometry)
  expect_equal(back$x,coord$x,tolerance=1e-8)
  expect_equal(back$y_down,coord$y_down,tolerance=1e-8)
  fractional <- data.frame(x=c(-1.37,3.123,7.721),y_down=c(.832,2.145,-.317))
  back <- coord_from_internal_directions(internal_directions_from_coord(fractional,geometry),geometry)
  expect_equal(back[c("x","y_down")],fractional,tolerance=1e-8)
  class(geometry) <- "geometry_spec"
  geometry$notes <- "Additional registry fields are accepted."
  expect_identical(internal_directions_from_coord(coord,geometry),calculated)
})

test_that("the coordinate path ignores external angles and drops obsolete spatial columns", {
  geometry <- radial_coord_geometry()
  directions <- data.frame(off_axis=c(.3,.3,.5),internal_azimuth=c(.2,.2,2),
    zenith=c(NA,1,2),azimuth=c(NA,2,3),x=-999,y_down=-999,
    col=999,row=999,y=999,relative_radius=999,image_angle=999,
    id=factor(c("a","b","c")),note=c(NA,"kept","kept"),row.names=c("z","y","x"))
  result <- coord_from_internal_directions(directions,geometry)
  expect_false(any(c("col","row","y","relative_radius","image_angle") %in% names(result)))
  expect_identical(result[c("zenith","azimuth","id","note")],
                   directions[c("zenith","azimuth","id","note")])
  expect_identical(rownames(result),rownames(directions))
  expect_identical(result$x[1],result$x[2])
  expect_identical(result$y_down[1],result$y_down[2])
  coord <- result
  coord$off_axis <- coord$internal_azimuth <- -999
  before <- coord
  evaluated <- internal_directions_from_coord(coord,geometry)
  expect_equal(evaluated$off_axis,directions$off_axis,tolerance=1e-11)
  expect_identical(evaluated[setdiff(names(coord),c("off_axis","internal_azimuth"))],
                   coord[setdiff(names(coord),c("off_axis","internal_azimuth"))])
  expect_identical(coord,before)
})

test_that("axis, angular limits, missing values and empty tables have explicit behavior", {
  geometry <- radial_coord_geometry()
  q <- data.frame(off_axis=c(0,pi/2,pi/2+.01,NA,.2,0),
                  internal_azimuth=c(NA,pi/2,pi/2,0,NA,2))
  out <- coord_from_internal_directions(q,geometry)
  expect_equal(out$x,c(20,50,NA,NA,NA,20))
  expect_equal(out$y_down,c(15,15,NA,NA,NA,15))
  coord <- data.frame(x=c(20,50,50.1,NA,20),y_down=c(15,15,15,15,NA))
  out <- internal_directions_from_coord(coord,geometry)
  expect_equal(out$off_axis,c(0,pi/2,NA,NA,NA))
  expect_equal(out$internal_azimuth,c(NA,pi/2,NA,NA,NA))
  geometry$max_off_axis_angle <- 180
  out <- coord_from_internal_directions(data.frame(off_axis=pi,internal_azimuth=1),geometry)
  expect_true(is.na(out$x) && is.na(out$y_down))
  expect_equal(nrow(coord_from_internal_directions(q[FALSE,],geometry)),0L)
  expect_equal(nrow(internal_directions_from_coord(coord[FALSE,],geometry)),0L)
  expect_equal(names(coord_from_internal_directions(q[FALSE,],geometry)),c("x","y_down",names(q)))
})

test_that("one Fibonacci table reaches indices in both calibrated image grids", {
  geometry <- radial_coord_geometry()
  geometry$normalization_radius <- 25
  uncal_g <- set_internal_geometry(c(40,30),geometry)
  other <- geometry
  other$center_of_vision <- c(15,20)
  cal_g <- set_internal_geometry(c(30,40),other)
  common_fov <- intersect_fov(extract_fov(uncal_g),
                                         extract_fov(cal_g))
  sampling_directions <- fibonacci_directions(2)
  sampling_directions <- select_directions_in_fov(sampling_directions,common_fov)
  expect_gt(nrow(sampling_directions),0)
  for (item in list(list(geometry=geometry,g=uncal_g),list(geometry=other,g=cal_g))) {
    coord <- coord_from_internal_directions(sampling_directions,item$geometry)
    indices <- indices_from_coord(coord$x,coord$y_down)
    expect_equal(nrow(indices),nrow(sampling_directions))
    expect_true(all(indices$col >= 1 & indices$col <= ncol(item$g)))
    expect_true(all(indices$row >= 1 & indices$row <= nrow(item$g)))
    expect_gt(sum(duplicated(indices)),0)
    with_indices <- add_indices(coord)
    expect_identical(with_indices[c("col","row")],indices,ignore_attr=TRUE)
    recovered <- internal_directions_from_coord(coord,item$geometry)
    expect_equal(recovered$off_axis,sampling_directions$off_axis,tolerance=1e-11)
    expect_equal(recovered$internal_azimuth,sampling_directions$internal_azimuth,tolerance=1e-12)
  }
})

test_that("coordinate and angle contracts and radial-only inverse are enforced", {
  geometry <- radial_coord_geometry()
  expect_error(coord_from_internal_directions(data.frame(zenith=.5,azimuth=.2),geometry),"off_axis")
  expect_error(coord_from_internal_directions(data.frame(off_axis=90,internal_azimuth=0),geometry),"radians")
  expect_error(coord_from_internal_directions(data.frame(off_axis=.2,internal_azimuth=-1),geometry),"radians")
  expect_error(internal_directions_from_coord(data.frame(x=Inf,y_down=0),geometry),"Inf")
  expect_error(internal_directions_from_coord(data.frame(x=1,y=1),geometry),"y_down")
  geometry$model <- "anisotropic_projection"
  geometry$parameters <- list(off_axis_coef=1,sin_internal_azimuth_coef=c(0,1),
                       cos_internal_azimuth_coef=c(-1,0),atan_scale=.7)
  expect_error(coord_from_internal_directions(data.frame(off_axis=.2,internal_azimuth=0),geometry),"radial_projection")
  expect_error(internal_directions_from_coord(data.frame(x=1,y_down=1),geometry),"dim")
})
