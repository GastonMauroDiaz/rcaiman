internal_geometry_specs <- function() {
  list(
    list(model = "radial_projection", center_of_vision = c(3, 2),
         parameters = 2 / pi, normalization_radius = 3.5,
         normalization_angle = 90),
    list(model = "anisotropic_projection", center_of_vision = c(3, 2),
         parameters = list(off_axis_coef = 1,
                           sin_internal_azimuth_coef = c(0, 2),
                           cos_internal_azimuth_coef = c(-1, 0),
                           atan_scale = 0.6))
  )
}

test_that("geometry lists accept classes and metadata without mutation", {
  for (geometry in internal_geometry_specs()) {
    spec <- geometry
    spec$id <- "example"
    spec$date <- as.Date("2026-09-24")
    spec$notes <- "Additional registry fields are not model inputs."
    spec$radiometry <- list(example = list(parameters = c(1, 2)))
    class(spec) <- "geometry_spec"
    before <- serialize(spec, NULL)

    expect_true(.check_internal_geometry(spec))
    ordinary <- set_internal_geometry(c(8, 6), geometry)
    registered <- set_internal_geometry(c(8, 6), spec)
    expect_identical(terra::values(registered), terra::values(ordinary))
    expect_identical(names(registered), names(ordinary))
    expect_true(terra::compareGeom(registered, ordinary))
    expect_identical(serialize(spec, NULL), before)
    expect_null(attr(registered, "geometry"))
  }
})

test_that("an omitted radial limit means 90 without changing the list", {
  geometry <- internal_geometry_specs()[[1]]
  before <- geometry
  default <- set_internal_geometry(c(8, 6), geometry)
  explicit <- geometry
  explicit$max_off_axis_angle <- 90
  expect_identical(terra::values(default),
                   terra::values(set_internal_geometry(c(8, 6), explicit)))
  expect_identical(geometry, before)
  expect_false("max_off_axis_angle" %in% names(geometry))
})

test_that("metadata names cannot stand in for required geometry fields", {
  geometry <- internal_geometry_specs()[[1]]
  geometry$center_of_vision_original <- geometry$center_of_vision
  geometry$center_of_vision <- NULL
  expect_error(set_internal_geometry(c(8, 6), geometry), "center_of_vision")
})
