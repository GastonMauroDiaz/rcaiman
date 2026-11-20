hsp_fixture <- function() test_path("fixtures", "hsp", "MA_2021_11_18")

test_that("real HSP manual inputs retain settings, marks and zenith references", {
  cases <- list(
    DSCN6280 = list(settings = c(60, 10, 10, 60), sun = c(531, 236),
      zenith_dn = 3921.23, n = 141L, first = c(384, 385, 3),
      last = c(612, 342, 3)))
  for (id in names(cases)) {
    expected <- cases[[id]]
    input <- hsp_read_manual_input(hsp_fixture(), id)
    expect_equal(unlist(input[c("weight", "max_points", "angle", "point_radius")],
      use.names = FALSE), expected$settings)
    expect_s3_class(input$sun_indices, "data.frame")
    expect_equal(input$sun_indices, data.frame(row = expected$sun[1], col = expected$sun[2]))
    expect_equal(input$zenith_dn, expected$zenith_dn)
    expect_equal(nrow(input$sky_indices), expected$n)
    expect_named(input$sky_indices, c("col", "row", "type"))
    expect_equal(unname(unlist(input$sky_indices[1, ])), expected$first)
    expect_equal(unname(unlist(input$sky_indices[expected$n, ])), expected$last)
    with_ext <- hsp_read_manual_input(hsp_fixture(), paste0(id, ".pgm"))
    expect_equal(with_ext$sky_indices, input$sky_indices)
    expect_equal(with_ext$sun_indices, input$sun_indices)
  }
})

test_that("HSP coefficients are read in CIE order with decimal commas", {
  expected <- c(2.69285557, -1.99973217, 13.16075166, -.0176585, .00000001)
  expect_equal(hsp_read_opt_sky_coef(hsp_fixture(), "DSCN6280"), expected)
  expect_equal(hsp_read_opt_sky_coef(hsp_fixture(), "DSCN6280.pgm"), expected)
})

test_that("both polygon types preserve the saved vertices on the image grid", {
  input <- hsp_read_manual_input(hsp_fixture(), "DSCN6280.pgm")
  r <- read_caim(file.path(hsp_fixture(), "manipulate", "DSCN6280.pgm"))
  expect_equal(dim(r), c(756, 756, 1))
  expected <- list(
    masks = rbind(c(476,254), c(467,245), c(468,260)),
    marked_areas = rbind(c(399,345), c(378,358), c(388,379), c(404,370),
      c(366,213), c(354,204), c(348,206), c(350,215), c(356,215)))
  counts <- c(masks = 1, marked_areas = 2)
  sorted_xy <- function(x) unname(x[order(x[, 1], x[, 2]), , drop = FALSE])
  for (field in names(expected)) {
    polygons <- input[[field]]
    expect_s4_class(polygons, "SpatVector")
    expect_equal(nrow(polygons), unname(counts[field]))
    expect_equal(terra::crs(polygons), terra::crs(r))
    vertices <- unique(terra::geom(polygons)[, c("x", "y"), drop = FALSE])
    # HSP indices map to pixel centers in the standardized image grid.
    xy <- cbind(expected[[field]][, 1] - .5,
                756 - expected[[field]][, 2] + .5)
    expect_equal(sorted_xy(vertices), sorted_xy(xy))
    expect_true(all(terra::expanse(polygons, transform = FALSE) > 0))
  }
})

test_that("empty masks are read from a temporary version of the Nikon project", {
  project <- tempfile("hsp-empty-masks-")
  dir.create(project)
  expect_true(file.copy(file.path(hsp_fixture(), "manipulate"), project,
                        recursive = TRUE))
  for (suffix in c("masks", "polygons")) {
    writeLines("0", file.path(project, "manipulate",
                             paste0("DSCN6280_", suffix, ".conf")))
  }
  input <- hsp_read_manual_input(project, "DSCN6280")
  expect_null(input$masks)
  expect_null(input$marked_areas)
  original <- hsp_read_manual_input(hsp_fixture(), "DSCN6280")
  expect_equal(nrow(original$masks), 1)
  expect_equal(nrow(original$marked_areas), 2)
})

test_that("HSP writers round-trip points and sun in a temporary project", {
  project <- tempfile("hsp-project-")
  dir.create(project)
  expect_true(file.copy(file.path(hsp_fixture(), "manipulate"), project,
                        recursive = TRUE))
  for (id in "DSCN6280") {
    input <- hsp_read_manual_input(project, id)
    selected <- input$sky_indices[c(1, 5, 10), ]
    selected$type <- 99
    sun <- input$sun_indices
    sun$row <- sun$row + 1
    sun$col <- sun$col - 1
    sun <- sun[, c("col", "row")]
    sun$zenith <- .5
    hsp_write_sky_indices(selected, project, paste0(id, ".pgm"))
    hsp_write_sun_indices(sun, project, paste0(id, ".pgm"))
    result <- hsp_read_manual_input(project, id)
    expect_equal(unname(as.matrix(result$sky_indices[, c("col", "row")])),
                 unname(as.matrix(selected[, c("col", "row")])))
    expect_equal(result$sky_indices$type, rep(3, 3))
    expect_equal(result$sun_indices[, c("col", "row")], sun[, c("col", "row")])
    expect_equal(readLines(file.path(project, "manipulate", paste0(id, "_sun.conf"))),
                 paste(sun$col, sun$row, sep = "."))
  }
  # Reading the original fixture again proves the test did not overwrite it.
  expect_equal(nrow(hsp_read_manual_input(hsp_fixture(), "DSCN6280")$sky_indices), 141)
})

test_that("all HSP helpers validate the image name", {
  expect_error(hsp_read_manual_input(hsp_fixture(), NA_character_), "img_name")
  expect_error(hsp_read_opt_sky_coef(hsp_fixture(), 1), "img_name")
  expect_error(hsp_write_sky_indices(data.frame(col = 1, row = 1),
    hsp_fixture(), NA_character_), "img_name")
  expect_error(hsp_write_sun_indices(data.frame(row = 1, col = 2), hsp_fixture(), 1), "img_name")
})

test_that("HSP writers enforce named numeric index tables", {
  project <- tempfile("hsp-contract-")
  dir.create(file.path(project, "manipulate"), recursive = TRUE)
  bad_tables <- list(c(col = 1, row = 2), matrix(c(1, 2), nrow = 1),
                     data.frame(col = 1))
  for (value in list(0, -1, 1.5, NA_real_, Inf, TRUE, "1", factor("1"))) {
    bad_tables[[length(bad_tables) + 1L]] <- data.frame(col = value, row = 1)
  }
  for (writer in list(hsp_write_sky_indices, hsp_write_sun_indices)) {
    for (indices in bad_tables) expect_error(writer(indices, project, "DSCN6280"))
  }
  expect_error(hsp_write_sun_indices(data.frame(col = c(1, 2), row = c(3, 4)),
    project, "DSCN6280"), "length 1")
  expect_error(hsp_write_sun_indices(data.frame(col = numeric(), row = numeric()),
    project, "DSCN6280"), "length 1")
  # Both integer storage modes are accepted; auxiliary columns are ignored.
  expect_invisible(hsp_write_sun_indices(data.frame(col = 236L, row = 531, note = "sun"),
    project, "DSCN6280"))
  expect_identical(readLines(file.path(project, "manipulate", "DSCN6280_sun.conf")),
                   "236.531")
  expect_invisible(hsp_write_sky_indices(data.frame(row = 531L, col = 236, note = "sky"),
    project, "DSCN6280"))
})
