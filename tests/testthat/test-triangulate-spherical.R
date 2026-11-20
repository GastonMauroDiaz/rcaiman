spherical_triangle_fixture <- function() {
  samples <- data.frame(
    zenith = c(20, 35, 55, 60, 80, 45, 75, 25) * pi / 180,
    azimuth = c(10, 80, 130, 180, 230, 280, 320, 350) * pi / 180,
    dn = seq(10, 80, 10))
  queries <- data.frame(
    zenith = c(30, 50, 65, 30, 120, 10) * pi / 180,
    azimuth = c(30, 170, 230, 310, 70, 210) * pi / 180,
    label = letters[1:6])
  rownames(queries) <- paste0("query_", 1:6)
  list(samples = samples, queries = queries)
}

test_that("refactoring preserves the previous triangulation estimates", {
  f <- spherical_triangle_fixture()
  result <- triangulate_spherical(f$samples, f$queries, "dn")
  # Values recorded from the preceding implementation on this nondegenerate set.
  expect_equal(result$dn, c(NA, 35.53079421981698, 50.34412984025131,
    71.00413631054465, NA, 26.47665303946306), tolerance = 1e-11)
  expect_identical(result[names(f$queries)], f$queries)
  expect_type(result$triangle_id, "integer")
  expect_identical(is.na(result$triangle_id), is.na(result$dn))
  expect_identical(result, triangulate_spherical(f$samples, f$queries, "dn"))
  expect_identical(result, triangulate_spherical(f$samples, f$queries, 3))
})

test_that("internal angles and extra index columns use the same calculation", {
  f <- spherical_triangle_fixture()
  expected <- triangulate_spherical(f$samples, f$queries, "dn")
  names(f$samples)[1:2] <- names(f$queries)[1:2] <- c("off_axis", "internal_azimuth")
  f$queries$col <- NA_integer_
  f$queries$row <- -1L
  result <- triangulate_spherical(f$samples, f$queries, "dn",
    angle_columns = c("off_axis", "internal_azimuth"))
  expect_identical(result$dn, expected$dn)
  expect_identical(result$triangle_id, expected$triangle_id)
  expect_identical(result[names(f$queries)], f$queries)
})

test_that("vertices and edge queries are repeatable and constant fields are preserved", {
  f <- spherical_triangle_fixture()
  at_vertices <- triangulate_spherical(f$samples, f$samples[1:2], "dn")
  expect_equal(at_vertices$dn, f$samples$dn, tolerance = 1e-10)
  xyz <- as.matrix(.spherical_neighbor_vectors(f$samples, c("zenith", "azimuth")))
  p <- .7 * xyz[6, ] + .3 * xyz[4, ]
  p <- p / sqrt(sum(p^2))
  edge <- as.data.frame(.cartesian_to_directions(matrix(p, nrow = 1)))
  names(edge) <- c("zenith", "azimuth")
  repeated <- edge[rep(1, 4), ]
  result <- triangulate_spherical(f$samples, repeated, "dn")
  expect_true(all(is.finite(result$dn)))
  expect_length(unique(result$dn), 1)
  expect_length(unique(result$triangle_id), 1)
  expect_identical(result, triangulate_spherical(f$samples, repeated, "dn"))
  f$samples$dn <- 7
  constant <- triangulate_spherical(f$samples, f$queries, "dn")
  expect_equal(constant$dn[!is.na(constant$dn)], rep(7, 4), tolerance = 1e-12)
})

test_that("missing observations are excluded and missing queries retain their positions", {
  f <- spherical_triangle_fixture()
  expected <- triangulate_spherical(f$samples, f$queries, "dn")
  extra <- f$samples[1:3, ]
  extra$dn[1] <- NA_real_
  extra$zenith[2] <- NA_real_
  extra$azimuth[3] <- NA_real_
  expect_identical(triangulate_spherical(rbind(f$samples, extra), f$queries, "dn"),
    expected)
  q <- f$queries
  q$zenith[2] <- NA_real_
  q$azimuth[3] <- NA_real_
  result <- triangulate_spherical(f$samples, q, "dn")
  expect_identical(result[names(q)], q)
  expect_true(all(is.na(result$dn[c(2, 3)])))
  expect_true(all(is.na(result$triangle_id[c(2, 3)])))
  expect_equal(result$dn[-c(2, 3)], expected$dn[-c(2, 3)])
})

test_that("a polar direction does not require azimuth", {
  f <- spherical_triangle_fixture()
  q <- data.frame(zenith = c(0, 0), azimuth = c(NA_real_, 0))
  result <- triangulate_spherical(f$samples, q, "dn")
  expect_true(all(is.finite(result$dn)))
  expect_identical(result$dn[1], result$dn[2])
  polar <- data.frame(zenith = 0, azimuth = NA_real_, dn = 100)
  samples <- rbind(f$samples, polar)
  at_pole <- triangulate_spherical(samples, q, "dn")
  expect_equal(at_pole$dn, c(100, 100), tolerance = 1e-10)
})

test_that("empty inputs and insufficient observations preserve output structure", {
  f <- spherical_triangle_fixture()
  empty <- triangulate_spherical(f$samples, f$queries[FALSE, ], "dn")
  expect_identical(empty[names(f$queries)], f$queries[FALSE, ])
  expect_identical(empty$dn, numeric(0))
  expect_identical(empty$triangle_id, integer(0))
  for (n in 0:2) {
    result <- triangulate_spherical(f$samples[seq_len(n), ], f$queries, "dn")
    expect_identical(result[names(f$queries)], f$queries)
    expect_true(all(is.na(result$dn)))
    expect_true(all(is.na(result$triangle_id)))
  }
})

test_that("direction and output-column contracts are explicit", {
  f <- spherical_triangle_fixture()
  q <- f$queries
  q$dn <- 1
  expect_error(triangulate_spherical(f$samples, q, "dn"), "reserved")
  q$dn <- NULL
  q$triangle_id <- 1L
  expect_error(triangulate_spherical(f$samples, q, "dn"), "reserved")
  expect_error(triangulate_spherical(f$samples, f$queries, "zenith"), "must differ")
  expect_error(triangulate_spherical(f$samples, f$queries, "dn",
    angle_columns = c("zenith", "zenith")), "different")
  bad <- f$samples
  bad$dn[1] <- Inf
  expect_error(triangulate_spherical(bad, f$queries, "dn"), "Inf")
  bad <- f$samples
  bad$zenith[1] <- 30
  expect_error(triangulate_spherical(bad, f$queries, "dn"), "radians")
})
