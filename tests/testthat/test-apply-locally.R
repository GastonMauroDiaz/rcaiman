local_fixture <- function() {
  list(samples = data.frame(dn = c(10, 20, 30)),
       queries = data.frame(id = c("b", "a", "c"), fallback = c(90, 80, 70),
                            row.names = c("q3", "q1", "q2")),
       neighbors = list(list(index = c(3L, 1L), distance = c(.25, .75)),
                        list(index = integer(0), distance = numeric(0)),
                        list(index = 2L, distance = 0)))
}

test_that("local application needs no geometry and retains query order", {
  f <- local_fixture()
  original <- f
  fun <- function(data, distance, i) {
    expect_s3_class(data, "data.frame")
    expect_identical(names(data), "dn")
    expect_identical(distance, f$neighbors[[i]]$distance)
    if (!nrow(data)) return(f$queries$fallback[i])
    expect_equal(data$dn, f$samples$dn[f$neighbors[[i]]$index])
    mean(data$dn)
  }
  got <- apply_locally(f$samples, f$queries, f$neighbors, fun)
  expect_identical(got[names(f$queries)], f$queries)
  expect_equal(got$output, c(20, 80, 20))
  expect_identical(f, original)
  # A single-column sample table and a repeated neighbor stay a data frame.
  f$neighbors[[1]] <- list(index = c(3, 1, 3), distance = c(0, 2, 0))
  got <- apply_locally(f$samples, f$queries, f$neighbors,
                       function(data, ...) nrow(data))
  expect_equal(got$output, c(3, 0, 1))
})

test_that("empty samples and queries retain their documented structure", {
  f <- local_fixture()
  empty <- rep(list(list(index = integer(), distance = numeric())), 3)
  got <- apply_locally(f$samples[FALSE, , drop = FALSE], f$queries, empty,
                       function(data, distance, i) {
                         expect_equal(nrow(data), 0)
                         expect_identical(distance, numeric())
                         f$queries$fallback[i]
                       })
  expect_equal(got$output, f$queries$fallback)
  got <- apply_locally(f$samples, f$queries[FALSE, ], list(),
                       function(...) stop("should not run"), parallel = TRUE)
  expect_identical(got$output, numeric())
  expect_identical(got[names(f$queries)], f$queries[FALSE, ])
})

test_that("only individually simple scalars are combined into a column", {
  f <- local_fixture()
  run <- function(results) apply_locally(f$samples, f$queries, f$neighbors,
                                         function(data, distance, i) results[[i]])
  expect_identical(run(list(c(a = 1), NA_real_, 2L))$output, c(1, NA_real_, 2))
  expect_identical(run(list(TRUE, NA, FALSE))$output, c(TRUE, NA, FALSE))
  expect_identical(run(list(1, "a", 2))$output, c("1", "a", "2"))
  # Flattening these would produce three elements and lose the empty result.
  results <- list(NULL, c(a = 1, b = 2), 3)
  expect_identical(run(results), results)
  results <- list(numeric(), list(value = 1), matrix(2))
  expect_identical(run(results), results)
  results <- list(as.Date("2026-10-27"), factor("a"), list(value = 3))
  expect_identical(run(results), results)
  expect_identical(run(rep(list(NULL), 3)), rep(list(NULL), 3))
  model <- lm(dn ~ seq_along(dn), data = f$samples)
  results <- list(model, NULL, model)
  expect_identical(run(results), results)
})

test_that("neighborhoods are validated before evaluation", {
  f <- local_fixture()
  run <- function(neighbors) apply_locally(f$samples, f$queries, neighbors,
                                           function(...) stop("callback called"))
  expect_error(run(f$neighbors[-1]), "one element per")
  for (index in list(0, -1, 1.5, NA_integer_, Inf, 4, TRUE)) {
    bad <- f$neighbors
    bad[[1]] <- list(index = index, distance = 0)
    expect_error(run(bad), "index")
  }
  for (distance in list(-1, NA_real_, Inf, "far", c(1, 2))) {
    bad <- f$neighbors
    bad[[1]] <- list(index = 1L, distance = distance)
    expect_error(run(bad), "distance")
  }
  bad <- f$neighbors
  bad[[1]] <- list(index = 1L)
  expect_error(run(bad), "index.*distance")
  expect_error(apply_locally(f$samples, transform(f$queries, output = 0),
                            f$neighbors, mean), "reserved")
  expect_error(apply_locally(1:3, f$queries, f$neighbors, mean), "data.frame")
  expect_error(apply_locally(f$samples, f$queries, f$neighbors, NULL), "function")
})

test_that("callback errors identify the query row including empty neighborhoods", {
  f <- local_fixture()
  expect_error(apply_locally(f$samples, f$queries, f$neighbors,
                             function(data, ...) {
                               if (!nrow(data)) stop("no samples")
                               mean(data$dn)
                             }), "query_data row 2: no samples")
})

test_that("nn_spherical neighborhoods work with both angular conventions", {
  sampling <- sky_grid_centers(15)
  sampling$dn <- seq_len(nrow(sampling))
  query <- sky_grid_centers(30)
  fun <- function(data, ...) if (nrow(data)) mean(data$dn) else NA_real_
  neighbors <- nn_spherical(sampling, query, k = 3, search_radius = 20,
                             parallel = FALSE)
  external <- apply_locally(sampling, query, neighbors, fun)
  names(sampling)[1:2] <- names(query) <- c("off_axis", "internal_azimuth")
  neighbors <- nn_spherical(sampling, query, k = 3, search_radius = 20,
                             angle_columns = names(query), parallel = FALSE)
  internal <- apply_locally(sampling, query, neighbors, fun)
  expect_equal(external$output, internal$output)
})

test_that("parallel evaluation preserves nested results and NULL", {
  skip_on_cran()
  skip_if_not_installed("doParallel")
  skip_if(is.na(parallel::detectCores()) || parallel::detectCores() < 2)
  on.exit(foreach::registerDoSEQ(), add = TRUE)
  f <- local_fixture()
  # The sequential tests cover scalar simplification. Here exercise the harder
  # serialization case: a callback returning nested objects and NULL.
  fun <- function(data, distance, i) {
    if (!nrow(data)) return(NULL)
    list(query = i, samples = data, distance = distance)
  }
  expected <- apply_locally(f$samples, f$queries, f$neighbors, fun)
  actual <- apply_locally(f$samples, f$queries, f$neighbors, fun,
                          parallel = TRUE, cores = 2)
  expect_identical(actual, expected)
})
