#' Apply a function locally using sampled points on the hemispherical domain
#'
#' @description
#' Evaluate a user-supplied function locally using neighboring points
#' defined on the hemispherical domain.
#'
#' @details
#' The `query_points` argument lets the user evaluate `fun` at a custom set of
#' locations. For each point in `query_points`, neighbors are searched within
#' `sampling_points`. Neighbor selection is based on
#' spherical distance, see [nn_spherical()]. The function `fun` is
#' then applied to the corresponding subset of neighboring sampling points.
#'
#' The input passed to `fun` is a **subset of rows of `sampling_points`**,
#' therefore, `fun` must extract and process the column(s) of interest, e.g.:
#' `function(df) mean(df$dn)`.
#'
#' The angular neighborhood is defined by two parameters:
#' \itemize{
#'   \item `k`: number of nearest neighbors used by `fun`.
#'   \item `search_radius`: maximum allowed spherical distance (deg).
#' }
#'
#' Argument `rule` controls how neighbor availability is interpreted:
#' \itemize{
#'   \item `"any"`: all neighbors found within `search_radius` are passed to `fun`
#'                 (even if fewer than `k`).
#'   \item `"all"`: if fewer than `k` neighbors lie within the radius, the query
#'                 point is assigned `NA`.
#' }
#'
#' @param query_points `data.frame` with columns `"zen"` and `"azi"`, as in the
#'   output of [add_directions]. Points at which `fun` should be evaluated.
#' @param neighbors a list where each element corresponds to one query point and contains :
#' \describe{
#'   \item{index}{
#'   indices of neighboring sampling points
#'   }
#'   \item{distance}{
#'   angular distances in radians
#'   }
#' }
#' Tipically the output of [nn_spherical()].
#' @param fun function. It receives a data frame, an integer vector of length
#'   equal to the number of rows in the data frame, and an integer vector of
#'   length one. The first object is a subset of sampling_points, the second is the
#'   angular distance of each point to the query point, and the third is the
#'   query-point row index. When the search retrieves no points, the data frame
#'   is empty. Users should include the ellipsis (...) in their function
#'   signature to gracefully handle any of these arguments that are not required
#'   for their specific calculation. See *Examples*.
#' @param parallel logical vector of length one. If `TRUE`, operations are
#'   executed in parallel.
#' @param cores numeric vector of length one. Must be positive and interpreted as
#'   integer. Number of CPU cores to use when
#'   `parallel = TRUE`. If `NULL`, the function uses all available
#'   cores detected internally. Values greater than the number of available
#'   cores are silently reduced to the maximum allowed. If the resulting number
#'   of cores is less than 2, `parallel` is silently set to `FALSE`.
#' @param leave_free numeric vector of length one. Must be positive and interpreted as
#'   integer. Number of CPU cores to leave
#'   free when `parallel = TRUE`. Ignored when `cores` is not `NULL`. Must be
#'   greater than or equal to 0. If `total_cores - leave_free` is less than 2,
#'   `parallel` is silently set to `FALSE`.
#' @param logical logical vector of length one. Internally it is passed unchanged
#'   to the argument `logical` in [parallel::detectCores()]. Ignored when
#'   `parallel = FALSE`.
#'
#' @inheritParams nn_spherical
#'
#' @note
#' This function assumes that `sky_points`, `query_points`, and `neighbors`
#' share the same geometry. No checks are performed to enforce this.
#'
#' @seealso
#' [sky_grid_centers()], [fibonacci_points()], and [raster_grid_in_circle()] for
#' producing `query_points`.
#'
#' @return
#' If `fun` returns a vector of length one, a `data.frame` equal to
#' `query_points` with one additional column, named `"output"`, containing the
#' result of applying `fun` at each query point.
#'
#' If `fun` returns any other object, it returns a list with as many components
#' as rows in `query_points`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(rcaiman.registry)
#' path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
#' my_entry <- read_registry_snapshot(path)
#' spec <- get_geometry_spec(my_entry,  "exif_01", "simple_method")
#' caim <- read_caim()
#' r <- caim$Cyan
#' z <- zenith_image(ncol(r), spec$parameters)
#' a <- azimuth_image(z)
#' rspec <- get_radiometry_spec(spec, "simple_method")
#' ff <- ff_simple_polynomial(z, rspec$parameters$`5.0`)
#' r <- apply_flat_field(r, ff)
#' m <- !is.na(z)
#'
#' # neighbors
#' sampling_points <- fibonacci_points(z, a, 1)
#' query_points <- fibonacci_points(z, a, 10)
#' neighbors <- nn_spherical(sampling_points, query_points,
#'                           k = nrow(sampling_points),
#'                           search_radius = 30)
#'
#' thr <- apply_by_direction(r, z, a, m)
#' bin <- binarize_with_thr(r, thr$v1)
#' plot(bin)
#'
#' # Inverse distance weighted interpolation ------------------------------
#'
#' ## See fit_cie_model() for details on the CSV file
#' path <- system.file("external/sky_points.csv",
#'                     package = "rcaiman")
#' sky_points <- read.csv(path)
#' sky_points <- sky_points[c("Y", "X")]
#' colnames(sky_points) <- c("row", "col")
#' head(sky_points)
#'
#' sky_points <- extract_dn(r, sky_points)
#' sky_points <- add_directions(sky_points, z, a)
#'
#' query_points <- raster_grid_in_circle(z, a, 10)
#'
#' neighbors <- nn_spherical(sky_points, query_points,
#'                           k = 10,
#'                           search_radius = 20,
#'                           rule = "any")
#'
#' eval_points <- apply_locally(
#'   sky_points,
#'   query_points,
#'   neighbors,
#'   fun = function(df, chi, ...) {
#'     b_k <- df$Cyan
#'     w_k <- 1 / chi^2
#'     sum(w_k * b_k) / sum(w_k)
#'   },
#'   parallel = FALSE
#' )
#'
#' head(eval_points)
#' summary(eval_points)
#' i <- !is.na(eval_points$output)
#' sky <- interpolate_planar(eval_points[i,],  "output", r,
#'                           k = 3,
#'                           p = 2,
#'                           search_radius = 20,
#'                           rigorous = TRUE)
#' plot(sky)
#' plot(r/sky)
#' }
apply_locally <- function(sampling_points,
                          query_points,
                          neighbors,
                          fun,
                          parallel = FALSE,
                          cores = NULL,
                          logical = TRUE,
                          leave_free = 1
) {

  .check_columns(sampling_points,
                 required_columns = c("zen", "azi"),
                 reserved_columns = NULL)
  .check_columns(query_points,
                 required_columns = c("zen", "azi"),
                 reserved_columns = "output")
  .assert_neighbors(neighbors)
  if (!is.function(fun))
    stop("`fun` must be a function.")
  .check_vector(parallel, "logical", 1)
  .check_vector(cores, "integerish", 1, allow_null = TRUE, sign = "positive")
  .check_vector(logical, "logical", 1)
  .check_vector(leave_free, "integerish", 1, sign = "nonnegative")

  if (parallel) {
    cores <- .cores(cores, logical, leave_free)
    if (cores < 2) parallel <- FALSE
  }

  n <- nrow(query_points)

  .process_one <- function(i) {
    neighbors_of_i <- neighbors[[i]]
    if (length(neighbors_of_i$index) == 0) return(NA)
    fun(sampling_points[neighbors_of_i$index,], neighbors_of_i$distance, i)
  }

  output <- if (parallel) {
    .with_cluster(cores, {
      # Only to avoid note from check, code is OK without this line.
      i <- NA

      foreach::foreach(i = seq_len(n), .combine = c) %dopar% {
        list(.process_one(i))
      }
    })
  } else {
    lapply(seq_len(n), .process_one)
  }

  foo <- unlist(output)
  output <- if (length(foo) == nrow(query_points)) {
    cbind(query_points, output = foo)
  } else {
    output
  }
  output
}
