#' Segment a hemisphere into quasi-equal-area cells
#'
#' @description
#' Segment a hemispherical view into quasi-equal-area cells by assigning each
#' pixel a grid-cell identifier.
#'
#' @details
#' Quasi-equal-area partitions can be constructed in two ways:
#' \describe{
#'   \item{Equal angle rings}{The zenith width (deg) of the rings is constant.
#'   Equal-area behavior is achieved by varying the number of azimuthal sectors
#'   per ring.}
#'   \item{Equal area rings}{Ring boundaries are computed in `cos(z)`,
#'   producing rings of similar projected area, later subdivided into a constant
#'   number of azimuthal sectors.}
#' }
#'
#' @param n_cells numeric vector of length one. Target number of cells for the
#'   segmentation. The final number of generated cells may differ slightly,
#'   since equal-area constraints take precedence over an exact count.
#' @param ring_mode character vector of length one. Ring construction method:
#'   `"equal_area_rings"` or `"equal_angle_rings"`. See *Details*.
#'
#' @inheritParams sky_grid_segmentation
#'
#' @return Numeric [terra::SpatRaster-class] with one layer containing integer
#'   cell identifiers. The object carries attribute `ring_mode`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' z <- zenith_image(500, lens())
#' a <- azimuth_image(500, lens())
#' seg <- sky_equalarea_segmentation(z, a, n_cells = 512)
#' plot(seg)
#' }
sky_equalarea_segmentation <- function(z, a, n_cells,
                                 ring_mode = "equal_angle_rings") {
  .check_r_z_a_m(NULL, z, a)
  .check_vector(n_cells, "numeric", 1, sign = "positive")
  .assert_choice(ring_mode, c("equal_angle_rings", "equal_area_rings"))


  # Determinar número de anillos
  R <- round(sqrt(n_cells))
  if (R < 1) R <- 1
  if (R > n_cells) R <- n_cells

  # ------------------------------------------------------------
  # 1) Calcular Ni (número de sectores por anillo)
  # ------------------------------------------------------------

  if (ring_mode == "equal_area_rings") {
    Ni <- rep(floor(n_cells / R), R)
    # rem <- n_cells - sum(Ni)
    # if (rem > 0) Ni[1:rem] <- Ni[1:rem] + 1

    # Cortes en cos(z)
    cos_breaks <- numeric(R + 1)
    cos_breaks[1] <- 1
    for (i in 1:R) cos_breaks[i + 1] <- cos_breaks[i] - Ni[i] / n_cells

  } else if (ring_mode == "equal_angle_rings") {
    delta_z <- 90 / R
    z_breaks <- seq(0, 90, by = delta_z)

    # Área relativa del anillo (proporcional a cos(z1) - cos(z2))
    ring_area <- numeric(R)
    for (i in 1:R) {
      z1 <- z_breaks[i]   * pi/180
      z2 <- z_breaks[i+1] * pi/180
      ring_area[i] <- abs(cos(z1) - cos(z2))
    }

    # Normalizar áreas a celdas
    Ni <- round(n_cells * ring_area / sum(ring_area))

    # diff_total <- n_cells - sum(Ni)
    # if (diff_total != 0) {
    #   adj <- sample(1:R, abs(diff_total), replace = TRUE)
    #   for (j in adj) Ni[j] <- Ni[j] + sign(diff_total)
    # }

    # Convertir a límites en cos(z)
    cos_breaks <- cos(z_breaks * pi/180)
  }

  # offsets para IDs
  offsets <- c(0, cumsum(Ni))

  res <- terra::rast(z)
  vz  <- terra::values(z)
  va  <- terra::values(a)
  cos_z <- cos(vz * pi/180)

  out_vals <- rep(NA, length(vz))

  for (i in 1:R) {

    upper <- cos_breaks[i]
    lower <- cos_breaks[i+1]

    if (i < R) {
      idx <- which(!is.na(vz) & cos_z <= upper & cos_z > lower)
    } else {
      idx <- which(!is.na(vz) & cos_z <= upper & cos_z >= lower)
    }

    if (length(idx) > 0) {
      phi_span <- 360 / Ni[i]
      phi <- va[idx] %% 360
      j <- floor(phi / phi_span) + 1
      out_vals[idx] <- offsets[i] + j
    }
  }

  terra::values(res) <- out_vals

  names(res) <- "Sky segments of quasi-equal area"
  attr(res, "ring_mode") <- ring_mode
  res
}
