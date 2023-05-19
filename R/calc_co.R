.calc_co <- function(bin, m, z, a, angle_width) {
  g <- sky_grid_segmentation(z, a, angle_width)
  g[!m] <- NA
  ds <- extract_feature(bin, g, return_raster = FALSE)

  ids <- .decode_label(as.numeric(names(ds)))
  rcl <- tapply(ids$sector_ID, ids$ring_ID, length)
  rcl <- data.frame(ring_id = names(rcl), no = rcl)
  .n <- ids$ring_ID
  for (i in seq_along(rcl$ring_id)) {
    .n[.n == rcl$ring_id[i]] <- rcl$no[i]
  }
  mx_angles <- ids$ring_ID * angle_width * pi/180
  mn_angles <- mx_angles - angle_width * pi/180
  sum(ds * ((cos(mn_angles) - cos(mx_angles))/.n))
}
