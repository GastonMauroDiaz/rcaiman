path <- system.file("external/DSCN4500.JPG", package = "rcaiman")
caim <- read_caim(path)

set.seed(1)
sampling_indices <- data.frame(
  col = runif(100, 1, ncol(caim)) %>% trunc(),
  row = runif(100, 1, nrow(caim)) %>% trunc()
)
m <- is.na(caim$Red)
names(m) <- "mask 10"
ratio <- caim$Blue/caim$Red
names(ratio) <- "ratio with sky_cie"

pts <- extract_dn(c(caim,  m, ratio), sampling_indices, use_window = FALSE)
pts <- cbind(coord_from_indices(pts$col, pts$row), pts)
write_points(pts, "data-raw/test.csv")

foo <- read_points("data-raw/test.csv")
head(foo)
