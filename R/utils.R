.radian2degree <- function(x) x * 180 / pi
.degree2radian <- function(x) x * pi / 180

.get_max <- function(r) max(r[], na.rm = TRUE)
.get_min <- function(r) min(r[], na.rm = TRUE)

.decode_label <- function(label) {
  sector_ID <- trunc(label / 1000)
  rings_ID <- label - sector_ID * 1000
  data.frame(sector_ID, rings_ID)
}

.calc_rmse <- function(x) sqrt(mean(x^2))

.calc_angular_distance <- function(z1, a1, z2, a2) {
  #https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  acos(pmax(pmin(cos(z1) * cos(z2) + sin(z1) * sin(z2) * cos(abs(a2 - a1)), 1), -1))
}

.extension <- function(file_name, new_extension = "tif") {
  file_name <- filenamer::as.filename(file_name)
  file_name <- filenamer::insert(file_name, ext = new_extension, replace = TRUE)
  as.character(file_name)
}
