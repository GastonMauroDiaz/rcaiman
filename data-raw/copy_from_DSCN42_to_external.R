from <- "data-raw/DSCN4248/sky_points/sky_points.csv"
to <- "inst/external/sky_points.csv"
file.copy(from, to, overwrite = TRUE)

from <- "data-raw/DSCN4248/sun_disk_mark/sun_disk_mark.csv"
to <- "inst/external/sun_disk_mark.csv"
file.copy(from, to, overwrite = TRUE)

from <- "data-raw/DSCN4248/img_l0/DSCN4248.tif"
to <- "inst/external/DSCN4248.tif"
file.copy(from, to, overwrite = TRUE)

from <- "data-raw/DSCN4248/sensor_data/DSCN4248.tsv"
to <- "inst/external/DSCN4248.tsv"
