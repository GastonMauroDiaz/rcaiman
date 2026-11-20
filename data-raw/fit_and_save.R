require(rcaiman.registry)
path <- system.file("external/DSCN4248.tif", package = "rcaiman")
metadata <- read.delim(gsub("tif", "tsv", path))
sensor_data <- read_caim(path)

# Get specs
path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
my_entry <- read_registry_snapshot(path)
raw <- get_signature(my_entry, "raw")
cygm <- get_band_spec(raw, "cygm")
raw_geometry <- get_geometry_spec(cygm, "simple_method")
raw_radiometry <- get_radiometry_spec(raw_geometry, "photometric_sphere")

# Represent image geometry
g <- set_internal_geometry(raw$dim, raw_geometry)
d <- direction_image(g)
m <- !is.na(g$off_axis)
w <- area_weight_image(g, raw_geometry)

# Read and preprocess the image
caim <- split_cfa(sensor_data, cfa_pattern = cygm$cfa_pattern)
caim
r <- caim$Cyan
ff <- trend_surface_image(r, raw_radiometry$parameters[[metadata$aperture]]$Cyan, m = m)
r <- apply_flat_field(r, ff)

path <- system.file("external/sky_points.csv",
                    package = "rcaiman")
sky_points <- read_points(path) %>% add_indices()

sky_points <- extract_dn(r, sky_points, use_window = TRUE)
sky_points <- add_directions(sky_points, d = d)
head(sky_points)
rr <- calc_relative_radiance(sky_points, "Cyan")

path <- system.file("external/sun_disk_mark.csv",
                    package = "rcaiman")
sun_indices <- read_points(path) %>% add_indices()
sun_angles <- add_directions(sun_indices, d = d)

set.seed(7)
model <- fit_cie_model(
  rr,
  sun_angles,
  general_sky_type = c("Overcast" , "Partly cloudy"),
  parallel = TRUE
)

saveRDS(model, "inst/external/cie_model.RDS")
