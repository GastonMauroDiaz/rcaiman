require(rcaiman.registry)
path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
my_entry <- read_registry_snapshot(path)
spec <- get_geometry_spec(my_entry,
                          "exif_01", "simple_method")
rspec <- get_radiometry_spec(spec, "spectral_bands")
caim <- read_caim_raw("data-raw/DSCN4606.NEF",
                      cfa_pattern = rspec$cfa_pattern)
caim <- crop_caim(caim,
                  snapped_zenith_coord = spec$snapped_zenith_coord,
                  horizon_radius = spec$horizon_radius)
caim <- normalize_minmax(caim) * (2^8 - 1)
caim <- round(caim)
caim
write_caim(caim, "inst/external/example", 8)
