require(rcaiman.registry)
path <- system.file("external/Nikon_Coolpix5700.FCE9_raw.yaml", package = "rcaiman")
my_entry <- read_registry_snapshot(path)
spec <- get_geometry_spec(my_entry,
                          "exif_01", "simple_method")
rspec <- get_radiometry_spec(spec, "spectral_bands")
caim <- read_caim_raw("data-raw/DSCN4606.NEF",
                      cfa_pattern = rspec$cfa_pattern)
side <- spec$horizon_radius*2
# TODO (2026-09-02): Temporary transition. Migrate the legacy
# `zenith_col_row` contract before treating this derivation as canonical.
upper_left_indices <- spec$zenith_col_row - spec$horizon_radius
caim <- crop_caim(caim,
                  upper_left_indices = upper_left_indices,
                  width = side,
                  height = side)
caim <- normalize_minmax(caim) * (2^8 - 1)
caim <- round(caim)
caim
write_caim(caim, "inst/external/example", 8)
