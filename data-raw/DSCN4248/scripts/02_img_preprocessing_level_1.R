# Split bands to allow image marking

require(rcaiman.registry.data)
rm(list = ls())
entry_Nikon_Coolpix5700.FCE9.CIEFAP
summary(entry_Nikon_Coolpix5700.FCE9.CIEFAP)

my_band_spec <- get_band_spec(entry_Nikon_Coolpix5700.FCE9.CIEFAP, "raw", "cygm")

sensor_data <- read_caim("img_l0/DSCN4248.tif")
caim <- split_cfa(sensor_data, cfa_pattern = my_band_spec$cfa_pattern)
write_caim(caim$Cyan, path = "img_l1/DSCN4248.tif", 8)

export_registry_snapshot(entry_Nikon_Coolpix5700.FCE9.CIEFAP,
                         "raw", "cygm", file = "img_l1/registry_entry.yaml")
