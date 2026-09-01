require(rcaiman.registry.data)
summary(entry_Nikon_Coolpix5700.FCE9.CIEFAP)
export_registry_snapshot(entry_Nikon_Coolpix5700.FCE9.CIEFAP,
                         "exif_jpeg", "delta_t", NULL, "wagner2001_method",
                         "inst/external/Nikon_Coolpix5700.FCE9_jpeg.yaml")
export_registry_snapshot(entry_Nikon_Coolpix5700.FCE9.CIEFAP,
                         "exif_raw", "simple_method", "spectral_bands", "simple_method",
                         "inst/external/Nikon_Coolpix5700.FCE9_raw_poly.yaml")
export_registry_snapshot(entry_Nikon_Coolpix5700.FCE9.CIEFAP,
                         "exif_01", "simple_method", "spectral_bands", "flat_field_reconstruction",
                         "inst/external/Nikon_Coolpix5700.FCE9_raw2.yaml")
