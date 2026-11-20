# To generate the sensor data and marks for inst/external

- **Project folder:** `DSCN4248`
- **Description:** To generate the sensor data and marks for inst/external
- **Created:** 2026-10-28
- **Report generated:** 2026-10-28T13:05:57

## Overview

- RAW files: 1
- Sensor images: 1
- Registered collections: 4
- Scripts: 2

## Registry Validation

| Collection | Registry | Status | Validated | Exact | Accepted | Rejected | Problems |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| - | - | No eligible collections | - | 0 | 0 | 0 | 0 |

## Project Structure

- `.Rproj.user/` - 21 file(s)
- `.rcaiman/` - 1 file(s) - reserved
- `img_l0/` - 2 file(s) - preprocessed_image_collection
- `img_l1/` - 3 file(s) - preprocessed_image_collection
- `raw_files/` - 1 file(s) - reserved
- `registry_drafts/` - 0 file(s) - reserved
- `registry_reviews/` - 0 file(s) - reserved
- `registry_validations/` - 0 file(s) - reserved
- `scripts/` - 2 file(s) - reserved
- `sensor_data/` - 3 file(s) - reserved
- `sky_points/` - 3 file(s) - point_collection
- `sun_disk_mark/` - 3 file(s) - point_collection

## Collections

| Folder | Class | Relevant files | Previous level | Image frame | Built from |
| --- | --- | ---: | --- | --- | --- |
| `img_l0` | preprocessed_image_collection | 1 | - | - | - |
| `img_l1` | preprocessed_image_collection | 1 | img_l0 | - | - |
| `sky_points` | point_collection | 1 | - | img_l1 | - |
| `sun_disk_mark` | point_collection | 1 | - | img_l1 | - |

## Processing Story

### `scripts/01_img_preprocessing_level_0.R`

Degrade sensor data to 8-bits to save space in packaging

| Collection | Class | Previous level | Image frame | Built from | Registry |
| --- | --- | --- | --- | --- | --- |
| `img_l0` | preprocessed_image_collection | - | - | - | - |

### `scripts/02_img_preprocessing_level_1.R`

Split bands to allow image marking

| Collection | Class | Previous level | Image frame | Built from | Registry |
| --- | --- | --- | --- | --- | --- |
| `img_l1` | preprocessed_image_collection | img_l0 | - | - | `img_l1/registry_entry.yaml` |

