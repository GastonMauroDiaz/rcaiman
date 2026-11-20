# HSP project fixture

Subset of the real MA_2021_11_18 project supplied by Gastón Díaz on
2026-11-11, including the subsequent mask edits saved that morning.
Source: `C:/Users/gdiaz/Documents/HSP/Projects/MA_2021_11_18/`.

The `manipulate` directory contains unchanged copies of the PGM image and
configuration files for DSCN6280 (Nikon E5700), about 1.15 MB in total.
RAW photographs, previews, and unrelated project entries are omitted.
This is a fixture for the rcaiman HSP helpers, not a complete HSP project.
`sha256.json` records the SHA-256 digest of every copied file.

- DSCN6280: 756 x 756 pixels, 141 sky marks, one triangle in `_masks.conf`
  and two polygons (four and five vertices) in `_polygons.conf`.
- Decimal commas occur in the sky coefficients and zenith reference.
- Coordinate dots separate column and row; they are not decimal points.

Tests read this fixture and copy it to a temporary directory before testing
writers or replacing both polygon files with `0` to test empty masks. Neither
the fixture nor the supplied original project is edited by the tests.
Nonempty polygons are checked against their actual saved vertices.
The Canon case was removed to keep the fixture small without losing these
test scenarios.
