# Testing mission: prevent incidental algorithmic changes within the intended
# use of the package. Prefer package builders when creating test fixtures.
# See tests/README.md for focused runs and criteria for adding tests.

library(testthat)
library(rcaiman)

test_check("rcaiman")
