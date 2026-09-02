# Testing mission: prevent incidental algorithmic changes within the intended
# use of the package. Prefer package builders when creating test fixtures.

library(testthat)
library(rcaiman)

test_check("rcaiman")
