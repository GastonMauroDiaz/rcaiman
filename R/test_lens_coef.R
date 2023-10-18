#' Test lens projection functions
#'
#' Test if a lens projection function will work between the 0-to-1 range.
#'
#' The package tolerate a number very close to 1 but not exactly 1 as long as it
#' is greater than 1. Therefore, when the test fails at this *"Test that
#' lens projection function does not predict values barely below one"*, the best
#' practice is to manually edit the last coefficient. For instance, changing it
#' from -0.0296 to -0.0295. See [testthat::expect_equal()] for further details.
#'
#' If it fails in *"Test that lens projection function works between the
#' 0-to-1 range"*, collecting data again might be necessary.
#'
#' @inheritParams zenith_image
#'
#' @family Lens Functions
#'
#' @export
#'
#' @return Returns `invisible(TRUE)` and print "Test passed" if all tests
#'   pass, otherwise throws an error.
#'
#' @examples
#' test_lens_coef(lens("Nikon_FCE9"))
#' test_lens_coef(2/pi)
test_lens_coef <- function(lens_coef) {
  testthat::test_that(
    "Test that lens projection function works between the 0-to-1 range",
    {
      testthat::expect_equal(calc_relative_radius(0, lens_coef) %>%
                               round(., 2), 0)
      testthat::expect_equal(calc_relative_radius(90, lens_coef) %>%
                               round(., 2), 1)
    }
  )
  testthat::test_that(
    "Test that lens projection function does not predict values barely below one",
    {
      testthat::expect_true(calc_relative_radius(90, lens_coef) >= 1)
    }
  )
}
