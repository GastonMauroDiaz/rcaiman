#' Test lens projection functions
#'
#' Test that lens projection function works between the 0-to-1 range.
#'
#' @inheritParams zenith_image
#'
#' @family Lens functions
#'
#' @export
#'
#' @examples
#' test_lens_coef(lens("Nikon_FCE9"))
#' test_lens_coef(2 / pi)
#' test_lens_coef(c(1.06065, -0.49054, 0.14044))
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
}
