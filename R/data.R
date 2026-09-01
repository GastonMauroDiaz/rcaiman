#' Set of 15 CIE Standard Skies
#'
#' @description
#' The Commission Internationale de l’Éclairage (CIE; International Commission
#' on Illumination) standard \insertCite{CIE2004}{rcaiman} defines 15
#' pre-calibrated sky luminance distributions, each described by a pair of
#' analytical functions:
#'
#' \describe{
#'  \item{*Gradation function*}{
#'  \deqn{\Phi(\theta) = 1 + a \cdot \exp\left(\frac{b}{\cos\theta}\right)}
#'  }
#'  \item{*Indicatrix function*}{
#'  \deqn{f(\chi) = 1 + c \cdot \left[ \exp(d \cdot \chi) - \exp\left(d \cdot \frac{\pi}{2}\right) \right] + e \cdot \cos^2\chi}
#'  }
#' }
#'
#' Combined, they can predict the relative radiance \eqn{\rho_R} in any sky
#' direction \eqn{(\theta, \phi)} as:
#'
#' \deqn{
#' \hat{\rho_R}^{\circ}(\theta, \phi) =
#' \frac{\Phi(\theta) \cdot f(\chi(\theta, \phi; \theta_\odot, \phi_\odot))}
#' {\Phi(0) \cdot f(\chi(\theta, \phi; 0, 0))}
#' }
#'
#' where \eqn{\theta} is the zenith angle, \eqn{\phi} is the azimuth angle,
#' and \eqn{\theta_\odot, \phi_\odot} are the zenith and azimuth of the sun disk.
#'
#' @format `data.frame` with 15 rows and 8 columns:
#' \describe{
#'   \item{`a`}{gradation function parameter.}
#'   \item{`b`}{gradation function parameter.}
#'   \item{`c`}{indicatrix function parameter.}
#'   \item{`d`}{indicatrix function parameter.}
#'   \item{`e`}{indicatrix function parameter.}
#'   \item{`indicatrix_group`}{factor with six categories and numerical tags.}
#'   \item{`general_sky_type`}{factor with three categories: `"Overcast"`,
#'     `"Clear"`, and `"Partly cloudy"`.}
#'   \item{`description`}{user-friendly description of the sky type.}
#' }
#'
#' @source \insertCite{Li2016;textual}{rcaiman}
#'
#' @references \insertAllCited{}
#'
"cie_table"

#' Equidistant projection coefficient
#'
#' @description
#' Coefficient for an ideal equidistant projection normalized at the horizon.
#'
#' @details
#' In an equidistant projection, also called polar projection, radial distance
#' in image space is linearly proportional to the zenith angle
#' \insertCite{Schneider2009}{rcaiman}:
#'
#' \deqn{r_{rel} = k \cdot \theta}
#'
#' where \eqn{r_{rel}} is the relative radius (radius expresed as a fraction of
#' the radius at 90 deg), and \eqn{\theta} is the zenith angle in radians and
#' \eqn{k} is a scalar.
#'
#' Enforcing the geometric constraint \eqn{r_{rel} = 1} at
#' \eqn{\theta = \pi/2} yields \eqn{k = 2 / \pi}.
#'
#' @format Numeric vector of length one.
#'
#' @references \insertAllCited{}
#'
#' @seealso [orthographic_coef], [equal_area_coef]
#'
#' @export
#'
#' @examples
#' zenith_image(diameter, equidistant_coef)
equidistant_coef <- 2/pi

#' Orthographic projection coefficients
#'
#' @description
#' Polynomial coefficients approximating an orthographic projection normalized
#' at the horizon.
#'
#' @details
#' In an orthographic projection, radial distance in image space is proportional
#' to the sine of the zenith angle \insertCite{Schneider2009}{rcaiman}:
#'
#' \deqn{r_{rel} = \sin(\theta)}
#'
#' where \eqn{r_{rel}} is the relative radius (radius expresed as a fraction of
#' the radius at 90 deg), and \eqn{\theta} is the zenith angle in
#' radians.
#'
#' In `rcaiman`, radial projection models are defined as polynomials (see
#' [calc_relative_radius()]):
#'
#' \deqn{r_{rel}(\theta) = \sum_{p=1}^{n} a_p \theta^p}
#'
#' where the coeficients are \eqn{a_1, a_2, ..., a_n)} and the coeficient at
#' position `p` multiplies \eqn{\theta^p}.
#'
#' The coefficients provided here were obtained by fitting this
#' polynomial to the orthographic projection over the interval \eqn{[0, \pi/2]},
#' and satisfy the normalization constraints enforced by
#' [validate_projection_function()].
#'
#' @format Numeric vector of length four.
#'
#' @references \insertAllCited{}
#'
#' @seealso [equidistant_coef], [equal_area_coef]
#'
#' @examples
#' z <- zenith_image(diameter, orthographic_coef)
#' @export
orthographic_coef <- c(0.9972, 0.01719, -0.2008, 0.028018)


#' Equal-area projection coefficients
#'
#' @description
#' Polynomial coefficients approximating an equal-area projection (also
#' called equisolid-angle projection), normalized at the horizon.
#'
#' @details
#' In an equal-area projection, radial distance in image space is
#' proportional to the sine of half the zenith angle
#' \insertCite{Schneider2009}{rcaiman}:
#'
#' \deqn{r_{rel} = \sqrt{2} \sin(\theta/2)}
#'
#' where \eqn{r_{rel}} is the relative radius (radius expresed as a fraction of
#' the radius at 90 deg), and \eqn{\theta} is the zenith angle in
#' radians.
#'
#' This formulation preserves solid angle, meaning that equal areas in the image
#' correspond to equal portions of the hemisphere.
#'
#' In `rcaiman`, radial projection models are defined as polynomials (see
#' [calc_relative_radius()]):
#'
#' \deqn{r_{rel}(\theta) = \sum_{p=1}^{n} a_p \theta^p}
#'
#' where the coeficients are \eqn{a_1, a_2, ..., a_n)} and the coeficient at
#' position `p` multiplies \eqn{\theta^p}.
#'
#' The coefficients provided here were obtained by fitting this
#' polynomial to the equal-area projection over the interval \eqn{[0, \pi/2]},
#' and satisfy the normalization constraints enforced by
#' [validate_projection_function()].
#'
#' @format Numeric vector of length four.
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @seealso [equidistant_coef], [orthographic_coef]
#'
#' @examples
#' z <- zenith_image(diameter, equal_area_coef)
equal_area_coef <- c(0.7070, 0.000913, -0.0312, 0.00135)
