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
#' {\Phi(0) \cdot f(\theta_\odot)}
#' }
#'
#' where \eqn{\theta} is the zenith angle, \eqn{\phi} is the azimuth angle,
#' and \eqn{\theta_\odot, \phi_\odot} are the zenith and azimuth of the sun disk.
#' The angular distance between zenith and the sun is \eqn{\theta_\odot},
#' which determines the indicatrix value in the normalization denominator.
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
#' Coefficient for an ideal equidistant projection normalized at an off-axis
#' angle of 90 deg.
#'
#' @details
#' In an equidistant projection, also called polar projection, radial distance
#' in image space is linearly proportional to the off-axis angle
#' \insertCite{Schneider2009}{rcaiman}:
#'
#' \deqn{r_{rel} = k \cdot \theta'}
#'
#' where \eqn{r_{rel}} is the relative radius (radius expressed as a fraction of
#' the radius at 90 deg), \eqn{\theta'} is the off-axis angle in radians, and
#' \eqn{k} is a scalar.
#'
#' Enforcing the geometric constraint \eqn{r_{rel} = 1} at
#' \eqn{\theta' = \pi/2} yields \eqn{k = 2 / \pi}.
#'
#' @format Numeric vector of length one.
#'
#' @references \insertAllCited{}
#'
#' @seealso [normalize_projection_coef()] to change the normalization
#'   reference; [orthographic_coef], [equal_area_coef].
#'
#' @export
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equidistant_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' set_internal_geometry(c(80, 60), geometry)
equidistant_coef <- 2/pi

#' Orthographic projection coefficients
#'
#' @description
#' Polynomial coefficients approximating an orthographic projection normalized
#' at an off-axis angle of 90 deg.
#'
#' @details
#' In an orthographic projection, radial distance in image space is proportional
#' to the sine of the off-axis angle \insertCite{Schneider2009}{rcaiman}:
#'
#' \deqn{r_{rel} = \sin(\theta')}
#'
#' where \eqn{r_{rel}} is the relative radius (radius expressed as a fraction of
#' the radius at 90 deg), and \eqn{\theta'} is the off-axis angle in
#' radians.
#'
#' The supplied coefficients fit the radial polynomial described in
#' [set_internal_geometry()] to the orthographic projection over off-axis
#' angles in \eqn{[0, \pi/2]}.
#'
#' @format Numeric vector of length five.
#'
#' @references \insertAllCited{}
#'
#' @seealso [normalize_projection_coef()] to change the normalization
#'   reference; [equidistant_coef], [equal_area_coef].
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = orthographic_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' set_internal_geometry(c(80, 60), geometry)
#' @export
orthographic_coef <- c(0.9998, 0.001840, -0.171600, 0.005626, 0.005848)


#' Equal-area projection coefficients
#'
#' @description
#' Polynomial coefficients approximating an equal-area projection (also
#' called equisolid-angle projection), normalized at an off-axis angle of
#' 90 deg.
#'
#' @details
#' In an equal-area projection, radial distance in image space is
#' proportional to the sine of half the off-axis angle
#' \insertCite{Schneider2009}{rcaiman}:
#'
#' \deqn{r_{rel} = \sqrt{2} \sin(\theta'/2)}
#'
#' where \eqn{r_{rel}} is the relative radius (radius expressed as a fraction of
#' the radius at 90 deg), and \eqn{\theta'} is the off-axis angle in
#' radians.
#'
#' This formulation preserves solid angle, meaning that equal areas in the image
#' correspond to equal portions of the hemisphere.
#'
#' The supplied coefficients fit the radial polynomial described in
#' [set_internal_geometry()] to the equal-area projection over off-axis
#' angles in \eqn{[0, \pi/2]}.
#'
#' @format Numeric vector of length four.
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @seealso [normalize_projection_coef()] to change the normalization
#'   reference; [equidistant_coef], [orthographic_coef].
#'
#' @examples
#' geometry <- list(
#'   model = "radial_projection",
#'   parameters = equal_area_coef,
#'   center_of_vision = c(40, 30),
#'   normalization_radius = 30,
#'   normalization_angle = 90
#' )
#' set_internal_geometry(c(80, 60), geometry)
equal_area_coef <- c(0.7070, 0.000913, -0.0312, 0.00135)
