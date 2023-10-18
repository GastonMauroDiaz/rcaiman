#' Model-based local thresholding
#'
#' Transform background digital number into threshold values
#'
#' This function transforms background digital numbers into threshold values by
#' means of the Equation 1 from \insertCite{Diaz2018;textual}{rcaiman}, which is
#' a linear function with the slope modified by a weighting parameter. This
#' simple function was found by studying canopy models, also known as targets,
#' which are perforated surfaces made of a rigid and dark material . These
#' models were backlighted with homogeneous lighting, photographed with a Nikon
#' Coolpix 5700 set to acquire in JPEG format, and those images were gamma back
#' corrected with a default gamma value equal to 2.2 (see [gbc()]). Results
#' shown that the optimal threshold value was linearly related with the
#' background digital number (see Figure 1 and Figure 7 from
#' \insertCite{Diaz2018;textual}{rcaiman}). This shifted the aim from finding
#' the optimal threshold, following \insertCite{Song2014}{rcaiman} method, to
#' obtaining the background DN as if the canopy was not there, as
#' \insertCite{Lang2010}{rcaiman} proposed.
#'
#' ## Working principle
#'
#' \insertCite{Diaz2018;textual}{rcaiman} observed the following linear
#' relationship between the background value, usually the sky digital number
#' (SDN), and the optimal threshold value (OTV):
#'
#' | \eqn{IV = a + b \cdot SDN}||||||||||||(Equation 1a)|
#' |:---|-|-|-|-|-|-|-|-|-|-|-|-----------:|
#'
#' | \eqn{OTV = a + b \cdot w \cdot SDN}||||||||||||(Equation 1b)|
#' |:---|-|-|-|-|-|-|-|-|-|-|-|-----------:|
#'
#' were IV is the initial value \insertCite{Wagner2001}{rcaiman}, which is the
#' boundary between SDN and the mixed pixels, i.e, the pixels that are neither
#' Gap or Non-gap \insertCite{Macfarlane2011}{rcaiman}, \eqn{a} and \eqn{b} are
#' the intercept and slope coefficients, respectively, and \eqn{w} is a
#' weighting parameter that takes into account that OTV is always lower than IV.
#' If SDN is calculated at the pixel level, a local thresholding method can be
#' applied by evaluating, pixel by pixel, if the below canopy digital number
#' (CDN) is greater than the OTV. Formally, If \eqn{CDN>OTV}, then assign _Gap_
#' class, else assign _Non-gap_ class.
#'
#' This conclusion drawn from an image processing point of view matches with
#' previous findings drawn from a radiometric measurement paradigm, which are
#' introduced next.
#'
#' \insertCite{Cescatti2007;textual}{rcaiman}) posed that cameras can be used as
#' a radiation measurement device if they are properly calibrated. This method,
#' denominated by the author as LinearRatio, seeks to obtain the transmittance
#' (T) as the ratio of below to above canopy radiation:
#'
#' | \eqn{T = CDN/SDN}||||||||||||(Equation 2)|
#' |:---|-|-|-|-|-|-|-|-|-|-|-|-----------:|
#'
#' were CDN is below canopy digital number (DN), i.e., the DN extracted from a
#' canopy hemispherical photograph.
#'
#' The LinearRatio method uses T as a proxy for gap fraction. It ideally
#' requires twin cameras, one below and the other above the canopy. In contrast,
#' \insertCite{Lang2010;textual}{rcaiman}) proposed to obtain SDN by manually
#' selecting pure sky pixels from canopy hemispherical photographs and
#' reconstructing the whole sky by subsequent modeling and interpolating---this
#' method is often referred to as LinearRatio single camera or LinearRatioSC.
#'
#' Equation 2 can be seen as a standardization of the distance between CDN and
#' SDN. With that in mind, it is useful to rewrite Equation 1b as an inequality
#' that can be evaluated to return a logical statement that is directly
#' translated into the desired binary classification:
#'
#' | \eqn{CDN > a + b \cdot w \cdot SDN}||||||||||||(Equation 3)|
#' |:---|-|-|-|-|-|-|-|-|-|-|-|-----------:|
#'
#' Then, combining Equation 2 and 3, we find that
#' \insertCite{Diaz2018;textual}{rcaiman} parameters can be applied to T:
#'
#' | \eqn{CDN/SDN > a + b \cdot w \cdot SDN/SDN}||||||||||||(Equation 4a)|
#' |:---|-|-|-|-|-|-|-|-|-|-|-|-----------:|
#'
#' | \eqn{T > a + b \cdot w}||||||||||||(Equation 4b)|
#' |:---|-|-|-|-|-|-|-|-|-|-|-|-----------:|
#'
#' From Equation 2 it is evident that any bias introduced by the camera optical
#' and electronic system will be canceled during the calculation of T as long as
#' only one camera is involved. Therefore, After examining Equation 4b, we can
#' conclude that intercept 0 and slope 1 are theoretically correct.In addition,
#' the w parameter can be used to filter out mixed pixels and sample the most
#' pure sky pixels possible. The greater w, the greater the possibility of
#' selecting pure sky pixels.
#'
#' @note
#'
#' It is worth noting that Equation 1 was developed with 8-bit images, so
#' calibration of new coefficient should be done in the 0 to 255 domain since
#' that is what [thr_mblt()] expect, although the input `dn` should be
#' normalized. The latter was a design decision aiming to harmonize the whole
#' package, although it might sound counter intuitive.
#'
#' Nevertheless, new empirical calibrations may be unnecessary since the values
#' -7.8 `intercept` and 0.95 `slope` observed with back-gamma corrected JPEG
#' files produced with the Nikon Coolpix 5700 camera are sufficiently close to
#' the theoretical values to suggest that using them is a correct assumption. In
#' addition, instead of calibrating JPEG images, users should adopt raw file
#' acquisition ([read_caim_raw()]).
#'
#' To apply the weighting parameter (w) from Equation 1, just provide the
#' argument `slope` as \eqn{slope \times w}.
#'
#' @param dn Numeric vector or [SpatRaster-class]. Digital number of the
#'   background. These values should be normalized and, if they are extracted
#'   from a JPEG image, gamma back corrected.
#' @param intercept,slope Numeric vector of length one. These are linear
#'   function coefficients.
#'
#' @export
#'
#' @return An object of the same class and dimensions than `dn`.
#'
#' @family Binarization Functions
#' @seealso [normalize()], [gbc()], [apply_thr()] and [regional_thresholding()].
#'
#' @references \insertAllCited{}
#'
#' @examples
#' thr_mblt(gbc(125), -8, 1)
thr_mblt <- function (dn, intercept, slope) {
  stopifnot(length(intercept) == 1)
  stopifnot(class(intercept) == "numeric")
  stopifnot(length(slope) == 1)
  stopifnot(class(slope) == "numeric")

  dn <- dn * 255

  if (.get_max(dn) > 255) warning("\"dn\" values should be normalized")
  dn[dn > 255] <- 255

  thr <- intercept  + slope * dn
  thr[thr < 0] <- 0
  thr[is.na(thr)] <- 0

  thr / 255
}
