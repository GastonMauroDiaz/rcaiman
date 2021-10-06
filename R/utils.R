radian2degree <- function(x) x * 180 / pi
degree2radian <- function(x) x * pi / 180

.get_max <- function(r) max(r[], na.rm = TRUE)
.get_min <- function(r) min(r[], na.rm = TRUE)

.is_even <- function(x) round(x/2) == x/2

.is_whole <- function(x) round(x) == x

.is_integerish <- function(x) x == round(x)

.check_if_r_was_normalized <- function(r) {
  if (max(r[], na.rm = TRUE) > 1)
    warning("Please check if \"r\" was correctly normalized")
}

.this_requires_EBImage <- function() {
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop(paste("Package \"EBImage\" needed for this function to work.",
               "Please install it. Find instructions here:",
               "https://bioconductor.org/packages/release/bioc/html/EBImage.html"
               ),
         call. = FALSE)
  }
}
