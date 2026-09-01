.this_requires_EBImage <- function() {
  if (!requireNamespace("EBImage", quietly = TRUE)) {
    stop(paste("Package \"EBImage\" needed for this function to work.",
               "Please install it. Find instructions here:",
               "https://bioconductor.org/packages/release/bioc/html/EBImage.html"
    ),
    call. = FALSE)
  }
}

.this_requires_multimode <- function() {
  if (!requireNamespace("multimode", quietly = TRUE)) {
    stop(paste("Package \"multimode\" needed for this function to work.",
               "Please install it."),
         call. = FALSE)
  }
}
