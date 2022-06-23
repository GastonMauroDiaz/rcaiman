# return or change file extensions
# Author: Robert J. Hijmans
# Date : October 2008
# Version 1.0
# Licence GPL v3
# https://github.com/rspatial/raster/blob/6f90d642efb3573de3df77b67dd7799248151a2a/R/extension.R

extension <- function(filename, value=NULL, maxchar=10) {
  if (!is.null(value)) {
    extension(filename) <- value
    return(filename)
  }
  lfn <- nchar(filename)
  ext <- list()
  for (f in 1:length(filename)) {
    extstart <- -1
    for (i in lfn[f] : 2) {
      if (substr(filename[f], i, i) == ".") {
        extstart <- i
        break
      }
    }
    if (extstart > 0) {
      ext[f] <- substr(filename[f], extstart, lfn[f])
    } else {
      ext[f] <- ""
    }
  }
  ext <- unlist(ext)
  ext[nchar(ext) > maxchar] <- ''
  return(ext)
}


'extension<-' <- function(filename, value) {
  # value <- trim(value)
  value <- trimws(value)
  if (value != "" & substr(value, 1, 1) != ".") {
    value <- paste(".", value, sep="")
  }
  lfn <- nchar(filename)
  fname <- list()
  for (f in 1:length(filename)) {
    extstart <- -1
    for (i in lfn[f] : 2) {
      if (substr(filename[f], i, i) == ".") {
        extstart <- i
        break
      }
    }
    if (extstart > 0 & (lfn[f] - extstart) < 8) {
      fname[f] <- paste(substr(filename[f], 1, extstart-1), value, sep="")
    } else {
      fname[f] <- paste(filename[f], value, sep="")
    }
  }
  return( unlist(fname) )
}
