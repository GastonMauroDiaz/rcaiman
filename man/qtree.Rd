% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qtree.R
\name{qtree}
\alias{qtree}
\title{Perform approximate quad-tree-like segmentation}
\usage{
qtree(r, scale_parameter = 0.2)
}
\arguments{
\item{r}{\linkS4class{SpatRaster}.}

\item{scale_parameter}{Numeric vector of length one. Quad-tree is a top-down
method. This parameter controls the stopping condition. Therefore, it
allows controlling the size of the resulting segments. Ultimately, segments
sizes will depend on both this parameter and the heterogeneity of \code{r}.}
}
\value{
A single layer image of the class \linkS4class{SpatRaster} with integer
values.
}
\description{
This function performs an efficient hierarchical segmentation of the planar
space inspired by the quad-tree algorithm. Instead of applying recursive
subdivision cell-by-cell, it uses several predefined segmentation levels and
evaluates local heterogeneity to decide whether finer subdivisions are
justified and should be retained.
}
\details{
Segments at each level are organized in
such a way that each coarser cell could theoretically be subdivided into four
finer subcells. However, the process does not follow a strict top-down
recursive logic as in a canonical quad-tree.

The function computes a metric (\emph{delta}) at a maximum of five levels. \emph{delta}
is defined as the sum of the standard deviation of its subregions minus the
standard deviation of the parent region. If the \emph{delta} is larger than a
user-defined \code{scale_parameter}, then the finer segmentation level is retained
locally. This evaluation is applied globally at each level, not recursively
from cell to cell.

This implementation results in a segmentation that \emph{resembles} a quad-tree in
appearance, but does not guarantee structural consistency across levels. In
particular, small segments may appear nested within regions that were not
formally subdivided from above, and not all parents have exactly four
children. The benefit of this approach is a significant reduction in
computational cost, at the expense of formal consistency with the classic
quad-tree hierarchy.
}
\examples{
\dontrun{
caim <- read_caim() \%>\% normalize_minmax()
seg <- qtree(caim, scale_parameter = 0.05)
plot(caim$Blue)
plot(extract_feature(caim$Blue, seg))
plot(extract_feature(seg, seg, length))
}
}
