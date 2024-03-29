% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combineRGBMatrices.R
\name{combineRGBMatrices}
\alias{combineRGBMatrices}
\title{combineRGBMatrices}
\usage{
combineRGBMatrices(
  c1.mat,
  c2.mat,
  combine_method = "scaled_geomean",
  max_per_channel = 255
)
}
\arguments{
\item{c1.mat}{A numeric matrix with three columns corresponding to red, green, blue, with numerical
values ranging from 0 to 255 to be combined numerically with \code{c2.mat}.}

\item{c2.mat}{A numeric matrix with three columns corresponding to red, green, blue, with numerical
values ranging from 0 to 255 to be combined numerically with \code{c1.mat}.}

\item{combine_method}{Method of combining colors, can be \code{'scaled_geomean'},
\code{'standard'}/\code{'euclidean'}, \code{'negative_euclidean'}, \code{'mean'}, \code{'scaled_geomean'}
(default = \code{'scaled_geomean'})}

\item{max_per_channel}{Maximal color value per RGB channel (default 255).}
}
\value{
A 3-column RGB matrix of combined colors.
}
\description{
Given 2 different matrices of colors, combine the colors numerically. This is used in \code{makeTwoColorEncodeFunction()}.
}
\examples{

c1.mat <- matrix( c(255, 100, 0 ), ncol = 3 )
c2.mat <- matrix( c( 0, 50, 255 ), ncol = 3 )
c12.mat <- combineRGBMatrices( c1.mat, c2.mat, "euclidean" )

}
\seealso{
\code{\link[=makeTwoColorEncodeFunction]{makeTwoColorEncodeFunction()}}
}
