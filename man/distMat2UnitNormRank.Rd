% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distMat2UnitNormRank.R
\name{distMat2UnitNormRank}
\alias{distMat2UnitNormRank}
\alias{negDistMat2UnitNormRank}
\title{distMat2UnitNormRank negDistMat2UnitNormRank}
\usage{
distMat2UnitNormRank(mat, lower_is_closer = TRUE)

negDistMat2UnitNormRank(mat)
}
\arguments{
\item{mat}{A numerical matrix containing distances.}

\item{lower_is_closer}{Logical indicating that lower is closer.}
}
\value{
A symmetric matrix of ranks. The matrix has the attribute 'lower_is_closer' set to
'TRUE'.
}
\description{
Convert a symmetrical numerical matrix of distances to a matrix of scaled ranks (from 0 to 1).
Note: Only the lower side of the matrix is used. Data on the upper right are assumed to be redundant. These
functions are intended to convert a matrix of distance or similarity values into a proper form for applying
hierarchical clustering with the \code{gsnPareNetGenericHierarchic()} function.
}
\details{
The difference between \code{distMat2UnitNormRank()} and \code{negDistMat2UnitNormRank()} is that
\code{negDistMat2UnitNormRank()} takes only the \code{mat} argument, and negates it prior to calling
\code{distMat2UnitNormRank()}.
}
\section{Functions}{
\itemize{
\item \code{negDistMat2UnitNormRank()}: Takes the same parameter distMat2UnitNormRank, but multiplies the distance by -1 first.

}}
\examples{

# For log Fisher values, lower is closer and more significant.
mat.dist <- matrix( c( NA, -400, -600, NA, NA, -120, NA, NA, NA ), nrow = 3, ncol = 3 )
mat.scaledranks <- distMat2UnitNormRank(mat.dist)

# With metrics for which higher is closer/more similar, use
# negDistMat2UnitNormRank():
mat.jaccard <- matrix( c( NA, 0.2, 0.3, NA, NA, 0.1, NA, NA, NA ), nrow = 3, ncol = 3 )
mat.srjaccard <- negDistMat2UnitNormRank(mat.jaccard)

# This also works:
mat.srjaccard <- distMat2UnitNormRank(mat.jaccard, lower_is_closer=FALSE)

}
\seealso{
\code{\link{distMat2Rank}()}
}
