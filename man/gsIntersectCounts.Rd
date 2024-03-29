% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gsIntersectCounts}
\alias{gsIntersectCounts}
\title{gsIntersectCounts}
\usage{
gsIntersectCounts(gs1, gs2, bg_size)
}
\arguments{
\item{gs1}{A character vector representing gene symbols in a gene set.}

\item{gs2}{A character vector representing gene symbols in a second gene set.}

\item{bg_size}{An integer representing the size of the background, i.e. the total number of observable genes.}
}
\value{
A numeric vector of length 4 containing the following 4 elements:
\itemize{
\item{\code{a}: The number of genes in the background that are absent in gs1 and gs2.}
\item{\code{b}: The number of background genes in gs1 but not gs2.}
\item{\code{c}: The number of background genes in gs2 but not gs1.}
\item{\code{d}: The number of background genes in in both gs1 and gs2.}
}
}
\description{
For two character vectors representing two gene sets (gs1 and gs2) and a total number of background
observable genes (that may also be present in gs1 and or gs2 or neither), this function calculates the counts in a
2x2 contingency table for presence and absence of genes in one or both sets or neither. The output of this function
is used as the input for a Fisher test calculation by the GSNA package.
}
\details{
This version of the function may not be retained since it's not currently used. Two alternative versions of the
function in C++ that find the overlap between a \code{std::set<std::string>} and a character vector are used since those versions
are much faster.

NOTE: This function assumes that all genes in gs1 and gs2 are present in the background, so to use this properly, gs1
and gs2 must be filtered to include only genes present in the background.
}
\examples{

library( GSNA )

# We can extract 2 gene sets from the sample data:
Bai.gsc <- tmod2gsc( Bai_gsc.tmod )
M29994.gs = Bai.gsc[['M29994']]
M40825.gs = Bai.gsc[['M40825']]

# Get background gene cout:
bg_gene_count <- nrow( Bai_empty_expr_mat )

# Generate a vector containing the number of contents of
# the 2x2 contingency table:
counts.v <- gsIntersectCounts( gs1 = M29994.gs,
                               gs2 = M40825.gs,
                               bg_size = bg_gene_count )

}
