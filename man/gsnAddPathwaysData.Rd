% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsnAddPathwaysData.R
\name{gsnAddPathwaysData}
\alias{gsnAddPathwaysData}
\alias{gsnAddPathwayData}
\title{gsnAddPathwaysData}
\usage{
gsnAddPathwaysData(
  object,
  pathways_data,
  type = NULL,
  id_col = NULL,
  stat_col = NULL,
  sig_order = NULL,
  stat_col_2 = NULL,
  sig_order_2 = NULL,
  n_col = NULL
)

gsnAddPathwayData(
  object,
  pathways_data,
  type = NULL,
  id_col = NULL,
  stat_col = NULL,
  sig_order = NULL,
  stat_col_2 = NULL,
  sig_order_2 = NULL,
  n_col = NULL
)
}
\arguments{
\item{object}{A GSNData object.}

\item{pathways_data}{A data.frame containing the results of pathways analysis.}

\item{type}{(optional) A character vector of length 1 indicating the type of pathways data being added to the
GSNData object. This can be \code{'cerno'}, \code{'gsea'}, \code{'gsnora'}, or other arbitrary types. If not
explicitly indicated, the method attempts to examine the column names of the data.frame in order to determine
what kind of import to perform, then calls other methods for the actual import. For \code{'cerno'}, \code{'gsea'},
and \code{'gsnora'}, the actual import is performed by methods specifically designed for CERNO and GSEA import.
Otherwise a method for generic import is used.}

\item{id_col}{(optional) A character vector of length 1 indicating the name of the column used as a key for gene
sets/modules. This corresponds to the ID field of \code{tmod} objects, or the names of vectors in a list vectors
gene sets/modules, both of which can be used as a geneSetCollection argument in building gene set networks. In
the case of CERNO and GSEA data sets, there are preset values for \code{id_col}, but in the case of generic
import, the import method attempts to guess. If an ID cannot be inferred, then an error is thrown.}

\item{stat_col}{(optional) A character vector of length 1 indicating the name of the column used as a statistic
to evaluate the quality of pathways results. This is generally a \emph{p}-value of some sort. In the case of CERNO
and GSEA data sets, there are preset values for \code{stat}, but in the case of generic import, the import
method attempts to guess.}

\item{sig_order}{(optional) Either \code{'loToHi'} or \code{'hiToLo'} depending on the statistic used to
evaluate pathways results. For \emph{p}-values, this should be \code{'loToHi'}.}

\item{stat_col_2}{(optional) A character vector of length 1 indicating the name of the column used as a second
statistic to evaluate pathway result quality. Used in 2-color networks.}

\item{sig_order_2}{(optional) Either \code{'loToHi'} or \code{'hiToLo'} depending on \code{stat_col_2}. Used
in 2-color networks.}

\item{n_col}{(optional) Specifies the column containing the number of genes in the gene set. Generally, this is the number
of genes in the gene set that are attested in an expression data set.}
}
\value{
This returns a GSNData object containing imported pathways data.
}
\description{
Add pathways search data to a GSNData object.

\code{gsnAddPathwayData()} is a synonym of \code{\link{gsnAddPathwaysData}()}, included to support
old code. Use \code{\link{gsnAddPathwaysData}()} for new code. \code{gsnAddPathwayData()} will be deprecated in]
a future version.
}
\details{
Pathways data are used by the \code{assignSubnets()} function, which organizes subnets on the basis
of this statistic. If \code{sig_order} is \code{'loToHi'}, and the evaluation statistic (\code{'stat'}) is a
\emph{p}-value, then the first node in each subnet will be the node with the lowest \emph{p}-value, for example. This
ordering is not an absolute requirement.

This is provided to simplify workflows and facilitate imports that can identify and handle multiple types of
pathways data, but also the CERNO, GSEA, GSNORA, and generic import methods can be used directly
( \code{\link{gsnImportCERNO}}, \code{\link{gsnImportGSEA}}, \code{\link{gsnImportGSNORA}},
and \code{\link{gsnImportGenericPathways}}).

Notes: These import handlers perform checks on the provided pathways data to verify that
all gene set IDs in the genePresenceAbsence matrix are present in the ID column of the pathways data. An error
is thrown if all gene set IDs in the genePresenceAbsense are not present in the pathways ID column. On the other
hand, if there are gene set IDs present in the pathways data that are absent from the genePresenceAbsence matrix,
then these methods emit a warning.
}
\examples{

# In this example, we generate a gene set network from CERNO example
# data. We begin by subsetting the CERNO data for significant results:
sig_pathways.cerno <- subset( Bai_CiHep_DN.cerno, adj.P.Val <= 0.05 )

# Now create a gene set collection containing just the gene sets
# with significant CERNO results, by subsetting Bai_gsc.tmod using
# the gene set IDs as keys:
sig_pathways.tmod <- Bai_gsc.tmod[sig_pathways.cerno$ID]

# And obtain a background gene set from differential expression data:
background_genes <- toupper( rownames( Bai_CiHep_v_Fib2.de ) )

# Build a gene set network:
sig_pathways.GSN <-
   buildGeneSetNetworkJaccard(geneSetCollection = sig_pathways.tmod,
                              ref.background = background_genes )

# Now import the CERNO data:
sig_pathways.GSN <- gsnAddPathwaysData( sig_pathways.GSN,
                                        pathways_data = sig_pathways.cerno )

}
\seealso{
\code{\link{gsnImportCERNO}}
\code{\link{gsnImportGSEA}}
\code{\link{gsnImportGSNORA}}
\code{\link{gsnImportGenericPathways}}
}
