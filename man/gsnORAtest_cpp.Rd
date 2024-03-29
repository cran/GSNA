% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gsnORAtest_cpp}
\alias{gsnORAtest_cpp}
\title{gsnORAtest_cpp}
\usage{
gsnORAtest_cpp(l, bg, geneSetCollection)
}
\arguments{
\item{l}{(required) A character vector containing a list of gene identifiers. These are generally differentially
expressed genes either genes significantly up or significantly down, but they can also be a list of genes that came
out of a genetic screen, gene loci with differential chromatin accessibility generated by ATAC-Seq data, lists of genes
from GWAS, etc. The order of the genes is unimportant.}

\item{bg}{(required) A character vector containing a list of gene identifiers corresponding to the total background of
observable genes.}

\item{geneSetCollection}{(required) A list of gene sets, in which the gene sets are character vectors containing
gene symbols, and the list names are the corresponding gene set identifiers. NOTE: This must be a list, not a
\code{tmod} object. It is trivial to extract such a list from a \code{tmod} object, however. The
\code{$MODULES2GENES} field of the \code{tmod} object contains a suitable list.}
}
\value{
A data frame containing the results of overrepresentation analysis.
\itemize{
\item{\emph{ID}: the gene set identifiers.}
\item{\emph{a}: the number of genes observed in the background but not in \emph{l} or the queried gene set.}
\item{\emph{b}: the number of observed genes in \emph{l} but not the queried gene set.}
\item{\emph{c}: the number of observed genes in the queried gene set but not \emph{l} and}
\item{\emph{d}: the number of observed genes in both \emph{l} and the queried gene set, i.e. the overlap.}
\item{\emph{N}: the number of observed genes the queried gene set.}
\item{\emph{Enrichment}: The fold overrepresentation of genes in the overlap set \emph{d} calculated as:
\deqn{E = (d / (c+d)) / ((b+d)/(a+b+c+d))}
}
\item{\emph{P_2S}: 2-sided Fisher \emph{p}-value. (\emph{NOT} log-transformed.)}
\item{\emph{P_1S}: 1-sided Fisher \emph{p}-value. (\emph{NOT} log-transformed.)}
}
}
\description{
This function performs ORA analysis and returns a data.frame containing various statistics including
fold enrichment, and 1 and 2-tailed p-values. (see details)
}
\details{
This is the main workhorse function for the ORA test in the \code{GSNA} package, however, it performs
no filtering of the output data set, nor \emph{p}-value adjustment, and most users of the package will want to use
\code{gsnORAtest()} function instead, which calculates adjusted \emph{p}-values, filters the output data for
significance, and can include a \code{Title} field in the output data.frame.
}
\examples{

library(GSNA)

# From a differential expression data set, we can generate a
# subset of genes with significant differential expression,
# up or down. Here we will extract genes with significant
# negative differential expression with
# avg_log2FC < 0 and p_val_adj <= 0.05 from **Seurat** data:

sig_DN.genes <-
   toupper( rownames(subset( Bai_CiHep_v_Fib2.de,
                      avg_log2FC < 0  & p_val_adj < 0.05 )) )

# Using all the genes in the differential expression data set,
# we can obtain a suitable background:
bg <- rownames( Bai_CiHep_v_Fib2.de )

# Next we need a gene set collection in the form of a list of
# character vectors. We can convert the **Bai_gsc.tmod** object
# included in the sample data to such a list:
Bai.gsc <- tmod2gsc( Bai_gsc.tmod )

# Now, we can do a overrepresentation analysis search on this
# data using **Bai.gsc**:
sig_DN.gsnora <- gsnORAtest_cpp( l = sig_DN.genes,
                                 bg = bg,
                                 geneSetCollection = Bai.gsc )

}
\seealso{
\code{\link{gsnORAtest}}
}
