% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsnImportDAVID.R
\name{gsnImportDAVID}
\alias{gsnImportDAVID}
\title{gsnImportDAVID}
\usage{
gsnImportDAVID(
  object,
  pathways_data = NULL,
  filename = NULL,
  id_col = NULL,
  stat_col = NULL,
  sig_order = NULL,
  n_col = NULL,
  sep = "\\t"
)
}
\arguments{
\item{object}{A GSNData object.}

\item{pathways_data}{An (optional) data.frame containing the results of DAVID analysis. (Either this or the
\code{filename} argument must be set. Such a data.frame can be obtained by using the \code{read_david_data_file()}
function to parse a DAVID "Functional Annotation Chart" or "Functional Annotation Cluster" results text file with
the default options (output = "flat", redundant = FALSE, sep = "\\t").}

\item{filename}{An (optional) filename for data sets read from a text file containing DAVID results. This is ignored
if the \code{pathways_data} argument is set.}

\item{id_col}{(optional) A character vector of length 1 indicating the name of the column used as a key for gene
sets or modules. This is normally the \code{Term} field of DAVID data, which must be the same as the names of gene
sets in the gene set collection specified with the \code{geneSetCollection} argument used when building the gene
set network. By default this value is \code{'Term'}. The IDs must correspond to the names of the gene sets provided,
or an error will be thrown.}

\item{stat_col}{(optional) A character vector of length 1 indicating the name of the column used as a statistic
to evaluate the quality of pathways results. The function scans through possible \code{stat_col} values
("FDR", "Bonferroni", "Benjamini", "PValue" ), and uses the first one it finds.}

\item{sig_order}{(optional) Either \code{'loToHi'} (default) or \code{'hiToLo'} depending on the statistic used to
evaluate pathways results.}

\item{n_col}{(optional) Specifies the column containing the number of genes in the gene set. Generally, this is the number
of genes in the gene set that are attested in an expression data set (Defaults to 'Count', if that is present, otherwise}

\item{sep}{A separator for text file import, defaults to "\\t". Ignored if \code{filename} is not specified.}
}
\value{
This returns a GSNData object containing imported pathways data.
}
\description{
Add DAVID search data to a GSNData object, as generated by the the DAVID web application
(\href{https://david.ncifcrf.gov/}{https://david.ncifcrf.gov/}) output using either the
"Functional Annotation Chart" or "Functional Annotation Cluster" results output options. The data set
can be either in the form of a data.frame or specified as import from an output text file. (See Details below)
}
\details{
Note: An error is thrown if all gene set IDs in the genePresenceAbsense are not present in the GSEA \code{NAME} column.
However, if there are gene set IDs present in the pathways data that are absent from the \code{$genePresenceAbsence}
matrix, then this method emits a warning. It also checks for the standard GSEA data set column names, and if some are
missing, it will emit a warning.
}
\seealso{
\code{\link{gsnAddPathwaysData}}
\code{\link{gsnImportCERNO}}
\code{\link{gsnImportGSNORA}}
\code{\link{gsnImportGenericPathways}}
}
