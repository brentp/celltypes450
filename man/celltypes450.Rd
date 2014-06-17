\docType{package}
\name{celltypes450}
\alias{celltypes450}
\alias{celltypes450-package}
\title{celltypes450}
\description{
adjust for cell types on 450K methylation data.
}
\details{
currently the main function is

\code{adjust.beta}

It will accept the beta (0, 1) values as a matrix with rows
of probes and columns of samples and return a beta adjusted
for cell type mixture.

It defaults to the cell-types supplied by houseman
(included in this package)
}

