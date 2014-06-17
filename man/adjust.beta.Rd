\name{adjust.beta}
\alias{adjust.beta}
\title{Adjust Beta for cell mixture. Returns beta of the average cell-type.}
\usage{
adjust.beta(B, top_n = 500, mc.cores = 2, cell.coefs = NULL,
  est.only = FALSE)
}
\arguments{
  \item{B}{matrix of beta values (n_probes * n_samples)
  rownames must be the illumina cg id or the
  chrom:position}

  \item{top_n}{number of probes from cell.coefs to use}

  \item{mc.cores}{number of cores to use for
  parallelization}

  \item{cell.coefs}{path to tab-delimited file containind
  DMR-probes for the cell types of interest.
  system.file("extdata", "houseman-dmrs.txt",
  package="celltypes450")}

  \item{est.only}{if TRUE, only estimate the coefficients
  for the cell-types without actually adjusting.}
}
\value{
matrix of adjusted beta values (or data.frame of cell coefs
if) \code{est.only} is FALSE
}
\description{
Can also return the estimates of each cell when
est.only=TRUE
}

