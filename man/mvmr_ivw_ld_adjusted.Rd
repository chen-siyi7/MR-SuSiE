\name{mvmr_ivw_ld_adjusted}
\alias{mvmr_ivw_ld_adjusted}
\title{LD-adjusted causal-effect estimation}
\usage{
mvmr_ivw_ld_adjusted(dat, invalid_snps = integer(0),
  which_R = c("in", "ref"), me = TRUE, sigma_exp = NULL)
}
\arguments{
\item{dat}{A list containing \code{num_exposures}, \code{betaZX},
\code{betaZY}, \code{se_betaZY}, \code{seZX}, \code{corr_in}, and
\code{corr_ref}. Optional elements are \code{sigma_exp} and \code{overlap}.}
\item{invalid_snps}{Indices of variants selected as having direct effects.}
\item{which_R}{Use the in-sample or reference-panel LD matrix.}
\item{me}{Apply the corrected-score measurement-error adjustment.}
\item{sigma_exp}{Cross-exposure sampling correlation matrix.}
}
\value{A list with causal coefficients, standard errors, p-values, adjusted
q-values, direct-effect estimates, and identification diagnostics.}
\description{Fits an LD-aware corrected score while retaining selected direct
variant effects as nuisance directions.}
