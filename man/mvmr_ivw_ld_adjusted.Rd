\name{mvmr_ivw_ld_adjusted}
\alias{mvmr_ivw_ld_adjusted}
\title{LD-adjusted IVW comparison helper}
\description{Fits the nuisance-adjusted LD-aware IVW comparison estimator.}
\usage{
mvmr_ivw_ld_adjusted (dat, invalid_snps = integer(0), which_R = c("in", "ref"), me = TRUE, sigma_exp = NULL)
}
\arguments{
\item{dat}{List containing exposure and outcome coefficients, their standard errors, exposure count, LD matrices, and optional exposure-error covariance inputs.}
\item{invalid_snps}{Indices of the supplied invalid variants.}
\item{which_R}{Which LD matrix to use: corr_in or corr_ref.}
\item{me}{Whether to apply the exposure-error adjustment.}
\item{sigma_exp}{Optional exposure-error correlation matrix overriding the value in dat.}
}
\value{A list containing identification status, coefficient estimates, uncertainty summaries and numerical diagnostics.}
\details{This comparison helper is distinct from the covariance-aware posterior point map returned by the full MR-SuSiE model.}
