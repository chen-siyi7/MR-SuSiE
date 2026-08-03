\name{mv_MR_susie_sel_sequential}
\alias{mv_MR_susie_sel_sequential}
\title{Sequential selection comparator}
\usage{
mv_MR_susie_sel_sequential(corr, betaZY, se_betaZY, gamma_hat,
  se_gamma = NULL, n2, L = 10, exposure_prior_mass = 0.5,
  coverage = 0.95, me_correct = FALSE, sigma_exp = NULL,
  sigma_res = NULL, n1 = NULL, purity = 0.5,
  ambiguous = c("unresolved", "invalid", "exposure"))
}
\arguments{
\item{corr}{Variant correlation or LD matrix.}
\item{betaZY}{Marginal variant-outcome coefficients.}
\item{se_betaZY}{Standard errors for \code{betaZY}.}
\item{gamma_hat}{Matrix of variant-exposure coefficients.}
\item{se_gamma}{Optional standard errors for \code{gamma_hat}.}
\item{n2}{Outcome GWAS sample size.}
\item{L}{Maximum number of SuSiE single effects.}
\item{exposure_prior_mass}{Prior mass assigned to the exposure block.}
\item{coverage}{Credible-set coverage.}
\item{me_correct}{Retained for interface matching.}
\item{sigma_exp}{Optional cross-exposure sampling correlation matrix.}
\item{sigma_res}{Optional first-stage residual covariance matrix.}
\item{n1}{Optional exposure GWAS sample size.}
\item{purity}{Minimum absolute within-set correlation.}
\item{ambiguous}{Compatibility argument.}
}
\value{A selection result analogous to \code{mv_MR_susie_sel_joint}.}
\description{Selects invalid variants first and exposures second. This is an
ordering comparator, not the proposed joint method.}
