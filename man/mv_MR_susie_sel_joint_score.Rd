\name{mv_MR_susie_sel_joint_score}
\alias{mv_MR_susie_sel_joint_score}
\title{Joint selection for a binary outcome}
\usage{
mv_MR_susie_sel_joint_score(corr, z, gamma_hat, L = 10,
  exposure_prior_mass = 0.5, coverage = 0.95, purity = 0.5,
  n_eff = NULL, ambiguous = c("unresolved", "invalid", "exposure"))
}
\arguments{
\item{corr}{Variant correlation or LD matrix.}
\item{z}{Signed marginal outcome score or Wald statistics.}
\item{gamma_hat}{Matrix of variant-exposure coefficients.}
\item{L}{Maximum number of SuSiE single effects.}
\item{exposure_prior_mass}{Prior mass assigned to the exposure block.}
\item{coverage}{Credible-set coverage.}
\item{purity}{Minimum absolute within-set correlation.}
\item{n_eff}{Optional effective sample size.}
\item{ambiguous}{How mixed variant-exposure credible sets are classified.}
}
\value{A selection result with exposure, invalid-variant, posterior, and
ambiguity diagnostics.}
\description{Uses the SuSiE RSS likelihood for signed score statistics without
reconstructing a continuous outcome.}
