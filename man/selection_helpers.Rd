\name{selection_helpers}
\alias{mv_MR_susie_sel_joint}
\alias{mv_MR_susie_sel_joint_score}
\alias{mv_MR_susie_sel_sequential}
\title{Role-localization helpers}
\description{Lower-level Gaussian and binary score selectors, and a sequential-selection comparison helper.}
\usage{
mv_MR_susie_sel_joint (corr, betaZY, se_betaZY, gamma_hat, se_gamma = NULL, n2, L = 10, exposure_prior_mass = 0.5, 
    coverage = 0.95, me_correct = FALSE, sigma_exp = NULL, sigma_res = NULL, n1 = NULL, purity = 0.5, 
    ambiguous = c("unresolved", "invalid", "exposure"))

mv_MR_susie_sel_joint_score (corr, z, gamma_hat, L = 10, exposure_prior_mass = 0.5, coverage = 0.95, purity = 0.5, 
    n_eff = NULL, ambiguous = c("unresolved", "invalid", "exposure"))

mv_MR_susie_sel_sequential (corr, betaZY, se_betaZY, gamma_hat, se_gamma = NULL, n2, L = 10, exposure_prior_mass = 0.5, 
    coverage = 0.95, me_correct = FALSE, sigma_exp = NULL, sigma_res = NULL, n1 = NULL, purity = 0.5, 
    ambiguous = c("unresolved", "invalid", "exposure"))
}
\arguments{
\item{corr}{Signed LD correlation matrix.}
\item{betaZY}{Marginal variant-outcome coefficients.}
\item{se_betaZY}{Standard errors of the outcome coefficients.}
\item{gamma_hat}{Joint variant-exposure coefficient matrix, with exposures in columns.}
\item{se_gamma}{Optional exposure coefficient standard-error matrix.}
\item{n2}{Outcome sample-size parameter for the Gaussian sufficient statistics.}
\item{L}{Maximum number of single-effect components.}
\item{exposure_prior_mass}{Total prior mass assigned to the exposure block.}
\item{coverage}{Requested credible-set probability.}
\item{me_correct}{Compatibility argument. Selection retains the coherent augmented Gram matrix.}
\item{sigma_exp}{Optional exposure-error correlation matrix.}
\item{sigma_res}{Optional residual covariance supplied for compatibility.}
\item{n1}{Optional exposure sample-size parameter supplied for compatibility.}
\item{purity}{Minimum absolute within-set correlation for retained credible sets.}
\item{ambiguous}{Handling of mixed exposure-variant sets. The common method uses unresolved roles.}
\item{z}{Signed logistic Wald statistics for score localization.}
\item{n_eff}{Effective sample size. NULL requests the large-sample score convention.}
}
\value{A list of selected exposures and invalid variants, role-specific and mixed credible sets, inclusion probabilities and convergence information.}
\details{These functions expose localization steps for compatibility and analysis comparisons. Use \code{make_mr_susie} for the complete common procedure with its point estimates, corrected-score inference and diagnostics. The sequential helper is a comparator that fits variant and exposure roles in successive stages.}
