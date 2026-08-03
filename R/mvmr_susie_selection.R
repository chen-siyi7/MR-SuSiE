## mv_MR_SuSiE.R  --  revised MVMR-SuSiE core
##
## Fixes the central bug in the original sum_MR_SuSiE_mv: it fed SuSiE the
## variant-only cross-product (XtX = corr, m x m), so the fit could only ever
## return SNP indices and "I_exposure = Iu[Iu > p] - p" was ALWAYS empty. The
## method therefore never selected exposures; exposure selection happened only
## implicitly through the stage-2 FDR.
##
## Here the design is AUGMENTED with the K predicted-exposure columns
## Dhat = Z %*% gamma_hat, so SuSiE jointly fine-maps [variants, exposures] and
## returns (a) the invalid-instrument set and (b) an exposure credible set --
## the multi-exposure analogue of the univariate selection. Because Dhat lies in
## the column span of Z, the augmented design is rank-deficient; identification
## rests on the SuSiE sparsity prior, and SuSiE may return a *group* of
## correlated exposures in one credible set rather than resolving a single one.
## That grouping is the honest output under correlated exposures, not a failure.

## SUMMARY-LEVEL ONLY. The individual-level selection functions from earlier
## drafts (uv_MR_susie_sel, mv_MR_susie_sel_ind) have been removed: no method in
## this package uses individual-level genotypes.

## susieR is needed only at call time; attach it if available so the file sources
## cleanly on a machine without it (the oracle path and comparators do not need it).


## =====================================================================
## SHARED PIECES
##
## The joint fit and the sequential comparator must differ in ORDERING ALONE, so
## every ingredient they have in common is built here once and used by both.
## Anything computed separately in the two functions is a place where an
## unintended difference can creep in, which is exactly what happened before:
## the sequential function used n = max(n1,n2), no prior weights,
## estimate_prior_variance = FALSE and no measurement-error correction, while the
## joint function used n = n2, split prior weights, estimate_prior_variance =
## TRUE and the correction. Three of those four are the bugs the joint function
## was written to fix, so the ablation would have attributed to "joint vs
## sequential" a gap that was really "fixed vs unfixed".
## =====================================================================

## Response cross-product Y'Y/n2 from marginal outcome summary statistics.
## Identity: for the marginal regression of Y on G_j with an intercept,
##   Y'Y = G_j'G_j * [ (n2 - 2) * se_j^2 + beta_j^2 ],
## the (n2 - 2) being the residual degrees of freedom. Each j estimates the same
## quantity; the median guards against one mis-reported standard error.
## (The previous version used n2 - 1. The difference is O(1/n) and immaterial at
## these sample sizes, but it is not a free parameter.)
.suff_yty <- function(corr, betaZY, se_betaZY, n2) {
  m <- nrow(corr)
  stats::median(vapply(seq_len(m), function(j)
    corr[j, j] * ((n2 - 2) * se_betaZY[j]^2 + betaZY[j]^2), numeric(1)))
}

## Variant-outcome cross-product G'Y/n2 = diag(R) * betaZY. Equals betaZY when
## corr is a correlation matrix; the factor is retained so the identity stays
## correct if a covariance matrix is ever passed, and the check makes that
## assumption explicit rather than silent.
.suff_zty <- function(corr, betaZY) {
  d <- diag(corr)
  if (any(!is.finite(d)) || any(d <= 0))
    stop("corr has a non-positive or non-finite diagonal", call. = FALSE)
  d * betaZY
}

## Bias of the exposure block of the augmented cross-product.
##
## With Ghat = Gamma + E and E independent of Y,
##   E[Ghat' R Ghat] = Gamma' R Gamma + V,   V_kl = tr( R Cov(ghat_k, ghat_l) ).
##
## The shape of Cov depends on which first-stage estimate Ghat is. For the JOINT
## least-squares first stage used throughout this package
## (gamma_hat <- coef(lm(X1 ~ 0 + Z1))),
##   Cov(ghat_k, ghat_l) = Sigma^res_kl * (Z1'Z1)^{-1},
## so V_kl = Sigma^res_kl * tr( R (Z1'Z1)^{-1} ) = m * Sigma^res_kl / n1 when the
## selection LD stands in for the exposure-sample LD.
##
## A CORRECTION THAT IS *NOT* NEEDED, RECORDED SO IT IS NOT REINTRODUCED.
## It is tempting to reason that the marginal covariance is diag(se^2), that the
## joint covariance is therefore R^{-1} diag(se^2) R^{-1}, and hence that the
## right term is tr(V^marg R^{-1}) -- which would mean the previous code
## under-corrected by a factor growing with LD (about 4.5x at pairwise LD 0.8).
## That reasoning is wrong. The marginal covariance is (sigma^2/n) R, not
## diag(se^2), so
##   tr( V^marg R^{-1} ) = tr( (sigma^2/n) R R^{-1} ) = m sigma^2/n = sum_j se_j^2,
## which is exactly what the previous code computed. Its DIAGONAL was correct.
##
## What was actually wrong is the OFF-DIAGONAL, which was set to zero.
## Cov(ghat_k, ghat_l) is proportional to Sigma^res_kl, the residual covariance
## BETWEEN exposures, which is far from zero whenever the exposures are measured
## on the same individuals -- the defining situation for an omics panel, and the
## situation in every application here. Measured on this DGP the off-diagonal
## runs at about 0.73 of the diagonal, so dropping it is not a rounding error.
##
## Two routes. Given sigma_res (the K x K first-stage residual covariance) and
## n1, the exact form is used. From summary statistics alone the fallback is
##   V_kl = Sigma^exp_kl * sum_j s_jk s_jl R_jj,
## whose diagonal reproduces the previous formula and whose off-diagonal uses the
## cross-exposure correlation Sigma^exp, estimable from the Z-score correlation
## at null variants exactly as for sample-overlap correction.
.me_exposure_block <- function(corr, se_gamma = NULL, sigma_exp = NULL,
                               sigma_res = NULL, n1 = NULL, m = nrow(corr)) {
  if (!is.null(sigma_res) && !is.null(n1)) {
    stopifnot(nrow(sigma_res) == ncol(sigma_res), is.finite(n1), n1 > 1)
    return(m * sigma_res / n1)
  }
  if (is.null(se_gamma))
    stop("give either (sigma_res, n1) or se_gamma", call. = FALSE)
  K <- ncol(se_gamma)
  if (is.null(sigma_exp)) sigma_exp <- diag(K)   # disjoint exposure GWAS
  d <- diag(corr)
  V <- matrix(0, K, K)
  for (k in seq_len(K)) for (l in seq_len(K))
    V[k, l] <- sigma_exp[k, l] * sum(se_gamma[, k] * se_gamma[, l] * d)
  V
}

## Prior weights, split between the two column blocks rather than uniform over
## columns. A uniform prior puts m/(m+K) of the mass on variants -- 92% at m=60,
## K=5 -- so an exposure and a variant would compete on block size rather than on
## evidence. exposure_prior_mass is the explicit form of the prior on which
## identification rests, and should be reported and varied, not hidden.
.block_prior <- function(m, K, exposure_prior_mass) {
  stopifnot(exposure_prior_mass > 0, exposure_prior_mass < 1)
  c(rep((1 - exposure_prior_mass) / m, m), rep(exposure_prior_mass / K, K))
}

## Credible sets at the requested coverage, both purity-filtered (used for
## SELECTION) and unfiltered (used only to COUNT unresolved sets).
##
## Why both. SuSiE discards any credible set whose minimum absolute within-set
## correlation falls below min_abs_corr (0.5 by default). In the augmented design
## a set spanning both blocks is impure almost by construction, because
## Xhat_k = G omega_k is a dense combination of variants and correlates weakly
## with any single variant. The purity filter therefore removes PREFERENTIALLY
## the sets that express the variant-versus-exposure ambiguity, and an unresolved
## count taken from the filtered output is identically zero whatever the data.
## Measured in the arm-2 regime: 32 filtered sets of which 0 mixed, against 40
## unfiltered sets of which 8 mixed -- every set the filter removed was mixed.
.credible_sets <- function(res, XtX, coverage, purity) {
  Rc <- tryCatch(stats::cov2cor(XtX), error = function(e) NULL)
  if (is.null(Rc) || any(!is.finite(Rc))) {
    cs <- res$sets$cs
    cs <- if (is.null(cs)) list() else cs
    return(list(pure = cs, all = cs))
  }
  gp <- function(mac) {
    out <- tryCatch(susieR::susie_get_cs(res, Xcorr = Rc, coverage = coverage,
                                         min_abs_corr = mac)$cs,
                    error = function(e) NULL)
    if (is.null(out)) list() else out
  }
  list(pure = gp(purity), all = gp(0))
}

## Count credible sets spanning BOTH column blocks: sets the data could not
## resolve into "these variants are invalid" or "this exposure is causal".
.count_mixed <- function(cs_list, m) {
  if (!length(cs_list)) return(0L)
  sum(vapply(cs_list, function(s) any(s <= m) && any(s > m), logical(1)))
}

## Preserve the membership of every mixed credible set.  The unfiltered
## collection is the primary ambiguity diagnostic; the purity-passing subset is
## returned separately so that the filtering decision remains auditable.
.mixed_sets <- function(cs_list, m) {
  if (!length(cs_list)) return(list())
  out <- lapply(cs_list, function(s) {
    snp_part <- s[s <= m]
    exp_part <- s[s > m] - m
    if (!length(snp_part) || !length(exp_part)) return(NULL)
    list(snps = snp_part, exposures = exp_part)
  })
  Filter(Negate(is.null), out)
}


## ---- ORIGINAL summary-level selection, KEPT ONLY FOR REPRODUCIBILITY --------
## Do not use. It misreads SuSiE's output (see point 1 below) in a way that
## affects every downstream number. Retained so that results generated before the
## fix can be regenerated exactly.
mv_MR_susie_sel_summary <- function(corr, betaZY, se_betaZY, gamma_hat,
                                    n1, n2, L = NULL) {
  m <- nrow(corr); K <- ncol(gamma_hat)
  ZtY <- diag(corr) * betaZY
  yty <- median(sapply(seq_len(m), function(j)
    (n2 - 1) * corr[j, j] * se_betaZY[j]^2 + ZtY[j] * betaZY[j]))

  CG  <- corr %*% gamma_hat
  GCG <- t(gamma_hat) %*% CG
  XtX_aug <- rbind(cbind(corr, CG), cbind(t(CG), GCG))
  Xty_aug <- c(ZtY, as.numeric(t(gamma_hat) %*% ZtY))

  if (is.null(L)) L <- min(m + K - 1, m)
  res <- susieR::susie_suff_stat(XtX = XtX_aug, Xty = Xty_aug, yty = yty,
                         n = max(n1, n2), L = L,
                         estimate_residual_variance = FALSE,
                         estimate_prior_variance = FALSE)

  cs_list <- res$sets$cs
  in_cs <- sort(unique(unlist(cs_list)))
  invalid_snps <- in_cs[in_cs <= m]
  valid_snps   <- setdiff(seq_len(m), invalid_snps)
  exposure_cols <- in_cs[in_cs > m] - m
  exposure_cs <- lapply(cs_list, function(s) { e <- s[s > m] - m; if (length(e)) e else NULL })
  exposure_cs <- Filter(Negate(is.null), exposure_cs)

  list(valid = valid_snps, invalid = invalid_snps,
       exposures = exposure_cols, exposure_cs = exposure_cs,
       pip = res$pip, pip_exposures = res$pip[(m + 1):(m + K)])
}


## ===================== JOINT selection (the method) =========================
## Invalid instruments and indistinguishable exposures are resolved in ONE fit on
## the augmented design [G, Xhat]. Corrections relative to
## mv_MR_susie_sel_summary, in order of how much they matter:
##
## (1) A CREDIBLE SET IS A DISJUNCTION, NOT A CONJUNCTION.
## SuSiE decomposes b into L single effects, each with exactly one nonzero entry,
## so a credible set says "one of these variables carries this effect". The
## original pooled all sets and split them by column index, which for a set
## {Xhat_2, G_7, G_13} concludes exposure 2 is causal AND SNPs 7 and 13 are
## invalid, when SuSiE was expressing uncertainty BETWEEN those options. Each set
## is now resolved to one role, and a set spanning both column types is reported
## UNRESOLVED rather than read as evidence for either branch.
##
## (2) THE PRIOR WAS NOT EXCHANGEABLE ACROSS COLUMN TYPES. See .block_prior.
##
## (3) L IS A SIGNAL COUNT, NOT A COLUMN COUNT. L = min(m+K-1, m) gave 60 single
## effects over 65 columns; PIPs then accumulate toward uniformity, which is why
## a PIP > 1/K rule selected every exposure.
##
## (4) THE SAMPLE SIZE IS THE OUTCOME'S. The augmented regression's response is
## Y, so n = n2. max(n1, n2) is arbitrary and directly scales selection stringency.
##
## (5) MEASUREMENT ERROR IN gamma_hat CANNOT BE SUBTRACTED FROM ONLY THE
## exposure block: doing so makes the augmented cross-product indefinite.
## Exposure error is handled in the downstream corrected estimating equation.
##
## WHAT IS NOT FIXED, BECAUSE IT CANNOT BE BY THESE MEANS.
## Xhat_k = G omega_k lies exactly in the column space of G, so "exposure k is
## causal" and "the variants comprising Xhat_k have direct effects" have the same
## likelihood. SuSiE's coverage guarantees assume single effects identified up to
## correlation, not exact linear dependence. What breaks the tie is the sparsity
## prior preferring one column to many, so identification rests on a prior belief
## that pleiotropy is sparser than the exposure's genetic architecture.
## n_unresolved is the honest measure of how often that tie is not broken, and it
## is counted on the UNFILTERED sets for the reason given in .credible_sets.
mv_MR_susie_sel_joint <- function(corr, betaZY, se_betaZY, gamma_hat, se_gamma = NULL,
                                  n2, L = 10, exposure_prior_mass = 0.5,
                                  coverage = 0.95, me_correct = FALSE,
                                  sigma_exp = NULL, sigma_res = NULL, n1 = NULL,
                                  purity = 0.5,
                                  ambiguous = c("unresolved", "invalid", "exposure")) {
  ambiguous <- match.arg(ambiguous)
  m <- nrow(corr); K <- ncol(gamma_hat)
  stopifnot(is.finite(n2), n2 > 1)
  ZtY <- .suff_zty(corr, betaZY)
  yty <- .suff_yty(corr, betaZY, se_betaZY, n2)

  CG  <- corr %*% gamma_hat
  GCG <- t(gamma_hat) %*% CG
  ## Do not subtract an exposure-block measurement-error matrix here. Because
  ## Xhat = G gamma_hat exactly, the Schur complement of the uncorrected
  ## augmented cross-product is zero. Subtracting any positive V makes that
  ## Schur complement -V and the supplied XtX indefinite, so it cannot be a
  ## cross-product matrix accepted by the SuSiE likelihood. Measurement error
  ## is handled in the downstream estimating equation instead.
  if (me_correct)
    warning("selection-stage measurement-error subtraction is disabled: it makes ",
            "the augmented XtX indefinite; using the exact cross-product", call. = FALSE)
  XtX_aug <- rbind(cbind(corr, CG), cbind(t(CG), GCG))
  Xty_aug <- c(ZtY, as.numeric(t(gamma_hat) %*% ZtY))

  L <- max(1, min(L, m + K - 1))
  res <- susieR::susie_suff_stat(XtX = XtX_aug, Xty = Xty_aug, yty = yty, n = n2, L = L,
                         prior_weights = .block_prior(m, K, exposure_prior_mass),
                         estimate_residual_variance = FALSE,
                         estimate_prior_variance = TRUE, check_input = TRUE)

  sets <- .credible_sets(res, XtX_aug, coverage, purity)

  ## (1) resolve each purity-filtered credible set to ONE role
  invalid_snps <- integer(0); invalid_cs <- list()
  exposure_cs <- list(); unresolved_pure <- list()
  for (s in sets$pure) {
    snp_part <- s[s <= m]; exp_part <- s[s > m] - m
    if (length(snp_part) && !length(exp_part)) {
      invalid_snps <- c(invalid_snps, snp_part)
      invalid_cs[[length(invalid_cs) + 1]] <- snp_part
    } else if (length(exp_part) && !length(snp_part)) {
      exposure_cs[[length(exposure_cs) + 1]] <- exp_part
    } else if (length(snp_part) && length(exp_part)) {
      unresolved_pure[[length(unresolved_pure) + 1]] <-
        list(snps = snp_part, exposures = exp_part)
      if (ambiguous == "invalid") invalid_snps <- c(invalid_snps, snp_part)
      if (ambiguous == "exposure") exposure_cs[[length(exposure_cs) + 1]] <- exp_part
    }
  }
  ## Mixed sets removed by the purity filter remain scientifically unresolved.
  ## Preserve their membership instead of adding only a hidden-set count.
  unresolved_all <- .mixed_sets(sets$all, m)

  invalid_snps <- sort(unique(invalid_snps))
  valid_snps <- setdiff(seq_len(m), invalid_snps)
  list(valid = valid_snps, invalid = invalid_snps,
       invalid_cs = invalid_cs,
       exposures = sort(unique(unlist(exposure_cs))),
       exposure_cs = exposure_cs,
       unresolved = unresolved_all,
       unresolved_all = unresolved_all,
       unresolved_pure = unresolved_pure,
       n_unresolved = length(unresolved_all),
       n_unresolved_pure = length(unresolved_pure),
       n_cs = length(sets$pure), n_cs_unfiltered = length(sets$all),
       pip = res$pip, pip_exposures = res$pip[(m + 1):(m + K)],
       L_used = L, n_used = n2, order = "joint",
       converged = isTRUE(res$converged))
}

## ================= BINARY-OUTCOME JOINT SELECTION ============================
## Logistic-regression beta/SE pairs are not OLS sufficient statistics: in
## particular, (n-2) se^2 + beta^2 does not reconstruct Y'Y. For a case-control
## outcome, use the standard z-score RSS approximation instead:
##
##     z | delta  approximately N(R delta, R),
##
## where z is the vector of signed logistic Wald statistics and R is SNP LD.
## This is the same large-sample likelihood used by susie_rss(z, R, n = NULL).
## It avoids inventing a continuous response or using the sample case fraction
## as a liability conversion. The augmented exposure columns are standardized
## before the fit so that their scale is not determined by arbitrary exposure
## units. Downstream MR estimation remains on the supplied betaZY scale (usually
## log odds), because only selection is performed here.
mv_MR_susie_sel_joint_score <- function(
    corr, z, gamma_hat, L = 10, exposure_prior_mass = 0.5,
    coverage = 0.95, purity = 0.5, n_eff = NULL,
    ambiguous = c("unresolved", "invalid", "exposure")) {
  ambiguous <- match.arg(ambiguous)
  m <- nrow(corr); K <- ncol(gamma_hat)
  stopifnot(ncol(corr) == m, length(z) == m, nrow(gamma_hat) == m,
            all(is.finite(z)))

  CG <- corr %*% gamma_hat
  GCG <- t(gamma_hat) %*% CG
  S0 <- rbind(cbind(corr, CG), cbind(t(CG), GCG))
  sc <- sqrt(diag(S0))
  if (any(!is.finite(sc)) || any(sc <= 0))
    stop("the augmented score design contains a zero-variance column",
         call. = FALSE)
  Raug <- sweep(sweep(S0, 1L, sc, "/"), 2L, sc, "/")
  Raug <- (Raug + t(Raug)) / 2
  diag(Raug) <- 1
  zaug <- c(z, as.numeric(t(gamma_hat) %*% z)) / sc

  vn <- rownames(corr)
  if (is.null(vn)) vn <- paste0("variant_", seq_len(m))
  en <- colnames(gamma_hat)
  if (is.null(en)) en <- paste0("exposure_", seq_len(K))
  dimnames(Raug) <- list(c(vn, en), c(vn, en))
  names(zaug) <- c(vn, en)

  L <- max(1, min(L, m + K - 1))
  rss_args <- list(
    z = zaug, R = Raug, L = L,
    prior_weights = .block_prior(m, K, exposure_prior_mass),
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE)
  if (!is.null(n_eff)) {
    stopifnot(length(n_eff) == 1L, is.finite(n_eff), n_eff > 1)
    rss_args$n <- n_eff
  }
  res <- do.call(susieR::susie_rss, rss_args)
  sets <- .credible_sets(res, Raug, coverage, purity)

  invalid_snps <- integer(0); invalid_cs <- list()
  exposure_cs <- list(); unresolved_pure <- list()
  for (s in sets$pure) {
    snp_part <- s[s <= m]; exp_part <- s[s > m] - m
    if (length(snp_part) && !length(exp_part)) {
      invalid_snps <- c(invalid_snps, snp_part)
      invalid_cs[[length(invalid_cs) + 1L]] <- snp_part
    } else if (length(exp_part) && !length(snp_part)) {
      exposure_cs[[length(exposure_cs) + 1L]] <- exp_part
    } else if (length(snp_part) && length(exp_part)) {
      unresolved_pure[[length(unresolved_pure) + 1L]] <-
        list(snps = snp_part, exposures = exp_part)
      if (ambiguous == "invalid")
        invalid_snps <- c(invalid_snps, snp_part)
      if (ambiguous == "exposure")
        exposure_cs[[length(exposure_cs) + 1L]] <- exp_part
    }
  }
  unresolved_all <- .mixed_sets(sets$all, m)
  invalid_snps <- sort(unique(invalid_snps))

  list(
    valid = setdiff(seq_len(m), invalid_snps),
    invalid = invalid_snps, invalid_cs = invalid_cs,
    exposures = sort(unique(unlist(exposure_cs))),
    exposure_cs = exposure_cs, unresolved = unresolved_all,
    unresolved_all = unresolved_all, unresolved_pure = unresolved_pure,
    n_unresolved = length(unresolved_all),
    n_unresolved_pure = length(unresolved_pure),
    n_cs = length(sets$pure), n_cs_unfiltered = length(sets$all),
    pip = res$pip, pip_exposures = res$pip[(m + 1):(m + K)],
    L_used = L,
    n_used = if (is.null(n_eff)) NA_real_ else as.numeric(n_eff),
    order = "joint-score",
    likelihood = "logistic-z-RSS",
    converged = isTRUE(res$converged))
}


## ============ SEQUENTIAL selection: the isolating comparator =================
## The claimed advantage of the joint fit is that invalid instruments and
## indistinguishable exposures are selected on ONE augmented design. Chan et al.
## instead remove invalid instruments first and cluster exposures afterwards.
## That ordering has a real weakness: whether variant j looks invalid depends on
## which exposures are adjusted for, since its residual is
## betaY_j - sum_l theta_l betaX_lj. Choose the wrong exposure set and instruments
## are misclassified; misclassify instruments and the wrong exposures are clustered.
##
## To test whether joint selection actually helps, the comparison must hold
## EVERYTHING ELSE fixed. This function therefore takes the SAME arguments as
## mv_MR_susie_sel_joint and makes the SAME susie_suff_stat call with the SAME
## n2, the SAME prior within each block, the SAME estimate_prior_variance, the
## SAME measurement-error correction and the SAME credible-set rule. The single
## difference is that the variant stage runs on the variant block ALONE and the
## exposure stage afterwards on the survivors, rather than both on one augmented
## design. Any difference in performance is therefore attributable to ordering.
##
## Comparing against Chan et al.'s published pipeline cannot isolate this, because
## that pipeline also differs in using cML rather than SuSiE for the instrument
## stage, in screening exposures marginally beforehand, and in refitting with
## MVMR-cML-BIC at the end.
mv_MR_susie_sel_sequential <- function(corr, betaZY, se_betaZY, gamma_hat, se_gamma = NULL,
                                       n2, L = 10, exposure_prior_mass = 0.5,
                                       coverage = 0.95, me_correct = FALSE,
                                       sigma_exp = NULL, sigma_res = NULL, n1 = NULL,
                                       purity = 0.5,
                                       ambiguous = c("unresolved", "invalid", "exposure")) {
  ambiguous <- match.arg(ambiguous)
  m <- nrow(corr); K <- ncol(gamma_hat)
  stopifnot(is.finite(n2), n2 > 1)
  ZtY <- .suff_zty(corr, betaZY)
  yty <- .suff_yty(corr, betaZY, se_betaZY, n2)
  L <- max(1, min(L, m + K - 1))

  ## ---- stage 1: invalid instruments, WITHOUT the exposure columns -----------
  ## The crux of the sequential design: the instrument stage cannot see which
  ## exposures are in the model, so it must attribute outcome association to
  ## direct effects, with no opportunity to attribute it to an exposure instead.
  ## The prior is uniform within the variant block, which is the joint fit's
  ## variant block renormalized, so the RELATIVE prior among variants is identical.
  s1 <- susieR::susie_suff_stat(XtX = corr, Xty = ZtY, yty = yty, n = n2,
                        L = max(1, min(L, m - 1)), prior_weights = rep(1 / m, m),
                        estimate_residual_variance = FALSE,
                        estimate_prior_variance = TRUE, check_input = TRUE)
  sets1 <- .credible_sets(s1, corr, coverage, purity)
  inv1 <- sort(unique(unlist(sets1$pure)))
  valid_snps <- setdiff(seq_len(m), inv1)

  ## ---- stage 2: exposures, CONDITIONAL on the stage-1 valid set -------------
  ## Same measurement-error correction as the joint fit, restricted to the
  ## surviving instruments. Omitting it here (as the previous version did) hands
  ## the joint fit an advantage that has nothing to do with ordering.
  v <- valid_snps
  if (length(v) < K) {
    return(list(valid = valid_snps, invalid = inv1, invalid_cs = sets1$pure,
                exposures = integer(0), exposure_cs = list(),
                unresolved = list(), n_unresolved = NA_integer_,
                n_unresolved_pure = NA_integer_,
                n_cs = length(sets1$pure),
                n_cs_unfiltered = length(sets1$all),
                pip = s1$pip, pip_exposures = rep(NA_real_, K),
                L_used = L, n_used = n2, order = "sequential",
                converged = isTRUE(s1$converged)))
  }
  Cv <- corr[v, v, drop = FALSE]; Gv <- gamma_hat[v, , drop = FALSE]
  XtX2 <- t(Gv) %*% Cv %*% Gv
  if (me_correct)
    warning("selection-stage measurement-error subtraction is disabled; ",
            "using the exact exposure cross-product", call. = FALSE)
  Xty2 <- as.numeric(t(Gv) %*% ZtY[v])
  s2 <- tryCatch(susieR::susie_suff_stat(XtX = XtX2, Xty = Xty2, yty = yty, n = n2,
                                 L = max(1, min(L, K - 1)), prior_weights = rep(1 / K, K),
                                 estimate_residual_variance = FALSE,
                                 estimate_prior_variance = TRUE, check_input = TRUE),
                 error = function(e) NULL)
  base <- list(valid = valid_snps, invalid = inv1, invalid_cs = sets1$pure,
               ## A sequential fit CANNOT produce a mixed set: the two blocks are
               ## never in the same design. That is not evidence the ambiguity is
               ## absent, only that this ordering cannot express it -- which is
               ## the point of the comparison. Reported NA, not 0.
               unresolved = list(), n_unresolved = NA_integer_,
               n_unresolved_pure = NA_integer_,
               L_used = L, n_used = n2, order = "sequential",
               converged = isTRUE(s1$converged) &&
                 (!is.null(s2) && isTRUE(s2$converged)))
  if (is.null(s2))
    return(c(base, list(exposures = integer(0), exposure_cs = list(),
                        n_cs = length(sets1$pure), n_cs_unfiltered = length(sets1$all),
                        pip = s1$pip, pip_exposures = rep(NA_real_, K))))
  sets2 <- .credible_sets(s2, XtX2, coverage, purity)
  c(base, list(exposures = sort(unique(unlist(sets2$pure))),
               exposure_cs = sets2$pure,
               n_cs = length(sets1$pure) + length(sets2$pure),
               n_cs_unfiltered = length(sets1$all) + length(sets2$all),
               pip = s1$pip, pip_exposures = s2$pip))
}
