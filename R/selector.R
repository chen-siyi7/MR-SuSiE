.suff_yty <- function(corr, betaZY, se_betaZY, n2) {
  m <- nrow(corr)
  stats::median(vapply(seq_len(m), function(j)
    corr[j, j] * ((n2 - 2) * se_betaZY[j]^2 + betaZY[j]^2), numeric(1)))
}

.suff_zty <- function(corr, betaZY) {
  d <- diag(corr)
  if (any(!is.finite(d)) || any(d <= 0))
    stop("corr has a non-positive or non-finite diagonal", call. = FALSE)
  d * betaZY
}

.block_prior <- function(m, K, exposure_prior_mass) {
  stopifnot(exposure_prior_mass > 0, exposure_prior_mass < 1)
  c(rep((1 - exposure_prior_mass) / m, m), rep(exposure_prior_mass / K, K))
}

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

