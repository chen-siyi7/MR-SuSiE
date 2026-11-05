mv_MR_susie_sel_sequential <-
function (corr, betaZY, se_betaZY, gamma_hat, se_gamma = NULL, n2, L = 10, exposure_prior_mass = 0.5, 
    coverage = 0.95, me_correct = FALSE, sigma_exp = NULL, sigma_res = NULL, n1 = NULL, purity = 0.5, 
    ambiguous = c("unresolved", "invalid", "exposure")) 
{
    ambiguous <- match.arg(ambiguous)
    m <- nrow(corr)
    K <- ncol(gamma_hat)
    stopifnot(is.finite(n2), n2 > 1)
    ZtY <- .suff_zty(corr, betaZY)
    yty <- .suff_yty(corr, betaZY, se_betaZY, n2)
    L <- max(1, min(L, m + K - 1))
    s1 <- susieR::susie_suff_stat(XtX = corr, Xty = ZtY, yty = yty, n = n2, L = max(1, min(L, m - 1)), 
        prior_weights = rep(1/m, m), estimate_residual_variance = FALSE, estimate_prior_variance = TRUE, 
        check_input = TRUE)
    sets1 <- .credible_sets(s1, corr, coverage, purity)
    inv1 <- sort(unique(unlist(sets1$pure)))
    valid_snps <- setdiff(seq_len(m), inv1)
    v <- valid_snps
    if (length(v) < K) {
        return(list(valid = valid_snps, invalid = inv1, invalid_cs = sets1$pure, exposures = integer(0), 
            exposure_cs = list(), unresolved = list(), n_unresolved = NA_integer_, n_unresolved_pure = NA_integer_, 
            n_cs = length(sets1$pure), n_cs_unfiltered = length(sets1$all), pip = s1$pip, pip_exposures = rep(NA_real_, 
                K), L_used = L, n_used = n2, order = "sequential", converged = isTRUE(s1$converged)))
    }
    Cv <- corr[v, v, drop = FALSE]
    Gv <- gamma_hat[v, , drop = FALSE]
    XtX2 <- t(Gv) %*% Cv %*% Gv
    if (me_correct) 
        warning("selection-stage measurement-error subtraction is disabled; ", "using the exact exposure cross-product", 
            call. = FALSE)
    Xty2 <- as.numeric(t(Gv) %*% ZtY[v])
    s2 <- tryCatch(susieR::susie_suff_stat(XtX = XtX2, Xty = Xty2, yty = yty, n = n2, L = max(1, min(L, 
        K - 1)), prior_weights = rep(1/K, K), estimate_residual_variance = FALSE, estimate_prior_variance = TRUE, 
        check_input = TRUE), error = function(e) NULL)
    base <- list(valid = valid_snps, invalid = inv1, invalid_cs = sets1$pure, unresolved = list(), n_unresolved = NA_integer_, 
        n_unresolved_pure = NA_integer_, L_used = L, n_used = n2, order = "sequential", converged = isTRUE(s1$converged) && 
            (!is.null(s2) && isTRUE(s2$converged)))
    if (is.null(s2)) 
        return(c(base, list(exposures = integer(0), exposure_cs = list(), n_cs = length(sets1$pure), 
            n_cs_unfiltered = length(sets1$all), pip = s1$pip, pip_exposures = rep(NA_real_, K))))
    sets2 <- .credible_sets(s2, XtX2, coverage, purity)
    c(base, list(exposures = sort(unique(unlist(sets2$pure))), exposure_cs = sets2$pure, n_cs = length(sets1$pure) + 
        length(sets2$pure), n_cs_unfiltered = length(sets1$all) + length(sets2$all), pip = s1$pip, pip_exposures = s2$pip))
}
