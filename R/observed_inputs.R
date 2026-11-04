r2_assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

r2_ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
    s <- svd(X)
    keep <- s$d > max(s$d) * tol
    if (!any(keep))
        return(matrix(0, ncol(X), nrow(X)))
    s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE])/s$d[keep])
}

r2_observed <- function(dat, selection_outcome = NULL, selection_se = NULL, n_select = NULL, selector_family = "gaussian",
    outcome_reuse = TRUE, gamma_mode = c("reference_joint", "supplied"), n_reference = NULL) {
    gamma_mode <- match.arg(gamma_mode)
    K <- dat$num_exposures
    m <- dat$m
    r2_assert(identical(dim(dat$betaZX), c(as.integer(m), as.integer(K))), "Exposure dimensions differ.")
    r2_assert(identical(dim(dat$seZX), dim(dat$betaZX)), "Exposure SE dimensions differ.")
    r2_assert(length(dat$betaZY) == m && length(dat$se_betaZY) == m, "Outcome dimensions differ.")
    R <- dat$corr_ref
    r2_assert(identical(dim(R), c(as.integer(m), as.integer(m))), "Reference LD dimensions differ.")
    allowed <- c("m", "num_exposures", "betaZX", "seZX", "betaZY", "se_betaZY", "sigma_exp", "overlap",
        "n1", "n2", "rsid", "traits", "reference_dosage_sd", "overlap_status")
    fit_dat <- dat[intersect(allowed, names(dat))]
    fit_dat$corr_ref <- R
    fit_dat$corr_in <- R
    bsel <- if (is.null(selection_outcome))
        dat$betaZY
    else selection_outcome
    ssel <- if (is.null(selection_se))
        dat$se_betaZY
    else selection_se
    r2_assert(length(bsel) == m && length(ssel) == m, "Selection summary dimensions differ.")
    gamma <- if (gamma_mode == "supplied")
        dat$gamma_hat
    else r2_ginv(R) %*% dat$betaZX
    r2_assert(identical(dim(gamma), dim(dat$betaZX)), "Selector dictionary dimensions differ.")
    list(dat = fit_dat, selection = list(beta = bsel, se = ssel, gamma_hat = gamma, n_eff = n_select,
        family = selector_family, outcome_reuse = outcome_reuse), reference = list(n = n_reference, diagonal = diag(R),
        scale = if (all(abs(diag(R) - 1) < 1e-10)) "correlation" else "genotype_covariance"), truth_use = "Observed inputs only; population LD is excluded or replaced by the supplied reference matrix.")
}

