## comparators.R -- UV and MV MR comparator methods for the head-to-head with
## (MV)MR-SuSiE. Three tiers:
##   (a) from scratch, base R, tested here: IVW (UV/MV), weighted median (UV/MV);
##   (b) glmnet-based MR-Lasso (UV/MV), if glmnet is installed;
##   (c) wrappers to the canonical MendelianRandomization package
##       (mr_median, mr_lasso, mr_mvivw, mr_mvmedian, mr_mvlasso), which are the
##       reference implementations reviewers will expect. Reimplementing the
##       lasso-based ones from memory risks silent deviation, so prefer (c) when
##       the package is installed; (a)/(b) are dependency-light fallbacks.
##
## All MV methods take dat (from gen_mv) and a valid-instrument index set, and
## return per-exposure list(beta, se, p, q). Exposure selection = which(q<0.05).

## ============ from-scratch, tested ============

if (!exists("MASS_ginv", mode = "function")) {
  MASS_ginv <- function(M, tol = sqrt(.Machine$double.eps)) {
    s <- svd(M)
    if (!length(s$d)) return(matrix(0, ncol(M), nrow(M)))
    keep <- s$d > max(s$d) * tol
    if (!any(keep)) return(matrix(0, ncol(M), nrow(M)))
    s$v[, keep, drop = FALSE] %*%
      (t(s$u[, keep, drop = FALSE]) / s$d[keep])
  }
}

## LD-aware MVMR after invalid-variant selection.
##
## Marginal summary associations obey
##   betaY = Bx theta + R[, I] alpha_I + error,
## where I is the selected support of direct (non-exposure-mediated) effects.
## Merely dropping rows I leaves R[V, I] alpha_I in the retained marginal
## associations whenever instruments are in LD. We instead keep all rows and
## include R[, I] as nuisance columns. With the correct I this is the
## summary-statistic analogue of conditioning on the invalid variants.
##
## The SNP--exposure estimates are noisy regressors. Merely adding their error
## variance to Omega changes the GLS weights but leaves attenuation in the
## normal equations because
##   E[Bhat' W (bY - Bhat theta)] = -E[E' W E] theta.
## We therefore use the corrected score
##   {H' W H - blockdiag(C(W),0)} zeta = H' W bY,
## where C_kl(W) = tr{W Cov(E_k,E_l)}. Under the separable summary-error
## approximation
##   Cov(E_k,E_l) = sigma_exp[k,l] diag(sX_k) R diag(sX_l),
## this is available from summary standard errors and LD. The same approximation
## gives the FULL cross-variant residual covariance, rather than only its
## diagonal. These corrections assume non-overlapping exposure and outcome GWAS;
## sample overlap needs the corresponding exposure--outcome covariance terms.
mvmr_ivw_ld_adjusted <- function(dat, invalid_snps = integer(0),
                                 which_R = c("in", "ref"), me = TRUE,
                                 sigma_exp = NULL) {
  which_R <- match.arg(which_R)
  K <- dat$num_exposures
  if (!is.null(dat$overlap) && is.finite(dat$overlap) && dat$overlap > 0)
    return(list(beta = rep(NA, K), se = rep(NA, K), p = rep(NA, K),
                q = rep(NA, K), identified = FALSE,
                overlap_supported = FALSE))
  if (is.null(sigma_exp)) sigma_exp <- dat$sigma_exp
  if (is.null(sigma_exp) && K > 1L)
    return(list(beta = rep(NA, K), se = rep(NA, K), p = rep(NA, K),
                q = rep(NA, K), identified = FALSE,
                overlap_supported = TRUE,
                exposure_error_correlation_available = FALSE))
  if (is.null(sigma_exp)) sigma_exp <- matrix(1, 1, 1)
  Rfull <- if (which_R == "in") dat$corr_in else dat$corr_ref
  gsd <- sqrt(pmax(diag(Rfull), .Machine$double.eps))
  R <- cov2cor(Rfull)
  I <- sort(unique(as.integer(invalid_snps)))
  I <- I[I >= 1L & I <= nrow(R)]
  Bx <- sweep(dat$betaZX, 1, gsd, "*")
  by <- dat$betaZY * gsd
  sy <- dat$se_betaZY * gsd
  sX <- sweep(dat$seZX, 1, gsd, "*")
  H <- cbind(Bx, if (length(I)) R[, I, drop = FALSE] else NULL)
  Hscale <- sqrt(colSums(H^2))
  Hscale[!is.finite(Hscale) | Hscale == 0] <- 1
  Hs <- sweep(H, 2, Hscale, "/")
  rank_H <- qr(Hs, tol = 1e-8)$rank
  condition_H <- tryCatch(kappa(Hs), error = function(e) Inf)
  DY <- diag(sy)
  DRD <- DY %*% R %*% DY
  solveW <- function(M) tryCatch(solve(M), error = function(e) MASS_ginv(M))
  Sg <- sigma_exp
  stopifnot(nrow(Sg) == K, ncol(Sg) == K)
  Sg <- (Sg + t(Sg)) / 2

  ## Full covariance of E %*% theta across variants:
  ## R * {(S_X diag(theta)) Sigma_exp (S_X diag(theta))'}.
  exposure_error_cov <- function(b) {
    Vtheta <- sweep(sX, 2, b, "*")
    R * (Vtheta %*% Sg %*% t(Vtheta))
  }

  ## Expected noisy-regressor contribution E[E' W E]. The elementwise form
  ## avoids K^2 dense m x m products:
  ## C = Sigma_exp * {S_X' (W * R) S_X}.
  error_information <- function(Oi) {
    C <- Sg * (t(sX) %*% (Oi * R) %*% sX)
    (C + t(C)) / 2
  }

  corrected_fit <- function(Oi, correct_me = TRUE) {
    A <- t(H) %*% Oi %*% H
    C <- matrix(0, K, K)
    if (correct_me) {
      C <- error_information(Oi)
      A[seq_len(K), seq_len(K)] <-
        A[seq_len(K), seq_len(K), drop = FALSE] - C
    }
    rhs <- as.numeric(t(H) %*% Oi %*% by)
    coef <- as.numeric(MASS_ginv(A) %*% rhs)
    list(coef = coef, A = A, C = C)
  }

  Oi <- solveW(DRD)
  fit <- corrected_fit(Oi, correct_me = FALSE)
  coef <- fit$coef
  beta <- coef[seq_len(K)]
  if (me) {
    for (it in seq_len(8)) {
      Omega <- DRD + exposure_error_cov(beta)
      Oi <- solveW(Omega)
      fit_new <- corrected_fit(Oi, correct_me = TRUE)
      coef_new <- fit_new$coef
      beta_new <- coef_new[seq_len(K)]
      converged <- all(is.finite(beta_new)) &&
        max(abs(beta_new - beta)) < 1e-10 * max(1, max(abs(beta_new)))
      fit <- fit_new
      coef <- coef_new
      beta <- beta_new
      if (converged) break
    }
  } else {
    Omega <- DRD
  }

  A <- fit$A
  ## Identification of theta does not require the nuisance coefficients to be
  ## unique.  Residualize the exposure columns on the selected LD directions
  ## and report their rank separately from the (possibly deficient) rank of H.
  Bx_res <- Bx
  if (length(I)) {
    Z <- R[, I, drop = FALSE]
    Bx_res <- Bx - Z %*% MASS_ginv(t(Z) %*% Oi %*% Z) %*%
      t(Z) %*% Oi %*% Bx
  }
  Bx_res_scale <- sqrt(colSums(Bx_res^2))
  Bx_res_scale[!is.finite(Bx_res_scale) | Bx_res_scale == 0] <- 1
  rank_causal <- qr(sweep(Bx_res, 2, Bx_res_scale, "/"),
                    tol = 1e-8)$rank
  if (rank_causal < K)
    return(list(beta = rep(NA, K), se = rep(NA, K), p = rep(NA, K),
                q = rep(NA, K), identified = FALSE, rank_H = rank_H,
                rank_causal = rank_causal, ncol_H = ncol(H),
                condition_H = condition_H, overlap_supported = TRUE))

  ## The causal block must remain positive after partialling out the selected
  ## nuisance directions. A non-positive block is a weak/nonidentified corrected
  ## score, not a matrix to rescue silently with a generalized inverse.
  Axx <- A[seq_len(K), seq_len(K), drop = FALSE]
  if (length(I)) {
    Axa <- A[seq_len(K), K + seq_along(I), drop = FALSE]
    Aaa <- A[K + seq_along(I), K + seq_along(I), drop = FALSE]
    Axx <- Axx - Axa %*% MASS_ginv(Aaa) %*% t(Axa)
  }
  ev <- eigen((Axx + t(Axx)) / 2, symmetric = TRUE, only.values = TRUE)$values
  strength_ok <- all(is.finite(ev)) && min(ev) > 1e-10 * max(1, max(ev))
  if (!strength_ok)
    return(list(beta = rep(NA, K), se = rep(NA, K), p = rep(NA, K),
                q = rep(NA, K), identified = FALSE, rank_H = rank_H,
                rank_causal = rank_causal, ncol_H = ncol(H),
                condition_H = condition_H,
                min_corrected_eigen = min(ev), overlap_supported = TRUE))

  resid <- as.numeric(by - H %*% coef)
  ## Redundant selected LD columns do not consume additional residual degrees
  ## of freedom; use the effective design rank rather than its raw column count.
  phi <- max(1, as.numeric(t(resid) %*% Oi %*% resid) /
               max(1, nrow(H) - rank_H))
  ## Working inverse-information covariance. This is not the sandwich covariance
  ## of the nonlinear corrected score: the full score Jacobian and score variance
  ## are not yet implemented. It also omits selection randomness, LD-reference
  ## uncertainty, and higher-order weak-IV terms.
  V <- MASS_ginv(A) * phi
  Axx_inv <- MASS_ginv(Axx)
  conditional_information <- 1 / pmax(diag(Axx_inv), .Machine$double.eps)
  names(conditional_information) <- colnames(dat$betaZX)
  corrected_condition <- tryCatch(kappa(Axx), error = function(e) Inf)
  se <- sqrt(pmax(diag(V)[seq_len(K)], 0))
  p <- 2 * pnorm(-abs(beta / se))
  list(beta = beta, se = se, p = p, q = p.adjust(p, "fdr"),
       direct = if (length(I)) coef[K + seq_along(I)] else numeric(0),
       invalid_adjusted = I, identified = TRUE, rank_H = rank_H,
       rank_causal = rank_causal, ncol_H = ncol(H),
       condition_H = condition_H,
       min_corrected_eigen = min(ev), me_corrected_score = me,
       corrected_eigenvalues = ev,
       corrected_condition = corrected_condition,
       conditional_information = conditional_information,
       cross_variant_me = me, overlap_supported = TRUE,
       exposure_error_correlation_available = TRUE)
}

