cmp_mvmr_egger_ld <- function(dat, snps = seq_len(dat$m), which_R = "in") {
  K <- dat$num_exposures
  na <- list(beta = rep(NA, K), se = rep(NA, K), p = rep(NA, K), q = rep(NA, K))
  if (length(snps) < K + 2) return(na)
  Rfull <- if (which_R == "in") dat$corr_in else dat$corr_ref
  R <- cov2cor(Rfull)[snps, snps, drop = FALSE]
  Bx <- cbind(1, dat$betaZX[snps, , drop = FALSE])
  by <- dat$betaZY[snps]; sy <- dat$se_betaZY[snps]
  DY <- diag(sy)
  Oi <- tryCatch(solve(DY %*% R %*% DY), error = function(e) MASS_ginv(DY %*% R %*% DY))
  A <- t(Bx) %*% Oi %*% Bx
  b <- as.numeric(MASS_ginv(A) %*% (t(Bx) %*% Oi %*% by))
  resid <- as.numeric(by - Bx %*% b)
  phi <- max(1, as.numeric(t(resid) %*% Oi %*% resid) / max(1, length(snps) - K - 1))
  V <- MASS_ginv(A) * phi
  beta <- b[-1]; se <- sqrt(pmax(diag(V)[-1], 0))
  p <- 2 * pnorm(-abs(beta / se))
  list(beta = beta, se = se, p = p, q = p.adjust(p, "fdr"), intercept = b[1])
}

cmp_mv_iwas <- function(dat, snps = seq_len(dat$m), which_R = "in", sigma2 = NULL) {
  K <- dat$num_exposures
  na <- list(beta = rep(NA, K), se = rep(NA, K), p = rep(NA, K), q = rep(NA, K))
  if (length(snps) < K + 1) return(na)
  Rfull <- if (which_R == "in") dat$corr_in else dat$corr_ref
  R <- cov2cor(Rfull)[snps, snps, drop = FALSE]
  by <- dat$betaZY[snps]; sy <- dat$se_betaZY[snps]

  W <- tryCatch(MASS_ginv(R) %*% dat$betaZX[snps, , drop = FALSE],
                error = function(e) NULL)
  if (is.null(W)) return(na)

  A <- t(W) %*% R %*% W
  Ai <- MASS_ginv(A)
  beta <- as.numeric(Ai %*% (t(W) %*% by))

  resid <- as.numeric(by - R %*% W %*% beta)
  s2 <- if (!is.null(sigma2)) sigma2 else
    max(sum(resid^2 / pmax(sy^2, 1e-300)) / max(1, length(snps) - K), 1) * mean(sy^2)
  V <- Ai * s2
  se <- sqrt(pmax(diag(V), 0)); p <- 2 * pnorm(-abs(beta / se))
  list(beta = beta, se = se, p = p, q = p.adjust(p, "fdr"))
}

