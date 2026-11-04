make_data <- function(rep) {
    set.seed(master + rep)
    Gamma <- 0.05 * matrix(rnorm(M * K), M, K) %*% chol(Sg)
    U <- sort(sample.int(M, 8))
    C <- sort(c(sample(U, 4), sample(setdiff(1:M, U), 4)))
    u <- c <- numeric(M)
    u[U] <- runif(8, 0.5, 1.5)
    u[U] <- u[U]/sqrt(mean(u[U]^2))
    gc <- sweep(Gamma[C, , drop = FALSE], 2, colMeans(Gamma[C, , drop = FALSE]), "-")
    c[C] <- as.numeric(gc %*% rep(1/sqrt(K), K))
    c[C] <- c[C]/sqrt(mean(c[C]^2))
    mix <- sqrt(0.5) * u + sqrt(0.5) * c
    I <- which(mix != 0)
    scale <- sqrt(0.0016/mean(mix[I]^2))
    u <- sqrt(0.5) * u * scale
    c <- sqrt(0.5) * c * scale
    direct <- u + c
    ref <- matrix(rnorm(config$n_ref * M), config$n_ref, M) %*% chol(Rtrue)
    Rref <- cor(ref)
    BX <- Rtrue %*% Gamma + Ltrue %*% matrix(rnorm(M * K), M, K)/sqrt(config$n_exposure)
    mu <- as.numeric(Rtrue %*% (Gamma %*% config$theta + direct))
    bys <- mu + as.numeric(Ltrue %*% rnorm(M))/sqrt(config$n_select)
    bye <- mu + as.numeric(Ltrue %*% rnorm(M))/sqrt(config$n_estimate)
    dat <- list(m = M, num_exposures = K, betaZX = BX, seZX = matrix(1/sqrt(config$n_exposure), M, K),
        betaZY = bye, se_betaZY = rep(1/sqrt(config$n_estimate), M), corr_in = Rtrue, corr_ref = Rref,
        sigma_exp = diag(K), overlap = 0)
    list(dat = dat, select = bys, Gamma = Gamma, direct = direct, U = U, C = C, uhp = u, chp = c, truth_invalid = which(direct !=
        0), mu = mu)
}

