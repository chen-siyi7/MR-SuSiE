generate <- function(seed, theta, scenario, prop_invalid, n = N, keep_raw = FALSE) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    m <- 20L
    K <- 4L
    omega <- matrix(runif(m * K, 0, 0.22), m, K)
    eta <- xi <- numeric(m)
    zero_relevance <- integer(0)
    support <- sort(sample.int(m, round(m * prop_invalid)))
    if (scenario == "S1")
        eta[support] <- rnorm(length(support), 0, 0.2)
    if (scenario == "S2")
        eta[support] <- rnorm(length(support), 0.1, 0.2)
    if (scenario == "S3") {
        eta[support] <- rnorm(length(support), 0, 0.2)
        xi[support] <- runif(length(support), 0, 0.1)
    }
    if (scenario == "S4") {
        zero_relevance <- sort(sample.int(m, 6))
        omega[zero_relevance, ] <- 0
        eta[support] <- rnorm(length(support), 0.1, 0.1)
    }
    G <- matrix(rbinom(2 * n * m, 2, 0.3), 2 * n, m)
    U <- as.numeric(G %*% xi) + rnorm(2 * n)
    residual_correlation <- (1 - RHO) * diag(K) + matrix(RHO, K, K)
    E <- matrix(rnorm(2 * n * K), 2 * n, K) %*% chol(residual_correlation)
    X <- G %*% omega + 0.25 * U + E
    Y <- as.numeric(X %*% theta + G %*% eta) + U + rnorm(2 * n)
    center <- function(x) sweep(as.matrix(x), 2, colMeans(as.matrix(x)), "-")
    Gx <- center(G[seq_len(n), , drop = FALSE])
    Gy <- center(G[n + seq_len(n), , drop = FALSE])
    Xx <- center(X[seq_len(n), , drop = FALSE])
    Yy <- as.numeric(center(Y[n + seq_len(n)]))
    xx <- colSums(Gx^2)
    yy <- colSums(Gy^2)
    BX <- crossprod(Gx, Xx)/xx
    SX <- sqrt(pmax(outer(rep(1, m), colSums(Xx^2)) - BX^2 * xx, 0)/(n - 2)/xx)
    by <- as.numeric(crossprod(Gy, Yy))/yy
    sy <- sqrt(pmax(sum(Yy^2) - by^2 * yy, 0)/(n - 2)/yy)
    corr <- diag(0.42, m)
    dat <- list(m = m, num_exposures = K, n1 = n, n2 = n, betaZX = BX, seZX = SX, betaZY = by, se_betaZY = sy,
        corr_in = corr, corr_ref = corr, gamma_hat = BX, sigma_exp = cor(Xx), overlap = 0)
    truth <- list(theta = theta, omega = omega, eta = eta, xi = xi, a = eta + xi, true_invalid = which(eta +
        xi != 0), zero_relevance = zero_relevance, gamma_total = omega + 0.25 * xi)
    list(data = dat, truth = truth, seed = seed, raw = if (keep_raw) list(Gx = Gx, Gy = Gy, Xx = Xx,
        Yy = Yy) else NULL)
}

