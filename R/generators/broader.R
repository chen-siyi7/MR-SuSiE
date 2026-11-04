ld_matrix <- function(m, kind, rho) {
    if (kind == "independent")
        return(diag(m))
    if (kind == "ar1")
        return(rho^abs(outer(seq_len(m), seq_len(m), "-")))
    if (kind == "block") {
        block <- (seq_len(m) - 1L)%/%10L
        rr <- rho * outer(block, block, "==")
        diag(rr) <- 1
        return(rr)
    }
    stop("Unknown LD kind")
}

unit_rms <- function(x) {
    z <- which(x != 0)
    if (length(z))
        x[z] <- x[z]/sqrt(mean(x[z]^2))
    x
}

generate_data <- function(g, rep) {
    m <- g$m
    K <- 6L
    theta <- if (g$causal_configuration == "null")
        rep(0, K)
    else c(0.2, 0.1, 0.3, 0.4, 0, 0)
    seed <- as.integer(g$seed_base + rep)
    set.seed(seed)
    if (g$reused_central_data) {
        old <- readRDS(file.path(here, "inputs", "reused_central", sprintf("rep_%03d.rds", rep)))
        d <- old$generated
        d$theta <- theta
        d$seed <- old$seed
        d$shared_factor <- numeric(m)
        d$origin <- "previous independently validated central reconstruction"
        return(d)
    }
    Rtrue <- ld_matrix(m, g$ld_kind, g$ld_rho)
    Lt <- t(chol(Rtrue))
    Sg <- matrix(g$genetic_rho, K, K)
    diag(Sg) <- 1
    Gamma <- g$gamma_scale * matrix(rnorm(m * K), m, K) %*% chol(Sg)
    qn <- max(2L, as.integer(round(m * g$prop_invalid)))
    nr <- as.integer(round(qn * g$support_overlap))
    U <- sort(sample.int(m, qn))
    C <- sort(c(if (nr) sample(U, nr) else integer(), if (qn > nr) sample(setdiff(1:m, U), qn - nr) else integer()))
    u <- c <- shared <- numeric(m)
    if (g$mechanism %in% c("balanced_uhp", "unbalanced_uhp", "combined_balanced", "combined_unbalanced")) {
        u[U] <- runif(qn, 0.5, 1.5)
        if (g$mechanism %in% c("balanced_uhp", "combined_balanced"))
            u[U] <- u[U] * sample(c(-1, 1), qn, TRUE)
        u <- unit_rms(u)
    }
    if (g$mechanism %in% c("sparse_chp", "combined_balanced", "combined_unbalanced")) {
        gc <- sweep(Gamma[C, , drop = FALSE], 2, colMeans(Gamma[C, , drop = FALSE]), "-")
        c[C] <- as.numeric(gc %*% rep(1/sqrt(K), K))
        c <- unit_rms(c)
    }
    if (g$mechanism == "shared_factor_chp") {
        shared[C] <- as.numeric(scale(rnorm(qn), center = TRUE, scale = FALSE))
        shared <- unit_rms(shared)
        Gamma <- Gamma + outer(shared, rep(0.25 * g$gamma_scale, K))
        c <- shared
    }
    if (g$mechanism == "dense_aligned_chp")
        c <- unit_rms(as.numeric(Gamma %*% rep(1/sqrt(K), K)))
    if (g$mechanism %in% c("combined_balanced", "combined_unbalanced")) {
        u <- sqrt(1 - g$chp_fraction) * u
        c <- sqrt(g$chp_fraction) * c
    }
    direct <- u + c
    I <- which(direct != 0)
    if (length(I)) {
        sc <- sqrt(g$pleio_m2/mean(direct[I]^2))
        u <- sc * u
        c <- sc * c
        direct <- u + c
    }
    if (is.finite(g$n_ref)) {
        rr <- if (is.finite(g$ref_rho))
            ld_matrix(m, g$ld_kind, g$ref_rho)
        else Rtrue
        ref <- matrix(rnorm(g$n_ref * m), g$n_ref, m) %*% chol(rr)
        Rref_raw <- cor(ref)
    }
    else Rref_raw <- Rtrue
    Rref <- (Rref_raw + g$ridge * diag(m))/(1 + g$ridge)
    BX <- Rtrue %*% Gamma + Lt %*% matrix(rnorm(m * K), m, K)/sqrt(g$n_exposure)
    mu <- as.numeric(Rtrue %*% (Gamma %*% theta + direct))
    bys <- mu + as.numeric(Lt %*% rnorm(m))/sqrt(g$n_select)
    bye <- mu + as.numeric(Lt %*% rnorm(m))/sqrt(g$n_estimate)
    dat <- list(m = m, num_exposures = K, betaZX = BX, seZX = matrix(1/sqrt(g$n_exposure), m, K), betaZY = bye,
        se_betaZY = rep(1/sqrt(g$n_estimate), m), corr_in = Rtrue, corr_ref = Rref, sigma_exp = diag(K),
        overlap = 0)
    list(dat = dat, select = bys, Gamma = Gamma, direct = direct, U = U, C = C, uhp = u, chp = c, truth_invalid = I,
        mu = mu, theta = theta, seed = seed, shared_factor = shared, reference_raw = Rref_raw, origin = "new broad-grid reconstruction")
}

