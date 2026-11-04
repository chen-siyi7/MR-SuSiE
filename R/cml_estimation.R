nc4 <- function(x) any(grepl("did not converge|not converg|non.?converg", x, ignore.case = TRUE))

capture4 <- function(fun) {
    w <- character()
    err <- ""
    start <- proc.time()[3]
    ans <- tryCatch(withCallingHandlers(fun(), warning = function(z) {
        w <<- c(w, conditionMessage(z))
        invokeRestart("muffleWarning")
    }), error = function(z) {
        err <<- conditionMessage(z)
        NULL
    })
    list(result = ans, error = err, warnings = unique(w), elapsed = unname(proc.time()[3] - start))
}

native_models4 <- function(cml) {
    if (!isTRUE(cml$completed) || !length(cml$screened) || is.null(cml$step3))
        return(list(models = list(), top = integer(), clusters = list()))
    al <- cml$step3$alpha
    sc <- as.integer(cml$screened)
    cl <- lapply(seq_len(nrow(al)), function(i) which(al[i, ] > 1/length(sc)))
    rows <- which(lengths(cl) > 0)
    cl <- cl[rows]
    if (!length(cl))
        return(list(models = list(), top = integer(), clusters = list()))
    models <- list(integer())
    for (ids in cl) {
        nxt <- list()
        for (old in models) for (k in ids) {
            v <- sort(unique(c(old, sc[k])))
            nxt[[paste(v, collapse = ",")]] <- v
        }
        models <- unname(nxt)
    }
    top <- sort(unique(vapply(seq_along(cl), function(i) {
        ids <- cl[[i]]
        m <- max(al[rows[i], ids])
        min(sc[ids[al[rows[i], ids] == m]])
    }, integer(1))))
    stopifnot(any(vapply(models, function(x) identical(x, top), logical(1))))
    list(models = models, top = top, clusters = lapply(cl, function(ids) sc[ids]))
}

model_data4 <- function(dat, ids) {
    pval <- 2 * pnorm(-abs(dat$betaZX[, ids, drop = FALSE]/dat$seZX[, ids, drop = FALSE]))
    iv <- which(rowSums(pval < 5e-08) > 0)
    p <- length(ids)
    if (length(iv) < max(4L, p + 1L))
        stop("insufficient_model_instruments")
    bx <- dat$betaZX[iv, ids, drop = FALSE]
    sx <- dat$seZX[iv, ids, drop = FALSE]
    by <- matrix(dat$betaZY[iv], ncol = 1)
    sy <- dat$se_betaZY[iv]
    rho <- diag(p + 1L)
    rho[seq_len(p), seq_len(p)] <- dat$sigma_exp[ids, ids, drop = FALSE]
    si <- MVMRcML::invcov_mvmr(se_bx = sx, se_by = sy, rho_mat = rho)
    list(bx = bx, sx = sx, by = by, sy = sy, rho = rho, si = si, iv = iv, n = min(dat$n1, dat$n2), K = 0:(length(iv) -
        p - 1L))
}

bic4 <- function(z) {
    fit <- MVMRcML::MVmr_cML(b_exp = z$bx, b_out = z$by, se_bx = z$sx, Sig_inv_l = z$si, n = z$n, K_vec = z$K,
        random_start = 1L, maxit = 100L, thres = 1e-04)
    if (!identical(as.integer(fit$Converge), 0L))
        stop("native_BIC_nonconvergence")
    se <- as.numeric(MVMRcML::MVcML_SdTheta(b_exp = z$bx, b_out = z$by, Sig_inv_l = z$si, theta = fit$BIC_theta,
        zero_ind = setdiff(seq_len(nrow(z$bx)), as.integer(fit$BIC_invalid))))
    list(beta = as.numeric(fit$BIC_theta), se = se, raw = fit)
}

