test_that("adjusted estimator returns a complete diagnostic object", {
  R <- diag(4)
  bx <- matrix(c(0.20, 0.02, 0.16, -0.01,
                 0.01, 0.18, 0.03, 0.15), 4, 2, byrow = TRUE)
  dat <- list(
    num_exposures = 2L, betaZX = bx,
    betaZY = as.numeric(bx %*% c(0.4, -0.2)),
    se_betaZY = rep(0.02, 4), seZX = matrix(0.01, 4, 2),
    corr_in = R, corr_ref = R, sigma_exp = diag(2), overlap = 0
  )
  fit <- mvmr_ivw_ld_adjusted(dat, me = FALSE)
  expect_true(fit$identified)
  expect_length(fit$beta, 2)
  expect_equal(fit$beta, c(0.4, -0.2), tolerance = 1e-6)
})
