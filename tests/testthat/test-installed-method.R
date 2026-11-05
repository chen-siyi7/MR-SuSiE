test_that("installed scalar and joint fits reproduce saved numerical examples", {
  previous <- getwd()
  working <- tempfile("mr_susie_installed_")
  dir.create(working)
  setwd(working)
  on.exit(setwd(previous), add=TRUE)
  model <- make_mr_susie()
  expect_true(startsWith(model$root, system.file(package="MRSuSiE")))
  for (name in c("scalar_strong_0001", "joint_strong_0001")) {
    x <- jsonlite::read_json(system.file("examples", paste0(name,".json"), package="MRSuSiE"), simplifyVector=TRUE)
    m <- nrow(x$B)
    k <- ncol(x$B)
    fit <- model$fit(x$B, matrix(x$sx,m,k), x$y, rep(x$sy,m), x$R, x$Sigma,
      x$selection_y, rep(x$sy,m), 10000, rep("Outcome SD per exposure SD",k))
    expected <- read.csv(system.file("examples", paste0(name,"_expected.csv"), package="MRSuSiE"))
    expect_equal(as.matrix(data.frame(point=fit$point, fit$score_inference)), as.matrix(expected), tolerance=1e-9, ignore_attr=TRUE)
    expect_equal((fit$score_inference$lower+fit$score_inference$upper)/2, fit$score_inference$center, tolerance=1e-12)
    expect_identical(fit$method, if (k==1L) "UVMR-SuSiE" else "MVMR-SuSiE")
    expect_length(list.files(model$work_dir, recursive=TRUE), 0L)
  }
})

test_that("support-union example is available through the installed package", {
  x <- jsonlite::read_json(system.file("examples", "support_union_example.json", package="MRSuSiE"), simplifyVector=TRUE)
  fit <- do.call(union_ar, x)
  expect_true(is.list(fit$sets))
  expect_true(is.finite(fit$joint_null_p))
})
