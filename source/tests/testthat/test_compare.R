set.seed(123)

LLarr <- example_loglik_array()
LLarr2 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 0.5), dim = dim(LLarr))
LLarr3 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 1), dim = dim(LLarr))
w1 <- suppressWarnings(waic(LLarr))
w2 <- suppressWarnings(waic(LLarr2))

test_that("loo_compare throws appropriate errors", {
  w3 <- suppressWarnings(waic(LLarr[,, -1]))
  w4 <- suppressWarnings(waic(LLarr[,, -(1:2)]))

  expect_error(loo_compare(2, 3), "must be a list if not a 'loo' object")
  expect_error(
    loo_compare(w1, w2, x = list(w1, w2)),
    "If 'x' is a list then '...' should not be specified"
  )
  expect_error(loo_compare(w1, list(1, 2, 3)), "class 'loo'")
  expect_error(loo_compare(w1), "requires at least two models")
  expect_error(loo_compare(x = list(w1)), "requires at least two models")
  expect_error(
    loo_compare(w1, w3),
    "All models must have the same number of observations, but models have inconsistent observation counts: 'model1' (32), 'model2' (31)",
    fixed = TRUE
  )
  expect_error(
    loo_compare(w1, w2, w3),
    "All models must have the same number of observations, but models have inconsistent observation counts: 'model1' (32), 'model2' (32), 'model3' (31)",
    fixed = TRUE
  )
  expect_error(
    loo_compare(x = list("Model A" = w1, "Model B" = w2, "Model C" = w3)),
    "All models must have the same number of observations, but models have inconsistent observation counts: 'Model A' (32), 'Model B' (32), 'Model C' (31)",
    fixed = TRUE
  )
})

test_that("loo_compare throws appropriate warnings", {
  w3 <- w1
  w4 <- w2
  class(w3) <- class(w4) <- c("kfold", "loo")
  attr(w3, "K") <- 2
  attr(w4, "K") <- 3
  expect_warning(
    loo_compare(w3, w4),
    "Not all kfold objects have the same K value"
  )

  class(w4) <- c("psis_loo", "loo")
  attr(w4, "K") <- NULL
  expect_warning(loo_compare(w3, w4), "Comparing LOO-CV to K-fold-CV")

  w3 <- w1
  w4 <- w2
  attr(w3, "yhash") <- "a"
  attr(w4, "yhash") <- "b"
  expect_warning(loo_compare(w3, w4), "Not all models have the same y variable")

  set.seed(123)
  w_list <- lapply(1:25, function(x) {
    suppressWarnings(waic(LLarr + rnorm(1, 0, 0.1)))
  })
  expect_warning(
    loo_compare(w_list),
    "Difference in performance potentially due to chance. See McLatchie and Vehtari (2023) for details.",
    fixed = TRUE
  )

  w_list_short <- lapply(1:4, function(x) {
    suppressWarnings(waic(LLarr + rnorm(1, 0, 0.1)))
  })
  expect_no_warning(loo_compare(w_list_short))
})


comp_colnames <- c(
  "model",
  "elpd_diff",
  "se_diff",
  "p_worse",
  "diag_diff",
  "diag_elpd",
  "elpd_waic",
  "se_elpd_waic",
  "p_waic",
  "se_p_waic",
  "waic",
  "se_waic"
)

test_that("loo_compare returns expected results (2 models)", {
  comp1 <- loo_compare(w1, w1)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "data.frame")
  expect_equal(colnames(comp1), comp_colnames)
  expect_equal(comp1$model, c("model1", "model2"))
  expect_equal(comp1$elpd_diff, c(0, 0), ignore_attr = TRUE)
  expect_equal(comp1$se_diff, c(0, 0), ignore_attr = TRUE)
  expect_equal(comp1$p_worse, c(NA_real_, NA_real_), ignore_attr = TRUE)
  expect_snapshot_value(comp1, style = "serialize")
  expect_snapshot(print(comp1))

  comp2 <- loo_compare(w1, w2)
  expect_s3_class(comp2, "compare.loo")
  expect_equal(colnames(comp2), comp_colnames)
  expect_equal(comp2$p_worse, c(NA, 1))
  expect_equal(comp2$diag_diff, c("", "N < 100"))
  expect_equal(comp2$diag_elpd, c("", ""))
  expect_snapshot_value(comp2, style = "serialize")
  expect_snapshot(print(comp2))
  expect_snapshot(print(comp2, p_worse = FALSE))
  expect_snapshot(print(comp2, simplify = FALSE))
  expect_snapshot(print(comp2, simplify = FALSE, p_worse = FALSE))

  # specifying objects via ... and via arg x gives equal results
  expect_equal(comp2, loo_compare(x = list(w1, w2)))

  # custom naming works
  comp3 <- loo_compare(x = list("A" = w2, "B" = w1))
  expect_equal(comp3$model, c("B", "A"))
})

test_that("loo_compare returns expected result (3 models)", {
  w3 <- suppressWarnings(waic(LLarr3))
  comp1 <- loo_compare(w1, w2, w3)

  expect_equal(colnames(comp1), comp_colnames)
  expect_equal(comp1$model, c("model1", "model2", "model3"))
  expect_equal(comp1$p_worse, c(NA, 1, 1))
  expect_equal(comp1$diag_diff, c("", "N < 100", "N < 100"))
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "data.frame")
  expect_snapshot_value(comp1, style = "serialize")
  expect_snapshot(print(comp1))

  # specifying objects via '...' gives equivalent results (equal
  # except rownames) to using 'x' argument
  expect_equal(comp1, loo_compare(x = list(w1, w2, w3)), ignore_attr = TRUE)
})

test_that("loo_compare with simplify=FALSE returns expected result", {
  LL <- example_loglik_array()
  loo1 <- loo(LL)
  loo2 <- loo(LL + 1)
  loo3 <- loo(LL + 2)
  comp <- loo_compare(loo1, loo2, loo3)
  expect_snapshot(print(comp, simplify = FALSE))
})

# Cluster-level joint sign-flip test ---------------------------------------

make_cluster_loos <- function(n_clusters = 6, per_cluster = 4, S = 40,
                              delta2 = NULL, delta3 = NULL, seed = 7) {
  set.seed(seed)
  N <- n_clusters * per_cluster
  ids <- sprintf("tx%03d", seq_len(N))
  cluster <- rep(sprintf("city%d", seq_len(n_clusters)), each = per_cluster)
  base <- matrix(rnorm(S * N, -2, 0.5), S, N)
  make_ll <- function(delta) {
    m <- base + matrix(rep(delta, each = S), S, N)
    colnames(m) <- ids
    m
  }
  zero <- rep(0, N)
  if (is.null(delta2)) {
    delta2 <- replace(zero, seq_len(n_clusters - 1) * per_cluster, -3)
    delta2[delta2 == 0] <- 2
  }
  if (is.null(delta3)) {
    delta3 <- delta2 / 2
  }
  loos <- lapply(list(zero, delta2, delta3), function(d) {
    suppressWarnings(loo(make_ll(d / S * 4)))
  })
  names(loos) <- c("deployed", "cand1", "cand2")
  list(loos = loos, cluster = cluster, ids = ids)
}

brute_force_cluster_test <- function(D) {
  K <- nrow(D)
  M <- ncol(D)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), K)))
  denom <- sqrt(colSums(D^2))
  safe_denom <- ifelse(denom == 0, 1, denom)
  obs <- colSums(D) / safe_denom
  obs[denom == 0] <- 0
  stat_mat <- sweep(signs %*% D, 2L, safe_denom, "/")
  stat_mat[, denom == 0] <- 0
  # floating-point ties count as "at least as large"
  tol <- 1e-10 * max(1, sqrt(K))
  thr <- obs - tol
  cand <- which(denom > 0 | colSums(D^2) > 0)
  raw <- colMeans(
    stat_mat[, cand, drop = FALSE] >=
      matrix(thr[cand], nrow(signs), length(cand), byrow = TRUE)
  )
  adj <- colMeans(
    apply(stat_mat[, cand, drop = FALSE], 1, max) >=
      matrix(thr[cand], nrow(signs), length(cand), byrow = TRUE)
  )
  names(raw) <- names(adj) <- colnames(D)[cand]
  list(stat = obs, p_raw = raw, p_adj = adj, candidates = cand)
}

test_that("cluster test matches brute-force enumeration", {
  dat <- make_cluster_loos()
  comp <- loo_compare(dat$loos[[1]], dat$loos[[2]], dat$loos[[3]],
                      clusters = dat$cluster)
  ct <- attr(comp, "cluster_test")
  ref <- brute_force_cluster_test(ct$diffs)

  expect_equal(ct$n_clusters, 6L)
  expect_equal(ct$n_sign_flips, 64L)
  expect_equal(unname(ct$stat[ct$candidates]), unname(ref$stat[ref$candidates]))
  expect_equal(unname(ct$p_raw[ct$candidates]), unname(ref$p_raw))
  expect_equal(unname(ct$p_adj[ct$candidates]), unname(ref$p_adj))
  # the reference model is not a candidate
  expect_true(is.na(ct$stat[1L]))
  expect_true(is.na(ct$p_raw[1L]))
  expect_true(is.na(ct$p_adj[1L]))
  # ties and the identity pattern are counted, so p-values are >= 1 / 2^K
  expect_true(all(ct$p_raw[ct$candidates] >= 1 / 2^6))
  expect_true(all(ct$p_adj[ct$candidates] >= 1 / 2^6))
})

test_that("cluster test is invariant to row shuffling", {
  dat <- make_cluster_loos()
  comp1 <- loo_compare(dat$loos[[1]], dat$loos[[2]], dat$loos[[3]],
                       clusters = dat$cluster)
  ct1 <- attr(comp1, "cluster_test")

  perm <- sample(length(dat$ids))
  shuffled <- lapply(dat$loos, function(x) {
    x$pointwise <- x$pointwise[perm, , drop = FALSE]
    x
  })
  comp2 <- loo_compare(shuffled[[1]], shuffled[[2]], shuffled[[3]],
                       clusters = dat$cluster[perm])
  ct2 <- attr(comp2, "cluster_test")
  expect_equal(ct1$p_raw, ct2$p_raw)
  expect_equal(ct1$p_adj, ct2$p_adj)

  # a data frame keyed by transaction ID given in reverse row order, with the
  # models also shuffled, gives the same result
  df <- data.frame(
    transaction_id = rev(dat$ids),
    city = rev(dat$cluster),
    stringsAsFactors = FALSE
  )
  comp3 <- loo_compare(shuffled[[1]], shuffled[[2]], shuffled[[3]],
                       clusters = df)
  ct3 <- attr(comp3, "cluster_test")
  expect_equal(ct1$p_raw, ct3$p_raw)
  expect_equal(ct1$p_adj, ct3$p_adj)

  # named vector keyed by transaction ID
  nv <- setNames(rev(dat$cluster), rev(dat$ids))
  comp4 <- loo_compare(shuffled[[1]], shuffled[[2]], shuffled[[3]],
                       clusters = nv)
  ct4 <- attr(comp4, "cluster_test")
  expect_equal(ct1$p_raw, ct4$p_raw)
  expect_equal(ct1$p_adj, ct4$p_adj)

  # an unnamed vector is positional (first model's row order). With row names
  # on every model the other models are still aligned by ID, so shuffling only
  # the candidate models' rows does not change the result.
  comp5 <- loo_compare(dat$loos[[1]], shuffled[[2]], shuffled[[3]],
                       clusters = dat$cluster)
  ct5 <- attr(comp5, "cluster_test")
  expect_equal(ct1$p_raw, ct5$p_raw)
  expect_equal(ct1$p_adj, ct5$p_adj)
})

test_that("cluster test counts exact ties (stable arithmetic)", {
  # example_loglik_array() gives identical pointwise differences for all
  # observations, so the observed statistic equals the maximum statistic for
  # the identity pattern and many other patterns. It must still get p >= 1/2^K.
  LL <- example_loglik_array()
  l1 <- loo(LL)
  l2 <- loo(LL + 1)
  l3 <- loo(LL + 2)
  rn <- sprintf("tx%02d", seq_len(nrow(l1$pointwise)))
  add_ids <- function(x) {
    rownames(x$pointwise) <- rn
    x
  }
  l1 <- add_ids(l1)
  l2 <- add_ids(l2)
  l3 <- add_ids(l3)
  city <- rep(c("A", "B", "C", "D"), length.out = length(rn))
  comp <- suppressWarnings(loo_compare(l1, l2, l3, clusters = city))
  ct <- attr(comp, "cluster_test")
  expect_equal(ct$n_sign_flips, 16L)
  expect_true(all(ct$p_joint[!is.na(comp$p_joint)] >= 1 / 16))
  expect_true(any(comp$p_joint == 1 / 16))
})

test_that("cluster test works for waic and kfold objects", {
  LL <- example_loglik_array()
  ids <- sprintf("tx%02d", seq_len(dim(LL)[3]))
  LL2 <- array(rnorm(prod(dim(LL)), c(LL), 0.4), dim = dim(LL))
  dimnames(LL) <- dimnames(LL2) <- list(NULL, NULL, ids)
  cl <- rep(c("a", "b", "c", "d", "e"), length.out = length(ids))

  w1 <- suppressWarnings(waic(LL))
  w2 <- suppressWarnings(waic(LL2))
  cm <- loo_compare(w1, w2, clusters = cl)
  expect_true(all(c("stat_joint", "p_joint", "p_joint_adj") %in% colnames(cm)))
  expect_equal(sum(!is.na(cm$p_joint)), 1)

  k1 <- loo(LL)
  k2 <- loo(LL2)
  class(k1) <- class(k2) <- c("kfold", "loo")
  attr(k1, "K") <- attr(k2, "K") <- 10
  cmk <- suppressWarnings(loo_compare(k1, k2, clusters = cl))
  expect_equal(cm$p_joint, cmk$p_joint)
})

test_that("cluster test handles all-zero cluster differences", {
  dat <- make_cluster_loos()
  comp <- loo_compare(dat$loos[[1]], dat$loos[[1]], dat$loos[[3]],
                      clusters = dat$cluster)
  ct <- attr(comp, "cluster_test")
  expect_identical(unname(ct$stat[2L]), 0)
  expect_identical(unname(ct$p_raw[2L]), 1)
  expect_identical(unname(ct$p_adj[2L]), 1)
})

test_that("cluster test returns raw and adjusted p-values", {
  dat <- make_cluster_loos()
  comp <- loo_compare(dat$loos[[1]], dat$loos[[2]], dat$loos[[3]],
                      clusters = dat$cluster)
  expect_true(all(c("stat_joint", "p_joint", "p_joint_adj") %in% colnames(comp)))
  ct <- attr(comp, "cluster_test")
  m <- match(comp$model, names(ct$stat))
  expect_equal(comp$stat_joint, unname(ct$stat[m]))
  expect_equal(comp$p_joint, unname(ct$p_raw[m]))
  expect_equal(comp$p_joint_adj, unname(ct$p_adj[m]))
  # adjusted p-values are never smaller than raw ones
  expect_true(all(comp$p_joint_adj >= comp$p_joint, na.rm = TRUE))
  expect_output(print(comp), "stat_joint")
})

test_that("cluster test supports choosing the reference model", {
  dat <- make_cluster_loos()
  base <- loo_compare(dat$loos, clusters = dat$cluster)
  d_cand1 <- attr(base, "cluster_test")$diffs[, "cand1"]

  comp <- loo_compare(dat$loos, clusters = dat$cluster, reference = "cand1")
  ct <- attr(comp, "cluster_test")
  expect_equal(ct$reference, "cand1")
  expect_true(is.na(ct$stat["cand1"]))
  expect_equal(
    unname(ct$stat["deployed"]),
    sum(-d_cand1) / sqrt(sum(d_cand1^2))
  )

  # the reference can also be given by index or as a loo object
  by_index <- loo_compare(dat$loos, clusters = dat$cluster, reference = 2L)
  expect_equal(attr(by_index, "cluster_test")$reference, "cand1")
  by_obj <- loo_compare(dat$loos, clusters = dat$cluster,
                        reference = dat$loos[["cand1"]])
  expect_equal(attr(by_obj, "cluster_test")$reference, "cand1")
})

test_that("cluster test does not change the old comparison output", {
  LL <- example_loglik_array()
  w1 <- suppressWarnings(waic(LL))
  w2 <- suppressWarnings(waic(LL + 1))
  old <- loo_compare(w1, w2)
  expect_identical(colnames(old), comp_colnames)
  expect_null(attr(old, "cluster_test"))

  # with ordinary row names (from a log-lik matrix) the numbers are unchanged
  M1 <- matrix(rnorm(prod(dim(LL)), c(LL), 0.01), ncol = dim(LL)[3])
  M2 <- matrix(rnorm(prod(dim(LL)), c(LL + 1), 0.01), ncol = dim(LL)[3])
  a <- suppressWarnings(loo(M1))
  b <- suppressWarnings(loo(M2))
  plain <- loo_compare(a, b)
  city <- rep(sprintf("c%d", 1:4), length.out = ncol(M1))
  with_test <- loo_compare(a, b, clusters = city)
  expect_equal(plain$model, with_test$model)
  expect_equal(plain$elpd_diff, with_test$elpd_diff)
  expect_equal(plain$se_diff, with_test$se_diff)
  expect_equal(plain$p_worse, with_test$p_worse)
  expect_identical(
    colnames(plain),
    c("model", "elpd_diff", "se_diff", "p_worse", "diag_diff", "diag_elpd",
      "elpd_loo", "se_elpd_loo", "p_loo", "se_p_loo", "looic", "se_looic")
  )
  expect_true(all(c("stat_joint", "p_joint", "p_joint_adj") %in%
                    colnames(with_test)))
})

test_that("cluster test validates its inputs", {
  dat <- make_cluster_loos()
  x <- dat$loos[[1]]
  y <- dat$loos[[2]]

  expect_error(
    loo_compare(x, y, clusters = dat$cluster, cluster = dat$cluster),
    "at most one of"
  )
  # too many clusters
  many <- rep(sprintf("c%02d", 1:19), length.out = nrow(x$pointwise))
  expect_error(loo_compare(x, y, clusters = many), "at most 18 clusters")
  # wrong length (positional specification)
  expect_error(
    loo_compare(x, y, clusters = dat$cluster[-1]),
    "one cluster label per observation"
  )
  # missing cluster labels
  bad <- dat$cluster
  bad[1] <- NA
  expect_error(loo_compare(x, y, clusters = bad), "missing cluster labels")
  # duplicate transaction IDs
  dup <- data.frame(id = c(dat$ids[1], dat$ids), city = c("a", dat$cluster))
  expect_error(loo_compare(x, y, clusters = dup), "duplicate transaction IDs")
  # unknown transaction IDs
  unknown <- data.frame(id = paste0(dat$ids, "x"), city = dat$cluster)
  expect_error(loo_compare(x, y, clusters = unknown), "are found in the pointwise")
  # invalid reference
  expect_error(
    loo_compare(x, y, clusters = dat$cluster, reference = "nope"),
    "name of one of the models"
  )
  # reference without clusters
  expect_error(loo_compare(x, y, reference = 1L),
               "only used together with 'clusters'")
  # IDs are required on the pointwise estimates when keyed by ID
  no_ids <- y
  rownames(no_ids$pointwise) <- NULL
  df <- data.frame(id = rev(dat$ids), city = rev(dat$cluster))
  expect_error(loo_compare(x, no_ids, clusters = df), "row names")
  # positional specification is still allowed (first model's row order)
  expect_silent(loo_compare(x, no_ids, clusters = dat$cluster))
  # wrong type
  expect_error(loo_compare(x, y, clusters = matrix(dat$cluster, ncol = 1)),
               "data frame")
})

test_that("cluster test is not available for subsampled loo", {
  LL <- example_loglik_array()
  w1 <- suppressWarnings(waic(LL))
  w2 <- suppressWarnings(waic(LL + 1))
  ss1 <- w1
  ss2 <- w2
  class(ss1) <- c("psis_loo_ss", class(ss1))
  class(ss2) <- c("psis_loo_ss", class(ss2))
  expect_error(
    loo_compare(ss1, ss2, clusters = rep("a", nrow(w1$pointwise))),
    "not available for subsampled"
  )
})

# Tests for deprecated compare() ------------------------------------------

test_that("compare throws deprecation warnings", {
  expect_warning(loo::compare(w1, w2), "Deprecated")
  expect_warning(loo::compare(w1, w1, w2), "Deprecated")
})

test_that("compare returns expected result (2 models)", {
  expect_warning(comp1 <- loo::compare(w1, w1), "Deprecated")
  expect_equal(comp1[1:2], c(elpd_diff = 0, se = 0))

  expect_warning(comp2 <- loo::compare(w1, w2), "Deprecated")
  expect_equal(round(comp2[1:2], 3), c(elpd_diff = -4.057, se = 0.088))
  expect_s3_class(comp2, "old_compare.loo")

  # specifying objects via ... and via arg x gives equal results
  expect_warning(comp_via_list <- loo::compare(x = list(w1, w2)), "Deprecated")
  expect_equal(comp2, comp_via_list)
})

test_that("compare returns expected result (3 models)", {
  w3 <- suppressWarnings(waic(LLarr3))
  expect_warning(comp1 <- loo::compare(w1, w2, w3), "Deprecated")

  expect_equal(
    colnames(comp1),
    c(
      "elpd_diff",
      "se_diff",
      "elpd_waic",
      "se_elpd_waic",
      "p_waic",
      "se_p_waic",
      "waic",
      "se_waic"
    )
  )
  expect_equal(rownames(comp1), c("w1", "w2", "w3"))
  expect_equal(comp1[1, 1], 0)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "matrix")
  expect_snapshot_value(comp1, style = "serialize")

  # specifying objects via '...' gives equivalent results (equal
  # except rownames) to using 'x' argument
  expect_warning(
    comp_via_list <- loo::compare(x = list(w1, w2, w3)),
    "Deprecated"
  )
  expect_equal(comp1, comp_via_list, ignore_attr = TRUE)
})

test_that("compare throws appropriate errors", {
  expect_error(
    suppressWarnings(loo::compare(w1, w2, x = list(w1, w2))),
    "should not be specified"
  )
  expect_error(suppressWarnings(loo::compare(x = 2)), "must be a list")
  expect_error(
    suppressWarnings(loo::compare(x = list(2))),
    "should have class 'loo'"
  )
  expect_error(
    suppressWarnings(loo::compare(x = list(w1))),
    "requires at least two models"
  )

  w3 <- suppressWarnings(waic(LLarr2[,, -1]))
  expect_error(
    suppressWarnings(loo::compare(x = list(w1, w3))),
    "same number of data points"
  )
  expect_error(
    suppressWarnings(loo::compare(x = list(w1, w2, w3))),
    "same number of data points"
  )
})
