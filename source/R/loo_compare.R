#' Model comparison
#'
#' @description Compare fitted models based on [ELPD][loo-glossary].
#'
#' @export
#' @param x An object of class `"loo"` or a list of such objects. If a list is
#'   used then the list names will be used as the model names in the output. See
#'   **Examples**.
#' @param ... Additional objects of class `"loo"`, if not passed in as a single
#'   list.
#' @param clusters Optional specification of clusters (for example cities) for
#'   a joint cluster-level sign-flip test of the pointwise ELPD differences. Can be a data frame with a transaction ID column (`id`,
#'   `transaction_id`, or the first column) and a cluster column (`cluster`,
#'   `city`, `group`, or the second column), a named vector of cluster labels
#'   where the names are the transaction IDs, or an unnamed vector of cluster
#'   labels with one element per observation. When the cluster specification is
#'   keyed by transaction ID (a data frame or a named vector), the models are
#'   aligned by transaction ID using the row names of the pointwise estimates,
#'   so the row order of each model's pointwise estimates does not matter. An
#'   unnamed cluster vector is interpreted in the row order of the **first**
#'   model; the other models are still aligned by transaction ID when they have
#'   unique row names. At most 18 clusters are allowed because all \eqn{2^K}
#'   sign patterns are enumerated. See **Details**.
#' @param cluster Alias for `clusters`. At most one of `clusters` and `cluster`
#'   may be specified.
#' @param reference Optional reference model for the cluster-level joint test.
#'   This is only used when `clusters` is supplied and does not change the
#'   usual `elpd_diff`/`se_diff` columns, which are always relative to the model
#'   with the best ELPD. Can be a model name (a name in `x` or supplied in
#'   `...`), a positive integer index in the order the models were supplied, or
#'   a `"loo"` object equal to one of the models. The joint test compares every
#'   other (candidate) model to this reference. The default (`NULL`) uses the
#'   first model supplied, which is convenient when the already deployed model
#'   is passed first (e.g., `loo_compare(deployed, candidate1, candidate2, ...)`).
#'
#' @return A data frame with class `"compare.loo"` that has its own
#'   print method. See the **Details** and **Examples** sections. When
#'   `clusters` is supplied, the data frame also has the columns `stat_joint`
#'   (the observed cluster-level statistic), `p_joint` (the raw sign-flip
#'   p-value), and `p_joint_adj` (the p-value adjusted for selecting the best
#'   candidate). The attribute `"cluster_test"` is a list holding the cluster
#'   labels, the number of clusters, the number of sign patterns, the reference
#'   model, the observed statistics, both p-values, and the \eqn{K \times M}
#'   matrix of cluster-level differences.
#'
#' @details
#'   When comparing two fitted models, we can estimate the difference in their
#'   expected predictive accuracy by the difference in
#'   [`elpd_loo`][loo-glossary] or `elpd_waic` (or multiplied by \eqn{-2}, if
#'   desired, to be on the deviance scale).
#'
#' ## `elpd_diff` and `se_diff`
#'   When using `loo_compare()`, the returned data frame will have one row per
#'   model and several columns of estimates. The values of
#'   [`elpd_diff`][loo-glossary] and [`se_diff`][loo-glossary] are computed by
#'   making pairwise comparisons between each model and the model with the
#'   largest ELPD (the model listed first). Therefore, the first `elpd_diff`
#'   value will always be `0` (i.e., the difference between the preferred model
#'   and itself) and the rest of the values will be negative.
#'
#'   To compute the standard error of the difference in [ELPD][loo-glossary] ---
#'   which should not be expected to equal the difference of the standard errors
#'   --- we use a paired estimate to take advantage of the fact that the same
#'   set of \eqn{N} data points was used to fit both models. These calculations
#'   should be most useful when \eqn{N} is large, because then non-normality of
#'   the distribution is not such an issue when estimating the uncertainty in
#'   these sums. These standard errors, for all their flaws, should give a
#'   better sense of uncertainty than what is obtained using the current
#'   standard approach of comparing differences of deviances to a Chi-squared
#'   distribution, a practice derived for Gaussian linear models or
#'   asymptotically, and which only applies to nested models in any case.
#'
#' ## `p_worse`, `diag_diff`, and `diag_elpd`
#'   The values in the `p_worse` column show the probability of each model
#'   having worse ELPD than the best model. These probabilities are computed
#'   with a normal approximation using the values from `elpd_diff` and
#'   `se_diff`. Sivula et al. (2025) present the conditions when the normal
#'   approximation used for SE and `se_diff` is good, and the column
#'   `diag_diff` contains possible diagnostic messages:
#'
#'   * `N < 100` (small data)
#'   * `|elpd_diff| < 4` (models make similar predictions)
#'
#'   If either of these diagnostic messages is shown, the error distribution is
#'   skewed or thick tailed and the normal approximation based on `elpd_diff`
#'   and `se_diff` is not well calibrated. In that case, the probabilities
#'   `p_worse` are likely to be too large. However, `elpd_diff` and `se_diff`
#'   will still be indicative of the differences and uncertainties (for example,
#'   if `|elpd_diff|` is many times larger than `se_diff` the difference is quite
#'   certain). In addition, if the model is not well specificed and there are
#'   outliers, the error distribution can also be skewed or thick tailed and the
#'   normal approximation is not well calibrated. Possible model misspecification
#'   and outliers can be diagnosed with usual predictive checking methods.
#'
#'   The column `diag_elpd` shows the PSIS-LOO Pareto k diagnostic for the
#'   pointwise ELPD computations for each model. If `K k_psis > 0.7` is shown,
#'   where `K` is the number of high Pareto k values in the PSIS
#'   computation, then there may be significant bias in `elpd_diff` favoring
#'   models with a large number of high Pareto k values.
#'
#' ## Warnings for many model comparisons
#'   If more than \eqn{11} models are compared, we internally recompute the model
#'   differences using the median model by ELPD as the baseline model. We then
#'   estimate whether the differences in predictive performance are potentially
#'   due to chance as described by McLatchie and Vehtari (2023). This will flag
#'   a warning if it is deemed that there is a risk of over-fitting due to the
#'   selection process. In that case users are recommended to avoid model
#'   selection based on LOO-CV, and instead to favor model averaging/stacking or
#'   projection predictive inference.
#'
#' ## Cluster-level joint sign-flip test
#'   The usual `elpd_diff` and `se_diff` columns assume that the pointwise ELPD
#'   contributions are independent. When observations are grouped (for example
#'   several transactions from the same city), this assumption can be badly
#'   violated. If `clusters` is supplied, `loo_compare()` additionally performs
#'   a joint sign-flip test at the cluster level.
#'
#'   The models are first aligned by transaction ID (the row names of the
#'   pointwise estimates). For each model \eqn{j} relative to the reference
#'   model, the pointwise differences are summed within each cluster \eqn{k},
#'   giving one cluster-level difference \eqn{D_{kj}}{D[k,j]}. The one-sided
#'   statistic for model \eqn{j} is
#'   \deqn{T_j = \frac{\sum_k D_{kj}}{\sqrt{\sum_k D_{kj}^2}},}{T[j] = sum_k D[k,j] / sqrt(sum_k D[k,j]^2),}
#'   and is defined to be `0` when all cluster differences are zero. A positive
#'   value favours model \eqn{j} over the reference model.
#'
#'   All \eqn{2^K} sign patterns \eqn{s \in \{-1,+1\}^K}{s in {-1,+1}^K} are
#'   enumerated, where \eqn{K} is the number of clusters. The **same** sign is
#'   applied to a given cluster for **all** candidate models, which preserves
#'   the dependence between candidates. The raw p-value is the fraction of sign
#'   patterns with a statistic at least as large as the observed one. The
#'   adjusted p-value compares the observed statistic for model \eqn{j} to the
#'   **maximum** statistic over all candidate models after each sign flip, and
#'   is the fraction of patterns for which this maximum is at least as large as
#'   the observed statistic. In both cases patterns equal to the observed one
#'   are counted, and the identity pattern is included in the reference
#'   distribution. The adjustment controls the family-wise error rate over the
#'   candidate models (Westfall and Young, 1993).
#'
#'   The full sign-flip-by-candidate matrix is never stored: because
#'   \eqn{\sqrt{\sum_k D_{kj}^2}} does not depend on the signs, each sign pattern
#'   only requires the matrix-vector product \eqn{s^\top D}{s'D}, and running
#'   counts and maxima are updated on the fly. With 18 clusters there are
#'   262,144 sign patterns.
#'
#'   The identity pattern uses the observed data, so when the observed statistic
#'   itself is the maximum the p-value is at least \eqn{1 / 2^K}. This is the
#'   usual resolution limit of an exhaustive sign-flip test.
#'
#'   The reference model defaults to the first model supplied (typically the
#'   model already in production). It can be changed with the `reference`
#'   argument. The reference model has no statistic of its own, so `stat_joint`,
#'   `p_joint`, and `p_joint_adj` are `NA` for its row; the other rows are the
#'   candidate models. Columns are only added when `clusters` is supplied, so
#'   the usual comparison output is unchanged otherwise.
#'
#' @seealso
#' * The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html) on
#'   the __loo__ website for answers to frequently asked questions.
#' @template loo-and-compare-references
#'
#' @examples
#' # very artificial example, just for demonstration!
#' LL <- example_loglik_array()
#' loo1 <- loo(LL)     # should be worst model when compared
#' loo2 <- loo(LL + 1) # should be second best model when compared
#' loo3 <- loo(LL + 2) # should be best model when compared
#'
#' comp <- loo_compare(loo1, loo2, loo3)
#' print(comp, digits = 2)
#' print(comp, simplify = FALSE) # full table
#'
#' # can use a list of objects with custom names
#' # the names will be used in the output
#' loo_compare(list("apple" = loo1, "banana" = loo2, "cherry" = loo3))
#'
#' # optional cluster-level joint sign-flip test. Here we pretend the 32
#' # observations belong to 4 cities and that the pointwise estimates carry
#' # transaction IDs in their row names.
#' rn <- sprintf("tx%02d", seq_len(nrow(loo1$pointwise)))
#' add_ids <- function(x) {
#'   rownames(x$pointwise) <- rn
#'   x
#' }
#' loo1 <- add_ids(loo1)
#' loo2 <- add_ids(loo2)
#' loo3 <- add_ids(loo3)
#' city <- rep(c("A", "B", "C", "D"), length.out = length(rn))
#' comp_joint <- loo_compare(loo1, loo2, loo3, clusters = city)
#' print(comp_joint)
#'
#' # a data frame keyed by transaction ID can be given in any row order
#' clusters <- data.frame(id = rev(rn), city = rev(city))
#' loo_compare(loo1, loo2, loo3, clusters = clusters)
#'
#' # by default the first model supplied is the reference; use 'reference'
#' # to compare the candidates against the deployed model instead
#' loo_compare(loo1, loo2, loo3, clusters = city, reference = "model1")
#'
#' \dontrun{
#' # works for waic (and kfold) too
#' loo_compare(waic(LL), waic(LL - 10))
#' }
#'
loo_compare <- function(x, ...) {
  UseMethod("loo_compare")
}

#' @rdname loo_compare
#' @export
loo_compare.default <- function(x, ..., clusters = NULL, cluster = NULL,
                                reference = NULL) {
  dots <- list(...)

  if (!is.null(clusters) && !is.null(cluster)) {
    stop("Specify at most one of 'clusters' and 'cluster'.", call. = FALSE)
  }
  if (is.null(clusters)) {
    clusters <- cluster
  }

  if (is.loo(x)) {
    loos <- c(list(x), dots)
  } else {
    if (!is.list(x) || !length(x)) {
      stop("'x' must be a list if not a 'loo' object.")
    }
    if (length(dots)) {
      stop("If 'x' is a list then '...' should not be specified.")
    }
    loos <- x
  }

  # if subsampling is used
  if (any(sapply(loos, inherits, "psis_loo_ss"))) {
    if (!is.null(clusters) || !is.null(reference)) {
      stop("The cluster sign-flip test is not available for subsampled ",
           "'loo' objects.", call. = FALSE)
    }
    return(loo_compare.psis_loo_ss_list(loos))
  }

  # run pre-comparison checks
  loo_compare_checks(loos)

  if (is.null(clusters) && !is.null(reference)) {
    stop("'reference' is only used together with 'clusters'.", call. = FALSE)
  }

  # normalize the cluster specification and align the models by transaction ID
  # (this also makes the usual comparison invariant to row shuffling). The
  # reference model index is resolved before alignment; alignment applies the
  # same sample ordering to every model, so the index stays valid.
  spec <- NULL
  ref_idx <- NULL
  if (!is.null(clusters)) {
    spec <- loo_cluster_spec(clusters)
    ref_idx <- loo_cluster_reference(loos, reference)
    aligned <- loo_cluster_align(loos, spec)
    loos <- aligned$loos
    spec <- aligned$spec
  }

  # the joint test is computed on the models as supplied (before ordering)
  joint <- NULL
  if (!is.null(spec)) {
    joint <- loo_compare_cluster_test(loos, spec, reference = ref_idx)
  }

  # compute elpd_diff and se_elpd_diff relative to best model
  comp <- loo_compare_matrix(loos)
  ord <- loo_compare_order(loos)
  rnms <- rownames(comp)
  diffs <- mapply(FUN = elpd_diffs, loos[ord[1]], loos[ord])
  colnames(diffs) <- rnms
  elpd_diff <- apply(diffs, 2, sum)
  se_diff <- apply(diffs, 2, se_elpd_diff)

  # compute probabilities that a model has worse elpd than the best model
  # using a normal approximation (Sivula et al., 2025)
  p_worse <- stats::pnorm(0, elpd_diff, se_diff)
  p_worse[elpd_diff == 0] <- NA

  comp_dat <- data.frame(
    model = rnms,
    elpd_diff = elpd_diff,
    se_diff = se_diff,
    p_worse = p_worse,
    diag_diff = diag_diff(nrow(diffs), elpd_diff),
    diag_elpd = diag_elpd(loos[ord])
  )
  if (!is.null(joint)) {
    comp_dat$stat_joint <- unname(joint$stat[ord])
    comp_dat$p_joint <- unname(joint$p_raw[ord])
    comp_dat$p_joint_adj <- unname(joint$p_adj[ord])
  }

  comp <- cbind(comp_dat, as.data.frame(comp))
  rownames(comp) <- NULL

  # run order statistics-based checks for many model comparisons
  loo_order_stat_check(loos, ord)

  class(comp) <- c("compare.loo", class(comp))
  if (!is.null(joint)) {
    attr(comp, "cluster_test") <- joint[c(
      "cluster", "n_clusters", "n_sign_flips", "reference",
      "stat", "p_raw", "p_adj", "diffs", "candidates"
    )]
  }
  comp
}

#' @rdname loo_compare
#' @export
#' @param digits For the print method only, the number of digits to use when
#'   printing.
#' @param p_worse For the print method only, should we include the normal
#'   approximation based probability of each model having worse performance than
#'   the best model? The default is `TRUE`.
#' @param simplify For the print method only, should the output be simplified to
#'   only include the model names, ELPD differences, and (when `p_worse = TRUE`)
#'   diagnostic columns? The default is `TRUE`. Set to `FALSE` to also print the
#'   available estimate columns (pointwise ELPD, LOOIC/WAIC, and their standard
#'   errors).
print.compare.loo <- function(x, ..., digits = 1, p_worse = TRUE,
                              simplify = TRUE) {
  if (inherits(x, "old_compare.loo")) {
    return(unclass(x))
  }
  if (!inherits(x, "data.frame")) {
    class(x) <- c(class(x), "data.frame")
  }
  if (!all(c("model", "elpd_diff", "se_diff") %in% colnames(x))) {
    print(as.data.frame(x))
    return(x)
  }
  base_cols <- c("model", "elpd_diff", "se_diff")
  diag_cols <- c("p_worse", "diag_diff", "diag_elpd")
  show_diag <- p_worse && "p_worse" %in% colnames(x)
  joint_cols <- c("stat_joint", "p_joint", "p_joint_adj")
  joint_cols <- intersect(joint_cols, colnames(x))

  estimate_cols <- setdiff(colnames(x), c(base_cols, diag_cols, joint_cols))
  estimate_cols <- estimate_cols[vapply(x[estimate_cols], is.numeric, logical(1))]

  cols <- c(
    base_cols,
    if (show_diag) diag_cols,
    joint_cols,
    if (!simplify) estimate_cols
  )
  cols <- intersect(cols, colnames(x))

  x2 <- x[, cols, drop = FALSE]

  fmt_cols <- setdiff(cols, c("model", "diag_diff", "diag_elpd", joint_cols))
  if (length(fmt_cols)) {
    if ("p_worse" %in% fmt_cols) {
      x2$p_worse <- .fr(x2$p_worse, digits = 2)
      fmt_cols <- setdiff(fmt_cols, "p_worse")
    }
    if (length(fmt_cols)) {
      x2[fmt_cols] <- .fr(x2[fmt_cols], digits)
    }
  }
  if (length(joint_cols)) {
    if ("stat_joint" %in% joint_cols) {
      x2[["stat_joint"]] <- .fr(x2[["stat_joint"]], digits = 2)
    }
    for (nm in intersect(c("p_joint", "p_joint_adj"), joint_cols)) {
      x2[[nm]] <- format(signif(x2[[nm]], 3), trim = TRUE)
    }
  }
  # Use `as.data.frame(x2)` here to drop "compare.loo"
  # so print() uses print.data.frame.
  print(as.data.frame(x2), quote = FALSE, row.names = FALSE)

  # show glossary for diagnostic flags
  has_diag <- any(nzchar(x[["diag_diff"]], keepNA = FALSE), na.rm = TRUE) ||
              any(nzchar(x[["diag_elpd"]], keepNA = FALSE), na.rm = TRUE)
  if (has_diag && p_worse) {
    message(
      "\nDiagnostic flags present.\n",
      "See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)\n",
      "or https://mc-stan.org/loo/reference/loo-glossary.html."
    )
  }
  invisible(x)
}


# internal ----------------------------------------------------------------

#' Compute pointwise elpd differences
#' @noRd
#' @param loo_a,loo_b Two `"loo"` objects.
elpd_diffs <- function(loo_a, loo_b) {
  pt_a <- loo_a$pointwise
  pt_b <- loo_b$pointwise
  elpd <- grep("^elpd", colnames(pt_a))
  pt_b[, elpd] - pt_a[, elpd]
}

#' Normalize a cluster specification
#'
#' Accepts a data frame with an ID column and a cluster column, a named vector
#' of cluster labels (names are the transaction IDs), or an unnamed vector with
#' one cluster label per observation.
#'
#' @noRd
#' @param clusters The user-supplied cluster specification.
#' @return A list with `ids` (transaction IDs or `NULL`), `labels` (one cluster
#'   label per observation), `cluster` (the unique cluster labels) and
#'   `n_clusters`.
loo_cluster_spec <- function(clusters) {
  id_names <- c(
    "id", "ID", "transaction_id", "transaction_ids", "transaction",
    "transactions", "obs_id", "obs_ids", "observation_id", "observations"
  )
  cluster_names <- c(
    "cluster", "clusters", "cluster_id", "cluster_ids", "city", "cities",
    "city_id", "city_ids", "group", "groups", "group_id", "group_ids"
  )

  if (is.data.frame(clusters)) {
    nms <- names(clusters)
    id_col <- which(nms %in% id_names)
    cl_col <- which(nms %in% cluster_names)
    if (!length(id_col)) {
      id_col <- 1L
    } else {
      id_col <- id_col[1L]
    }
    cl_col <- setdiff(cl_col, id_col)
    if (!length(cl_col)) {
      cl_col <- setdiff(seq_along(nms), id_col)[1L]
    } else {
      cl_col <- cl_col[1L]
    }
    if (is.na(cl_col)) {
      stop("'clusters' must have columns for the transaction ID and the ",
           "cluster.", call. = FALSE)
    }
    ids <- clusters[[id_col]]
    labels <- clusters[[cl_col]]
  } else if (is.atomic(clusters) && is.null(dim(clusters))) {
    if (!is.null(names(clusters)) && all(nzchar(names(clusters)))) {
      ids <- names(clusters)
      labels <- unname(clusters)
    } else {
      ids <- NULL
      labels <- clusters
    }
  } else {
    stop("'clusters' must be a data frame, a named vector of cluster labels, ",
         "or a vector with one cluster label per observation.", call. = FALSE)
  }

  if (anyNA(labels)) {
    stop("'clusters' contains missing cluster labels.", call. = FALSE)
  }
  if (!is.null(ids)) {
    if (anyNA(ids) || any(!nzchar(as.character(ids)))) {
      stop("'clusters' contains missing transaction IDs.", call. = FALSE)
    }
    ids <- as.character(ids)
    if (anyDuplicated(ids)) {
      stop("'clusters' contains duplicate transaction IDs.", call. = FALSE)
    }
  }
  labels <- as.character(labels)
  if (!is.null(ids) && length(ids) != length(labels)) {
    stop("'clusters' has different numbers of transaction IDs and cluster ",
         "labels.", call. = FALSE)
  }

  lev <- unique(labels)
  if (length(lev) > 18L) {
    stop("The cluster sign-flip test allows at most 18 clusters, but 'clusters'",
         " has ", length(lev), ". All 2^K sign patterns are enumerated.",
         call. = FALSE)
  }
  if (length(lev) < 1L) {
    stop("'clusters' does not contain any clusters.", call. = FALSE)
  }

  list(ids = ids, labels = labels, cluster = lev, n_clusters = length(lev),
       keyed = !is.null(ids))
}

#' Align models by transaction ID and match the cluster specification to them
#'
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param spec A normalized cluster specification from [loo_cluster_spec()].
#' @return A list with the (possibly reordered) `loos` and the updated `spec`.
loo_cluster_align <- function(loos, spec) {
  N <- nrow(loos[[1L]]$pointwise)

  if (!isTRUE(spec$keyed)) {
    # positional specification: the labels follow the row order of the first
    # model. If every model carries unique row names we can still align the
    # other models by ID, which makes the test invariant to shuffling their
    # rows.
    if (length(spec$labels) != N) {
      stop("'clusters' must have one cluster label per observation (", N,
           "), but has ", length(spec$labels), ".", call. = FALSE)
    }
    ids <- lapply(loos, function(lo) rownames(lo$pointwise))
    usable <- vapply(
      ids,
      function(rn) !is.null(rn) && !anyDuplicated(rn) && all(nzchar(rn)),
      logical(1)
    )
    if (!all(usable)) {
      return(list(loos = loos, spec = spec))
    }
    spec$ids <- ids[[1L]]
    spec$keyed <- TRUE
  }

  if (length(spec$ids) != N) {
    stop("'clusters' describes ", length(spec$ids), " observations by ",
         "transaction ID, but the models have ", N, " observations.",
         call. = FALSE)
  }

  loos <- lapply(loos, function(lo) {
    rn <- rownames(lo$pointwise)
    if (is.null(rn)) {
      stop("When 'clusters' identifies observations by transaction ID, all ",
           "models must have row names on their pointwise estimates.",
           call. = FALSE)
    }
    if (anyDuplicated(rn)) {
      stop("The pointwise row names (transaction IDs) of a model are not ",
           "unique.", call. = FALSE)
    }
    idx <- match(spec$ids, rn)
    if (anyNA(idx)) {
      stop("Not all transaction IDs in 'clusters' are found in the pointwise ",
           "estimates of every model.", call. = FALSE)
    }
    lo$pointwise <- lo$pointwise[idx, , drop = FALSE]
    lo
  })

  list(loos = loos, spec = spec)
}

#' Resolve the reference model for the cluster test
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param reference A model name, index, `"loo"` object, or `NULL`.
#' @return The index of the reference model in `loos`.
loo_cluster_reference <- function(loos, reference) {
  nms <- find_model_names(loos)
  if (is.null(reference)) {
    return(1L)
  }
  if (is.loo(reference)) {
    hit <- which(vapply(
      loos,
      function(x) isTRUE(all.equal(x$pointwise, reference$pointwise)),
      logical(1)
    ))
    if (!length(hit)) {
      stop("'reference' is not one of the models being compared.", call. = FALSE)
    }
    return(hit[1L])
  }
  if (is.character(reference) && length(reference) == 1L) {
    hit <- match(reference, nms)
    if (is.na(hit)) {
      stop("'reference' must be the name of one of the models.", call. = FALSE)
    }
    return(hit)
  }
  if (is.numeric(reference) && length(reference) == 1L &&
      !is.na(reference) && reference == as.integer(reference) &&
      reference >= 1L && reference <= length(loos)) {
    return(as.integer(reference))
  }
  stop("'reference' must be a model name, a model index, or a 'loo' object.",
       call. = FALSE)
}

#' Cluster-level joint sign-flip test
#'
#' Sums the pointwise ELPD differences within clusters and enumerates all
#' \eqn{2^K} sign patterns. Only running counts and maxima are kept, so the
#' full sign-pattern-by-candidate statistic matrix is never stored.
#'
#' @noRd
#' @param loos List of `"loo"` objects aligned by transaction ID.
#' @param spec Normalized cluster specification.
#' @param reference Reference model (see [loo_cluster_reference()]).
#' @return A list with the observed statistics, raw and adjusted p-values, and
#'   the cluster-level differences.
loo_compare_cluster_test <- function(loos, spec, reference = NULL) {
  nms <- find_model_names(loos)
  ref_idx <- loo_cluster_reference(loos, reference)
  lev <- spec$cluster
  K <- length(lev)
  M <- length(loos)

  # one elpd column per model (they are aligned already)
  elpd <- lapply(loos, function(lo) {
    pt <- lo$pointwise
    col <- grep("^elpd", colnames(pt))
    if (!length(col)) {
      stop("The models must have pointwise ELPD estimates.", call. = FALSE)
    }
    pt[, col[1L]]
  })
  ref_elpd <- elpd[[ref_idx]]

  D <- matrix(0, nrow = K, ncol = M, dimnames = list(lev, nms))
  for (j in seq_len(M)) {
    agg <- rowsum(elpd[[j]] - ref_elpd, group = spec$labels, reorder = TRUE)
    D[match(rownames(agg), lev), j] <- agg[, 1L]
  }

  denom <- sqrt(colSums(D^2))
  inv_denom <- ifelse(denom > 0, 1 / denom, 0)
  stat <- colSums(D) * inv_denom
  stat[denom == 0] <- 0

  candidates <- setdiff(seq_len(M), ref_idx)
  if (!length(candidates)) {
    stop("The cluster sign-flip test requires at least one candidate model in ",
         "addition to the reference model.", call. = FALSE)
  }

  # streaming enumeration of all 2^K sign patterns, in blocks small enough
  # that no full sign-pattern-by-candidate matrix is ever materialized.
  # The statistic is recomputed as a matrix-vector product, so values that are
  # mathematically equal to the observed statistic may differ in the last
  # bits. A tiny absolute tolerance (the normalized statistic is bounded by
  # sqrt(K)) makes ties count as "at least as large", as required.
  n_sign <- as.integer(2^K)
  block <- max(1L, min(n_sign, as.integer(65536 / max(K, 1L))))
  raw_count <- numeric(M)
  adj_count <- numeric(M)
  cand_stat <- stat[candidates]
  tol <- 1e-10 * max(1, sqrt(K))
  cand_thr <- cand_stat - tol
  nc <- length(candidates)

  starts <- seq.int(0L, n_sign - 1L, by = block)
  for (start in starts) {
    b <- as.integer(seq.int(start, min(start + block - 1L, n_sign - 1L)))
    nb <- length(b)
    S <- matrix(-1, nrow = nb, ncol = K)
    for (k in seq_len(K)) {
      S[, k] <- ifelse(bitwAnd(bitwShiftR(b, k - 1L), 1L) == 1L, 1, -1)
    }
    # s^T D for each pattern in the block (the denominators are sign-invariant)
    Tv <- sweep(S %*% D, 2L, inv_denom, "*")
    Tv_cand <- Tv[, candidates, drop = FALSE]
    ge <- Tv_cand >= matrix(cand_thr, nrow = nb, ncol = nc, byrow = TRUE)
    raw_count[candidates] <- raw_count[candidates] + colSums(ge)
    block_max <- matrixStats::rowMaxs(Tv_cand)
    adj_count[candidates] <- adj_count[candidates] +
      colSums(outer(block_max, cand_thr, ">="))
  }

  p_raw <- raw_count / n_sign
  p_adj <- adj_count / n_sign

  # the reference model is not a candidate, so it has no statistic
  stat_out <- setNames(rep(NA_real_, M), nms)
  p_raw_out <- setNames(rep(NA_real_, M), nms)
  p_adj_out <- setNames(rep(NA_real_, M), nms)
  stat_out[candidates] <- stat[candidates]
  p_raw_out[candidates] <- p_raw[candidates]
  p_adj_out[candidates] <- p_adj[candidates]

  list(
    diffs = D,
    cluster = lev,
    n_clusters = K,
    n_sign_flips = n_sign,
    reference = nms[ref_idx],
    reference_index = ref_idx,
    stat = stat_out,
    p_raw = p_raw_out,
    p_adj = p_adj_out,
    candidates = nms[candidates]
  )
}

#' Compute pointwise elpd differences
#' @noRd
#' @param loo_a,loo_b Two `"loo"` objects.
elpd_diffs <- function(loo_a, loo_b) {
  pt_a <- loo_a$pointwise
  pt_b <- loo_b$pointwise
  elpd <- grep("^elpd", colnames(pt_a))
  pt_b[, elpd] - pt_a[, elpd]
}

#' Compute standard error of the elpd difference
#' @noRd
#' @param diffs Vector of pointwise elpd differences
se_elpd_diff <- function(diffs) {
  N <- length(diffs)
  # As `elpd_diff` is defined as the sum of N independent components,
  # we can compute the standard error by using the standard deviation
  # of the N components and multiplying by `sqrt(N)`.
  sqrt(N) * sd(diffs)
}

#' Perform checks on `"loo"` objects before comparison
#' @noRd
#' @param loos List of `"loo"` objects.
#' @return Nothing, just possibly throws errors/warnings.
loo_compare_checks <- function(loos) {
  ## errors
  if (length(loos) <= 1L) {
    stop("'loo_compare' requires at least two models.", call.=FALSE)
  }
  if (!all(sapply(loos, is.loo))) {
    stop("All inputs should have class 'loo'.", call.=FALSE)
  }

  Ns <- vapply(loos, function(x) nrow(x$pointwise), integer(1))
  if (any(Ns != Ns[1L])) {
    stop(
      paste0(
        "All models must have the same number of observations, but models have inconsistent observation counts: ",
        paste(paste0("'", find_model_names(loos), "' (", Ns, ")"), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ## warnings

  yhash <- lapply(loos, attr, which = "yhash")
  yhash_ok <- sapply(yhash, function(x) { # ok only if all yhash are same (all NULL is ok)
    isTRUE(all.equal(x, yhash[[1]]))
  })
  if (!all(yhash_ok)) {
    .warn_diagnostic("Not all models have the same y variable. ('yhash' attributes do not match)",
                     class = "loo_warning_y_mismatch")
  }

  if (all(sapply(loos, is.kfold))) {
    Ks <- unlist(lapply(loos, attr, which = "K"))
    if (!all(Ks == Ks[1])) {
      .warn_diagnostic("Not all kfold objects have the same K value. ",
                       "For a more accurate comparison use the same number of folds. ",
                       class = "loo_warning_cv_mismatch")
    }
  } else if (any(sapply(loos, is.kfold)) && any(sapply(loos, is.psis_loo))) {
    .warn_diagnostic("Comparing LOO-CV to K-fold-CV. ",
                     "For a more accurate comparison use the same number of folds ",
                     "or loo for all models compared.",
                     class = "loo_warning_cv_mismatch")
  }
}


#' Find the model names associated with `"loo"` objects
#'
#' @export
#' @param x List of `"loo"` objects.
#' @return Character vector of model names the same length as `x.`
#'
find_model_names <- function(x) {
  stopifnot(is.list(x))
  out_names <- character(length(x))

  names1 <- names(x)
  names2 <- lapply(x, "attr", "model_name", exact = TRUE)
  names3 <- lapply(x, "[[", "model_name")
  names4 <- paste0("model", seq_along(x))

  for (j in seq_along(x)) {
    if (isTRUE(nzchar(names1[j]))) {
      out_names[j] <- names1[j]
    } else if (length(names2[[j]])) {
      out_names[j] <- names2[[j]]
    } else if (length(names3[[j]])) {
      out_names[j] <- names3[[j]]
    } else {
      out_names[j] <- names4[j]
    }
  }
  out_names
}


#' Compute the loo_compare matrix
#' @noRd
#' @param loos List of `"loo"` objects.
loo_compare_matrix <- function(loos){
  tmp <- sapply(loos, function(x) {
    est <- x$estimates
    setNames(c(est), nm = c(rownames(est), paste0("se_", rownames(est))))
  })
  colnames(tmp) <- find_model_names(loos)
  rnms <- rownames(tmp)
  comp <- tmp
  ord <- loo_compare_order(loos)
  comp <- t(comp)[ord, ]
  patts <- c("elpd", "p_", "^waic$|^looic$", "^se_waic$|^se_looic$")
  col_ord <- unlist(sapply(patts, function(p) grep(p, colnames(comp))),
                    use.names = FALSE)
  comp <- comp[, col_ord]
  comp
}

#' Computes the order of loos for comparison
#' @noRd
#' @param loos List of `"loo"` objects.
loo_compare_order <- function(loos){
  tmp <- sapply(loos, function(x) {
    est <- x$estimates
    setNames(c(est), nm = c(rownames(est), paste0("se_", rownames(est))))
  })
  colnames(tmp) <- find_model_names(loos)
  rnms <- rownames(tmp)
  ord <- order(tmp[grep("^elpd", rnms), ], decreasing = TRUE)
  ord
}

#' Perform checks on `"loo"` objects __after__ comparison
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param ord List of `"loo"` object orderings.
#' @return Nothing, just possibly throws errors/warnings.
loo_order_stat_check <- function(loos, ord) {

  ## breaks

  if (length(loos) <= 11L) {
    # procedure cannot be diagnosed for fewer than ten candidate models
    # (total models = worst model + ten candidates)
    # break from function
    return(NULL)
  }

  ## warnings

  # compute the elpd differences from the median model
  baseline_idx <- middle_idx(ord)
  diffs <- mapply(FUN = elpd_diffs, loos[ord[baseline_idx]], loos[ord])
  elpd_diff <- apply(diffs, 2, sum)

  # estimate the standard deviation of the upper-half-normal
  diff_median <- stats::median(elpd_diff)
  elpd_diff_trunc <- elpd_diff[elpd_diff >= diff_median]
  n_models <- sum(!is.na(elpd_diff_trunc))
  candidate_sd <- sqrt(1 / n_models * sum(elpd_diff_trunc^2, na.rm = TRUE))

  # estimate expected best diff under null hypothesis
  K <- length(loos) - 1
  order_stat <- order_stat_heuristic(K, candidate_sd)

  if (max(elpd_diff) <= order_stat) {
    # flag warning if we suspect no model is theoretically better than the baseline
    .warn_diagnostic("Difference in performance potentially due to chance. ",
                     "See McLatchie and Vehtari (2023) for details.",
                     class = "loo_warning_chance_difference")
  }
}

#' Returns the middle index of a vector
#' @noRd
#' @param vec A vector.
#' @return Integer index value.
middle_idx <- function(vec) floor(length(vec) / 2)

#' Computes maximum order statistic from K Gaussians
#' @noRd
#' @param K Number of Gaussians.
#' @param c Scaling of the order statistic.
#' @return Numeric expected maximum from K samples from a Gaussian with mean
#' zero and scale `"c"`
order_stat_heuristic <- function(K, c) {
  qnorm(p = 1 - 1 / (K * 2), mean = 0, sd = c)
}

#' Count number of high Pareto k values in PSIS-LOO and create diagnostic message
#' @noRd
#' @param loos Ordered list of loo objects.
#' @return Character vector of diagnostic messages.
diag_elpd <- function(loos) {
  sapply(loos, function(loo) {
    k <- loo$diagnostics[["pareto_k"]]
    if (is.null(k)) {
      out <- ""
    } else {
      S <- dim(loo)[1]
      khat_threshold <- ps_khat_threshold(S)
      K <- sum(k > khat_threshold)
      out <- ifelse(K == 0, "", paste0(K, " k_psis > ", round(khat_threshold, 2)))
    }
    out
  })
}

#' Create diagnostic for elpd differences
#' @noRd
#' @param N Number of data points.
#' @param elpd_diff Vector of elpd differences.
#' @return Character vector of diagnostic messages.
diag_diff <- function(N, elpd_diff) {
  if (N < 100) {
    diag_diff <- rep("N < 100", length(elpd_diff))
    diag_diff[elpd_diff == 0] <- ""
  } else {
    diag_diff <- rep("", length(elpd_diff))
    diag_diff[elpd_diff > -4 & elpd_diff != 0] <- "|elpd_diff| < 4"
  }
  diag_diff
}
