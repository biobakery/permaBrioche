#' Hajek Estimator for Repeated Measures (Subject-Respecting Split)
#'
#' @description
#' Implements a subject-level train/test split for a binary covariate and computes a
#' Hajek estimator on the test set. The distance-like outcome \eqn{Y} is constructed as
#' the (Bray-Curtis or Euclidean) distance of each test sample to the centroid of group 1
#' computed in the training set. Two core engines are provided:
#' (1) a constant-propensity version and (2) a logistic-regression propensity version.
#'
#' @param formula A model formula \code{D ~ X1 + X2 + ...}. \code{D} can be a
#'   \code{\link{dist}} object (it is checked for alignment but not used in the
#'   Hajek computation).
#' @param data A \code{data.frame} with rownames matching \code{bugs} and columns for
#'   the RHS variables and the \code{blocking_variable}.
#' @param bugs A \code{matrix} or \code{data.frame} with features (e.g., taxa) in
#'   columns and samples in rows. Rownames must match \code{data}.
#' @param sample_id Optional character scalar giving the column in \code{data}
#'   that contains sample identifiers matching the distance matrix labels.
#'   If \code{NULL}, rownames(data) are used.
#' @param blocking_variable Character; name of the subject/cluster column in \code{data}.
#' @param covariate_name Character; the name of the **binary** covariate used as treatment.
#' @param permutations Integer; number of permutations for the null (default 999).
#' @param split_ratio Proportion of subjects assigned to training (default 0.5).
#' @param numerical_metadata Optional character vector of numeric RHS covariate names
#'   to use for logistic propensity score estimation. If \code{NULL}, a constant e is used.
#' @param method One of \code{"bray"} (default) or \code{"euclidean"} for \eqn{Y}.
#'
#' @return A list with
#' \itemize{
#'   \item \code{observed}: observed Hajek estimate on the test set
#'   \item \code{null_dist}: vector of valid null estimates (NAs removed)
#'   \item \code{ratio}: \code{|observed / mean(null)|}
#'   \item \code{z_score}: \code{(observed - mean(null)) / sd(null)}
#'   \item \code{pval}: permutation p-value using absolute deviations
#' }
#'
#' @details
#' The subject-level split ensures no leakage: centroids are computed on train subjects,
#' test outcomes are evaluated against that centroid, and the Hajek estimator is computed
#' on test samples only. The permutation null respects subject blocking (via
#' \code{permute::how(blocks=...)}).
#'
#' If \code{numerical_metadata} is \code{NULL}, a constant propensity (\eqn{e}) equal
#' to the mean of \code{T} in the test set is used; otherwise a logit model
#' \code{T ~ numerical_metadata} is fitted on the training set and predicted on test.
#'
#' @examples
#' ## ---- minimal runnable example ----
#' set.seed(1)
#'
#' n_subj <- 6
#' reps   <- 2
#' N      <- n_subj * reps
#'
#' subject <- factor(rep(seq_len(n_subj), each = reps))
#' trt     <- factor(rep(c("A","B"), length.out = N))
#' X       <- rnorm(N)
#'
#' meta <- data.frame(
#'   subject = subject,
#'   trt = trt,
#'   X = X,
#'   row.names = paste0("id", seq_len(N))
#' )
#'
#' bugs <- matrix(
#'   rexp(N * 10),
#'   nrow = N,
#'   dimnames = list(rownames(meta), NULL)
#' )
#'
#' D <- vegan::vegdist(bugs, method = "bray")
#'
#' hajek_repeat_measures(
#'   D ~ trt,
#'   data = meta,
#'   bugs = bugs,
#'   blocking_variable = "subject",
#'   covariate_name = "trt",
#'   permutations = 49
#' )
#'
#' \donttest{
#'   ## logistic-propensity version
#'   hajek_repeat_measures(
#'     D ~ trt + X,
#'     data = meta,
#'     bugs = bugs,
#'     blocking_variable = "subject",
#'     covariate_name = "trt",
#'     numerical_metadata = "X",
#'     permutations = 99
#'   )
#' }
#'
#' @import vegan
#' @import permute
#' @importFrom stats as.formula aggregate binomial glm predict terms var
#' @importFrom dplyr full_join filter summarise across mutate pull
#' @importFrom tibble rownames_to_column
#' @export
hajek_repeat_measures <- function(formula,
                                  data,
                                  bugs,
                                  sample_id = NULL, 
                                  blocking_variable = "subject",
                                  covariate_name,
                                  permutations   = 999,
                                  split_ratio    = 0.5,
                                  numerical_metadata = NULL,
                                  method         = "bray") {
  data <- as.data.frame(data)
  
  # --- 0. Checks ---
  if (missing(covariate_name))
    stop("Please provide covariate_name (the binary treatment/covariate).")
  
  # --- 1. Parse formula and (optionally) extract distance object on LHS ---
  YVAR <- formula[[2]]
  lhs  <- eval(YVAR, environment(formula), globalenv())
  environment(formula) <- environment()
  D <- lhs  # not used in the core, but kept for symmetry/validation
  
  # --- 2. Align distance matrix samples with metadata ---
  d_labels <- attr(D, "Labels")
  if (is.null(d_labels))
    stop("Distance object has no Labels attribute")
  
  if (is.null(sample_id)) {
    # Use rownames(data)
    if (is.null(rownames(data)))
      stop("data has no rownames and sample_id is NULL")
    
    if (!all(d_labels %in% rownames(data)))
      stop("Not all distance labels are present in rownames(data). Consider specifying sample_id")
    
    data <- data[d_labels, , drop = FALSE]
    
  } else {
    # Use explicit sample_id column
    if (!sample_id %in% colnames(data))
      stop(sprintf("sample_id '%s' not found in data", sample_id))
    
    if (any(is.na(data[[sample_id]])))
      stop("sample_id column contains NA values")
    
    if (!all(d_labels %in% data[[sample_id]]))
      stop("Not all distance labels are present in data[[sample_id]]")
    
    if (anyDuplicated(data[[sample_id]]))
      stop("sample_id column must contain unique sample identifiers")
    
    data <- data[match(d_labels, data[[sample_id]]), , drop = FALSE]
    rownames(data) <- d_labels
  }
  
  # --- 3. Subset to blocking variable + RHS terms ---
  rhs_terms <- labels(terms(formula))
  data_sub  <- data[, c(blocking_variable, rhs_terms), drop = FALSE]
  
  if (any(is.na(data_sub)))
    stop("data must not have NA values for Hajek (or pre-filter before calling).")
  
  # --- 4. Define blocks (one level per subject / cluster) ---
  blocks <- as.factor(data_sub[[blocking_variable]])
  
  # --- 5. Identify which RHS vars are static within each block ---
  rhs_only <- data_sub[, rhs_terms, drop = FALSE]
  agg_res <- aggregate(
    rhs_only,
    list(block = data_sub[[blocking_variable]]),
    function(x) length(unique(x)) == 1
  )
  static_vars <- sapply(agg_res[, -1, drop = FALSE], all)  # logical named vector
  
  # --- 6. Split into permute_within (varying) and block_data_full (static) ---
  permute_within  <- rhs_only[, names(static_vars)[!static_vars], drop = FALSE]
  block_data_full <- rhs_only[, names(static_vars)[static_vars],  drop = FALSE]
  
  # Edge cases
  if (ncol(permute_within) == 0L) {
    permute_within <- data.frame(row.names = rownames(data_sub))
  }
  if (ncol(block_data_full) == 0L) {
    block_data_full <- as.data.frame(matrix(0, nrow = length(levels(blocks)), ncol = 0))
    rownames(block_data_full) <- levels(blocks)
  }
  
  # --- 7. Reduce block_data to one row per block level ---
  block_data <- block_data_full[!duplicated(blocks), , drop = FALSE]
  rownames(block_data) <- levels(blocks)
  
  # --- 8. Sanity check: covariate_name must appear in RHS ---
  if (!(covariate_name %in% colnames(permute_within) ||
        covariate_name %in% colnames(block_data))) {
    stop(sprintf("covariate_name '%s' not found among RHS variables in the formula.",
                 covariate_name))
  }
  
  # --- 9. Call the appropriate core engine ---
  if (is.null(numerical_metadata)) {
    res <- hajek_repeat_measures_core_without_propensity(
      D                  = D,
      permute_within     = permute_within,
      blocks             = blocks,
      block_data         = block_data,
      permutations       = permutations,
      split_ratio        = split_ratio,
      numerical_metadata = numerical_metadata,
      bugs               = bugs,
      method             = method,
      covariate_name     = covariate_name
    )
  } else {
    res <- hajek_repeat_measures_core(
      D                  = D,
      permute_within     = permute_within,
      blocks             = blocks,
      block_data         = block_data,
      permutations       = permutations,
      split_ratio        = split_ratio,
      numerical_metadata = numerical_metadata,
      bugs               = bugs,
      method             = method,
      covariate_name     = covariate_name
    )
  }
  
  heading <- sprintf(
    paste("Hajek estimator with subject-level split",
          "Blocked by %s",
          "Number of permutations: %d\n", sep = "\n"),
    blocking_variable, permutations
  )
  attr(res, "heading") <- paste0(heading, paste0(deparse(sys.call()), collapse = "\n"))
  
  res
}

#' @keywords internal
#' @noRd
#' @importFrom permute shuffle how
hajek_repeat_measures_core_without_propensity <- function(D,
                                                          permute_within,
                                                          blocks = NULL,
                                                          block_data,
                                                          permutations = 999,
                                                          split_ratio = 0.5,
                                                          numerical_metadata,
                                                          bugs,
                                                          method = "bray",
                                                          covariate_name) {
  # Assemble metadata aligned with blocks
  if (is.null(blocks)) stop("Need 'blocks' (subject ids per sample).")
  mtdat <- cbind(permute_within, block_data[blocks, , drop = FALSE])
  
  # Ensure bugs aligns with mtdat rownames
  if (!all(rownames(mtdat) == rownames(bugs)))
    stop("Row names of 'bugs' must match permuted metadata.")
  
  # Factorize covariate and collect levels
  mtdat[[covariate_name]] <- factor(mtdat[[covariate_name]])
  cov_levels <- levels(mtdat[[covariate_name]])
  if (length(cov_levels) != 2) stop("Binary covariate required.")
  
  subject_ids <- unique(blocks)
  
  # --- Train/test split by subject; ensure both levels appear in train ---
  train_contains_both_levels <- function(subject_train) {
    train_rows <- which(blocks %in% subject_train)
    lv <- unique(mtdat[[covariate_name]][train_rows])
    length(lv) == 2
  }
  max_tries <- 2000L
  train_subjects <- NULL
  for (try in seq_len(max_tries)) {
    n_train_subjects <- floor(length(subject_ids) * split_ratio)
    cand <- sample(subject_ids, n_train_subjects)
    if (train_contains_both_levels(cand)) { train_subjects <- cand; break }
  }
  if (is.null(train_subjects))
    stop("Cannot find subject-level split with both covariate levels in training.")
  
  train_indices <- which(blocks %in% train_subjects)
  test_indices  <- setdiff(seq_len(nrow(mtdat)), train_indices)
  
  # Identify feature columns
  bugs_cols <- colnames(bugs)
  
  # --- Compute group-1 centroid from train ---
  train_lvl1 <- train_indices[mtdat[[covariate_name]][train_indices] == cov_levels[1]]
  if (length(train_lvl1) == 0L) stop("No group-1 samples in training set.")
  group1_centroid <- colMeans(bugs[train_lvl1, bugs_cols, drop = FALSE])
  
  # --- Outcome Y on test ---
  if (method == "euclidean") {
    diffs <- sweep(bugs[test_indices, bugs_cols, drop = FALSE], 2L, group1_centroid, FUN = "-")
    Y <- rowSums(diffs * diffs)
  } else if (method == "bray") {
    Xtest <- bugs[test_indices, bugs_cols, drop = FALSE]
    den   <- rowSums(Xtest + matrix(group1_centroid, nrow = nrow(Xtest), ncol = length(group1_centroid), byrow = TRUE))
    num   <- rowSums(abs(Xtest - matrix(group1_centroid, nrow = nrow(Xtest), ncol = length(group1_centroid), byrow = TRUE)))
    Y     <- ifelse(den == 0, 0, num / den)
  } else {
    stop("Unsupported method. Use 'bray' or 'euclidean'.")
  }
  
  # --- Treatment T and constant propensity e on test ---
  Ttest <- ifelse(mtdat[[covariate_name]][test_indices] == cov_levels[1], 0, 1)
  if (all(Ttest == 0) || all(Ttest == 1))
    stop("Test set has only one treatment level; adjust split_ratio or data.")
  etest <- mean(Ttest)
  
  get_tau_hat <- function(T, Y, e) {
    treat_num <- sum(T * Y / e)
    treat_den <- sum(T / e)
    ctrl_num  <- sum((1 - T) * Y / (1 - e))
    ctrl_den  <- sum((1 - T) / (1 - e))
    (treat_num / treat_den) - (ctrl_num / ctrl_den)
  }
  
  observed_tau <- get_tau_hat(Ttest, Y, etest)
  
  # --- Permutation null ---
  ctrl <- how(blocks = blocks)
  null_taus <- replicate(permutations, {
    within.i <- shuffle(nrow(permute_within), control = ctrl)
    block.i  <- sample(seq_len(nrow(block_data)))
    mtdat_perm <- cbind(
      permute_within[within.i, , drop = FALSE],
      block_data[block.i, , drop = FALSE][blocks, , drop = FALSE]
    )
    
    # Train/test split stays the same (subject-respecting)
    # Recompute centroid and outcomes under permuted covariate
    cov_perm <- factor(mtdat_perm[[covariate_name]], levels = cov_levels)
    train_lvl1_p <- train_indices[cov_perm[train_indices] == cov_levels[1]]
    if (length(train_lvl1_p) == 0L) return(NA_real_)
    
    g1_centroid_p <- colMeans(bugs[train_lvl1_p, bugs_cols, drop = FALSE])
    
    if (method == "euclidean") {
      diffs_p <- sweep(bugs[test_indices, bugs_cols, drop = FALSE], 2L, g1_centroid_p, FUN = "-")
      Yp <- rowSums(diffs_p * diffs_p)
    } else {
      Xtest <- bugs[test_indices, bugs_cols, drop = FALSE]
      denp  <- rowSums(Xtest + matrix(g1_centroid_p, nrow = nrow(Xtest), ncol = length(g1_centroid_p), byrow = TRUE))
      nump  <- rowSums(abs(Xtest - matrix(g1_centroid_p, nrow = nrow(Xtest), ncol = length(g1_centroid_p), byrow = TRUE)))
      Yp    <- ifelse(denp == 0, 0, nump / denp)
    }
    
    Ttest_p <- ifelse(cov_perm[test_indices] == cov_levels[1], 0, 1)
    if (all(Ttest_p == 0) || all(Ttest_p == 1)) return(NA_real_)
    ep <- mean(Ttest_p)
    
    get_tau_hat(Ttest_p, Yp, ep)
  })
  
  null_taus_valid <- null_taus[!is.na(null_taus)]
  mean_null <- mean(null_taus_valid)
  sd_null   <- stats::sd(null_taus_valid)
  ratio     <- abs(observed_tau / mean_null)
  z_score   <- (observed_tau - mean_null) / sd_null
  pval      <- (sum(abs(null_taus_valid) >= abs(observed_tau)) + 1) /
    (length(null_taus_valid) + 1)
  
  list(observed = observed_tau, null_dist = null_taus_valid,
       ratio = ratio, z_score = z_score, pval = pval)
}

#' @keywords internal
#' @noRd
#' @importFrom stats glm predict
#' @importFrom permute shuffle how
hajek_repeat_measures_core <- function(D,
                                       permute_within,
                                       blocks = NULL,
                                       block_data,
                                       permutations = 999,
                                       split_ratio = 0.5,
                                       numerical_metadata,
                                       bugs,
                                       method = "bray",
                                       covariate_name) {
  if (is.null(blocks)) stop("Need 'blocks' (subject ids per sample).")
  
  # Metadata assembly
  mtdat <- as.data.frame(cbind(permute_within, block_data[blocks, , drop = FALSE]),
                         check.names = FALSE)
  
  # Sanity: covariate
  mtdat[[covariate_name]] <- factor(mtdat[[covariate_name]])
  cov_levels <- levels(mtdat[[covariate_name]])
  if (length(cov_levels) != 2) stop("Binary covariate required.")
  
  if (!all(rownames(mtdat) == rownames(bugs)))
    stop("Row names of 'bugs' must match permuted metadata.")
  bugs_cols <- colnames(bugs)
  
  subject_ids <- unique(blocks)
  
  # Helper: estimate e(T=1|X) via logistic regression on specified numeric metadata
  estimate_e <- function(df_train, df_test, Ttrain, num_cols) {
    if (length(num_cols) == 0) {
      # fallback: constant e on test
      e_test <- rep(mean(df_test$T), nrow(df_test))
    } else {
      # Build frames for glm: include only numeric metadata columns
      Xtr <- df_train[, num_cols, drop = FALSE]
      Xte <- df_test[,  num_cols, drop = FALSE]
      fit_e <- stats::glm(T ~ ., data = cbind(T = Ttrain, Xtr), family = binomial())
      e_test <- stats::predict(fit_e, newdata = Xte, type = "response")
    }
    eps <- 1e-3
    e_test[e_test < eps] <- eps
    e_test[e_test > 1 - eps] <- 1 - eps
    e_test
  }
  
  # Subject-level train/test split
  train_contains_both_levels <- function(subject_train) {
    train_rows <- which(blocks %in% subject_train)
    lv <- unique(mtdat[[covariate_name]][train_rows])
    length(lv) == 2
  }
  max_tries <- 2000L
  train_subjects <- NULL
  for (try in seq_len(max_tries)) {
    n_train_subjects <- floor(length(subject_ids) * split_ratio)
    cand <- sample(subject_ids, n_train_subjects)
    if (train_contains_both_levels(cand)) { train_subjects <- cand; break }
  }
  if (is.null(train_subjects))
    stop("Cannot find subject-level split with both covariate levels in training.")
  
  train_indices <- which(blocks %in% train_subjects)
  test_indices  <- setdiff(seq_len(nrow(mtdat)), train_indices)
  
  # Outcome Y (distance to group-1 centroid)
  train_lvl1 <- train_indices[mtdat[[covariate_name]][train_indices] == cov_levels[1]]
  if (length(train_lvl1) == 0L) stop("No group-1 samples in training set.")
  group1_centroid <- colMeans(bugs[train_lvl1, bugs_cols, drop = FALSE])
  
  if (method == "euclidean") {
    diffs <- sweep(bugs[test_indices, bugs_cols, drop = FALSE], 2L, group1_centroid, FUN = "-")
    Ytest <- rowSums(diffs * diffs)
  } else if (method == "bray") {
    Xtest <- bugs[test_indices, bugs_cols, drop = FALSE]
    den   <- rowSums(Xtest + matrix(group1_centroid, nrow = nrow(Xtest), ncol = length(group1_centroid), byrow = TRUE))
    num   <- rowSums(abs(Xtest - matrix(group1_centroid, nrow = nrow(Xtest), ncol = length(group1_centroid), byrow = TRUE)))
    Ytest <- ifelse(den == 0, 0, num / den)
  } else stop("Unsupported method. Use 'bray' or 'euclidean'.")
  
  # Treatment and e
  Ttrain <- ifelse(mtdat[[covariate_name]][train_indices] == cov_levels[1], 0, 1)
  Ttest  <- ifelse(mtdat[[covariate_name]][test_indices]  == cov_levels[1], 0, 1)
  if (all(Ttest == 0) || all(Ttest == 1))
    stop("Test set has only one treatment level; adjust split_ratio or data.")
  
  df_train <- mtdat[train_indices, , drop = FALSE]
  df_test  <- mtdat[test_indices,  , drop = FALSE]
  num_cols <- if (is.null(numerical_metadata)) character(0) else numerical_metadata
  e_test   <- estimate_e(df_train, df_test, Ttrain, num_cols)
  
  get_tau_hat <- function(T, Y, e) {
    treat_num <- sum(T * Y / e)
    treat_den <- sum(T / e)
    ctrl_num  <- sum((1 - T) * Y / (1 - e))
    ctrl_den  <- sum((1 - T) / (1 - e))
    (treat_num / treat_den) - (ctrl_num / ctrl_den)
  }
  
  observed_tau <- get_tau_hat(Ttest, Ytest, e_test)
  
  # Permutation null
  ctrl <- how(blocks = blocks)
  null_taus <- replicate(permutations, {
    within.i <- shuffle(nrow(permute_within), control = ctrl)
    block.i  <- sample(seq_len(nrow(block_data)))
    
    mtdat_perm <- as.data.frame(cbind(
      permute_within[within.i, , drop = FALSE],
      block_data[block.i, , drop = FALSE][blocks, , drop = FALSE]
    ), check.names = FALSE)
    
    cov_perm <- factor(mtdat_perm[[covariate_name]], levels = cov_levels)
    
    # Recompute centroid under permutation
    train_lvl1_p <- train_indices[cov_perm[train_indices] == cov_levels[1]]
    if (length(train_lvl1_p) == 0L) return(NA_real_)
    g1_centroid_p <- colMeans(bugs[train_lvl1_p, bugs_cols, drop = FALSE])
    
    if (method == "euclidean") {
      diffs_p <- sweep(bugs[test_indices, bugs_cols, drop = FALSE], 2L, g1_centroid_p, FUN = "-")
      Yp <- rowSums(diffs_p * diffs_p)
    } else {
      Xtest <- bugs[test_indices, bugs_cols, drop = FALSE]
      denp  <- rowSums(Xtest + matrix(g1_centroid_p, nrow = nrow(Xtest), ncol = length(g1_centroid_p), byrow = TRUE))
      nump  <- rowSums(abs(Xtest - matrix(g1_centroid_p, nrow = nrow(Xtest), ncol = length(g1_centroid_p), byrow = TRUE)))
      Yp    <- ifelse(denp == 0, 0, nump / denp)
    }
    
    Ttest_p <- ifelse(cov_perm[test_indices] == cov_levels[1], 0, 1)
    if (all(Ttest_p == 0) || all(Ttest_p == 1)) return(NA_real_)
    
    # Recompute propensity on permuted data
    df_train_p <- mtdat_perm[train_indices, , drop = FALSE]
    df_test_p  <- mtdat_perm[test_indices,  , drop = FALSE]
    Ttrain_p   <- ifelse(cov_perm[train_indices] == cov_levels[1], 0, 1)
    e_test_p   <- estimate_e(df_train_p, df_test_p, Ttrain_p, num_cols)
    
    get_tau_hat(Ttest_p, Yp, e_test_p)
  })
  
  null_taus_valid <- null_taus[!is.na(null_taus)]
  mean_null <- mean(null_taus_valid)
  sd_null   <- stats::sd(null_taus_valid)
  ratio     <- abs(observed_tau / mean_null)
  z_score   <- (observed_tau - mean_null) / sd_null
  pval      <- (sum(abs(null_taus_valid) >= abs(observed_tau)) + 1) /
    (length(null_taus_valid) + 1)
  
  list(observed = observed_tau,
       null_dist = null_taus_valid,
       ratio = ratio,
       z_score = z_score,
       pval = pval)
}

#' Hajek with Location and Dispersion (Subject-Respecting Split)
#'
#' Computes (i) the Hajek total effect on a distance-like outcome to the group-1
#' centroid and (ii) decompositions for location and dispersion in feature space,
#' using a subject-level train/test split that prevents leakage.
#'
#' @param D Unused in the computation (kept for API symmetry).
#' @param permute_within Data frame of variables that vary within blocks.
#' @param blocks Factor or vector of subject/cluster ids (one per sample).
#' @param block_data Data frame of static block-level variables (rows = block levels).
#' @param permutations Not used (kept for API symmetry; default 0).
#' @param split_ratio Proportion of subjects in training (default 0.5).
#' @param numerical_metadata Unused (kept for API symmetry).
#' @param bugs Matrix/data.frame of features, rows=samples, cols=features.
#' @param method "bray" (default) or "euclidean" for the Y outcome.
#' @param covariate_name Character; binary covariate used to form groups.
#'
#' @return A list with
#' \itemize{
#'   \item \code{observed_tau}: Hajek total effect (test set)
#'   \item \code{observed_loc}: location component (difference of squared distances to centroid)
#'   \item \code{observed_disp}: dispersion component (difference of traces of covariance)
#' }
#'
#' @examples
#' ## ---- minimal runnable sketch ----
#' set.seed(1)
#' n <- 4
#' permute_within <- data.frame(trt = factor(rep(c("A","B"), length.out = n)))
#' blocks <- factor(rep(1:2, each = 2))
#' block_data <- data.frame(Z = c(0, 1))
#' bugs <- matrix(rexp(n * 5), nrow = n)
#'
#' hajek_repeat_measures_loc_and_disp(
#'   D = NULL,
#'   permute_within = permute_within,
#'   blocks = blocks,
#'   block_data = block_data,
#'   bugs = bugs,
#'   covariate_name = "trt"
#' )
#' @export
hajek_repeat_measures_loc_and_disp <- function(D,
                                               permute_within,
                                               blocks = NULL,
                                               block_data,
                                               permutations = 0,
                                               split_ratio = 0.5,
                                               numerical_metadata,
                                               bugs,
                                               method = "euclidean",
                                               covariate_name) {
  if (is.null(blocks)) stop("Provide 'blocks' (subject ids per sample).")
  
  # --- Assemble metadata aligned to bugs ---
  if (nrow(permute_within) != nrow(bugs)) {
    stop("permute_within must have the same number of rows as bugs.")
  }
  
  # If permute_within has no rownames, inherit from bugs
  if (is.null(rownames(permute_within))) {
    rownames(permute_within) <- rownames(bugs)
  }
  
  # Assemble full metadata
  mtdat <- cbind(
    permute_within,
    block_data[blocks, , drop = FALSE]
  )
  
  # --- Enforce alignment between metadata and bugs ---
  
  # Case 1: permute_within had no rownames → assume row order matches bugs
  if (is.null(rownames(mtdat))) {
    rownames(mtdat) <- rownames(bugs)
  }
  
  # Case 2: exact match (fast path)
  if (identical(rownames(mtdat), rownames(bugs))) {
    # OK
  } else if (setequal(rownames(mtdat), rownames(bugs))) {
    # Case 3: same samples, different order → reorder metadata to bugs
    mtdat <- mtdat[rownames(bugs), , drop = FALSE]
  } else {
    # Case 4: incompatible
    stop(
      "Metadata rows do not align with 'bugs'. ",
      "Row names must either match exactly, be absent (row-order alignment), ",
      "or be a permutation of bugs row names."
    )
  }
  
  
  mtdat[[covariate_name]] <- factor(mtdat[[covariate_name]])
  cov_levels <- levels(mtdat[[covariate_name]])
  if (length(cov_levels) != 2) stop("Binary covariate required.")
  
  subject_ids <- unique(blocks)
  
  train_contains_both_levels <- function(subject_train) {
    train_rows <- which(blocks %in% subject_train)
    lv <- unique(mtdat[[covariate_name]][train_rows])
    length(lv) == 2
  }
  
  max_tries <- 2000L
  train_subjects <- NULL
  for (try in seq_len(max_tries)) {
    n_train_subjects <- floor(length(subject_ids) * split_ratio)
    cand <- sample(subject_ids, n_train_subjects)
    if (train_contains_both_levels(cand)) { train_subjects <- cand; break }
  }
  if (is.null(train_subjects))
    stop("Subject-level split failed; adjust split_ratio or data.")
  
  train_indices <- which(blocks %in% train_subjects)
  test_indices  <- setdiff(seq_len(nrow(mtdat)), train_indices)
  
  # Centroid from train (level 1)
  train_lvl1 <- train_indices[mtdat[[covariate_name]][train_indices] == cov_levels[1]]
  if (length(train_lvl1) == 0L) stop("No group-1 samples in training set.")
  group1_centroid <- colMeans(bugs[train_lvl1, , drop = FALSE])
  
  # Outcomes on test
  if (method == "euclidean") {
    diffs <- sweep(bugs[test_indices, , drop = FALSE], 2L, group1_centroid, FUN = "-")
    Y <- rowSums(diffs * diffs)
  } else if (method == "bray") {
    Xtest <- bugs[test_indices, , drop = FALSE]
    den   <- rowSums(Xtest + matrix(group1_centroid, nrow = nrow(Xtest),
                                    ncol = ncol(Xtest), byrow = TRUE))
    num   <- rowSums(abs(Xtest - matrix(group1_centroid, nrow = nrow(Xtest),
                                        ncol = ncol(Xtest), byrow = TRUE)))
    Y     <- ifelse(den == 0, 0, num / den)
  } else stop("Unsupported method.")
  
  # Treatment on test
  Ttest <- ifelse(mtdat[[covariate_name]][test_indices] == cov_levels[1], 0, 1)
  if (all(Ttest == 0) || all(Ttest == 1))
    stop("Test set has only one treatment level; adjust split_ratio or data.")
  
  # Constant propensity score in test
  etest <- mean(Ttest)
  
  # Hajek + location + dispersion
  get_tau_hat <- function(T, Y, e) {
    treat_num <- sum(T * Y / e)
    treat_den <- sum(T / e)
    ctrl_num  <- sum((1 - T) * Y / (1 - e))
    ctrl_den  <- sum((1 - T) / (1 - e))
    (treat_num / treat_den) - (ctrl_num / ctrl_den)
  }
  
  tau_hat <- get_tau_hat(Ttest, Y, etest)
  
  # Location & dispersion (feature space; unweighted diagnostics)
  x0 <- bugs[test_indices[Ttest == 0], , drop = FALSE]
  x1 <- bugs[test_indices[Ttest == 1], , drop = FALSE]
  
  if (nrow(x0) < 2 || nrow(x1) < 2) {
    tau_loc  <- NA_real_
    tau_disp <- NA_real_
  } else {
    mu0 <- colMeans(x0); mu1 <- colMeans(x1)
    tr_Sigma0 <- sum(apply(x0, 2, stats::var))
    tr_Sigma1 <- sum(apply(x1, 2, stats::var))
    # location: difference in squared distance to reference centroid
    tau_loc  <- sum((mu1 - group1_centroid)^2) - sum((mu0 - group1_centroid)^2)
    # dispersion: difference in total variance (trace)
    tau_disp <- tr_Sigma1 - tr_Sigma0
  }
  
  list(
    observed_tau  = tau_hat,
    observed_loc  = tau_loc,
    observed_disp = tau_disp
  )
}
