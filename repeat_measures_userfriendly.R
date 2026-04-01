#repeat_measures_userfriendly.R

# ---------- PERMANOVA (repeat measures; macOS-safe) ----------
PERMANOVA_repeat_measures_core <- function(
    D, permute_within, blocks = NULL, block_data,
    permutations = 999,
    metadata_order = c(names(permute_within), names(block_data)),
    na.rm = FALSE) {
  
  if (class(D) != "dist") stop("D must be a dist object")
  
  if (!missing(block_data) && is.null(blocks)){
    stop("blocks must be given if block_data is present")
  } else if (is.null(blocks)) {
    blocks     <- rep(1, nrow(permute_within))
    block_data <- as.data.frame(matrix(0, nrow = 1, ncol = 0))
  } else if (length(unique(blocks)) == 1) {
    warning("blocks only contains one unique value")
  }
  
  if (length(intersect(names(permute_within), names(block_data))) > 0)
    stop("metadata is repeated across permute_within and block_data")
  
  if (length(setdiff(metadata_order, union(names(permute_within), names(block_data)))) > 0)
    stop("metadata_order contains metadata not in permute_within and block_data")
  
  ord <- rownames(as.matrix(D))
  if (length(ord) != nrow(permute_within) || length(blocks) != length(ord))
    stop("blocks, permute_within, and D are not the same size")
  
  if (is.null(rownames(permute_within))) {
    warning("permute_within has no rownames - can't verify sample orders")
  } else if (!all(ord == rownames(permute_within))) {
    stop("rownames do not match between permute_within and D")
  }
  
  if (any(is.na(blocks))) stop("NAs are not allowed in blocks")
  
  if (is.factor(blocks)) {
    if (any(!(levels(blocks) %in% rownames(block_data))))
      stop("not all block levels are contained in block_data")
    block_data <- block_data[match(levels(blocks), rownames(block_data)), , drop = FALSE]
    blocks     <- as.numeric(blocks)
  } else if (is.numeric(blocks)) {
    if (blocks < 1 || max(blocks) > nrow(block_data))
      stop("Numeric blocks has indices out of range")
  } else if (is.character(blocks)) {
    if (is.null(rownames(block_data)) || !all(blocks %in% rownames(block_data)))
      stop("blocks does not match the rownames of block_data")
    blocks <- match(blocks, rownames(block_data))
  } else {
    stop("blocks must be a numeric, factor, or character vector")
  }
  
  na.removed <- 0
  if (any(is.na(permute_within)) || any(is.na(block_data))) {
    if (na.rm) {
      n_prerm <- length(blocks)
      hasna   <- (rowSums(is.na(block_data)) > 0) |
        (sapply(split(rowSums(is.na(permute_within)) > 0, blocks), mean) == 1)
      block_data <- block_data[!hasna,, drop = FALSE]
      keep   <- !hasna[blocks]
      blocks <- cumsum(!hasna)[blocks]
      
      blocks <- blocks[keep]
      permute_within <- permute_within[keep,, drop = FALSE]
      D      <- as.matrix(D)[keep, keep]
      
      keep   <- rowSums(is.na(permute_within)) == 0
      blocks <- blocks[keep]
      permute_within <- permute_within[keep,, drop = FALSE]
      D      <- as.dist(D[keep, keep])
      
      if (length(blocks) < ncol(permute_within) + ncol(block_data)) {
        stop(sprintf("After omitting samples, samples (%d) < metadata (%d)",
                     length(blocks), ncol(permute_within) + ncol(block_data)))
      }
      na.removed <- n_prerm - length(blocks)
    } else {
      stop("Some metadata is NA! adonis does not support any NA in the metadata")
    }
  }
  
  mtdat <- cbind(permute_within, block_data[blocks,,drop=FALSE])
  ad    <- adonis2(D ~ ., permutations = 0, data = mtdat[, metadata_order, drop=FALSE])
  R2    <- ad$R2; names(R2) <- rownames(ad)
  
  nullsamples <- matrix(NA, nrow = length(R2), ncol = permutations)
  for (i in seq_len(permutations)) {
    within.i <- shuffle(nrow(permute_within), control = how(blocks=blocks))
    block.i  <- sample(seq_len(nrow(block_data)))
    mtdat.i  <- cbind(
      permute_within[within.i,,drop=FALSE],
      block_data[block.i,,drop=FALSE][blocks,,drop=FALSE]
    )
    perm.ad <- adonis2(D ~ ., permutations = 0, data = mtdat.i[, metadata_order, drop=FALSE])
    nullsamples[,i] <- perm.ad$R2
  }
  
  n <- length(R2)
  R2[n-1]           <- 1 - R2[n-1]
  nullsamples[n-1,] <- 1 - nullsamples[n-1,]
  
  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)
  P[n] <- NA
  ad$`Pr(>F)` <- P
  if (na.rm) ad$na.removed <- na.removed
  ad
}

PERMANOVA_repeat_measures <- function(formula,
                                      data,
                                      blocking_variable = "subject",
                                      permutations = 999,
                                      na.rm = FALSE) {
  # --- 1. Parse formula and extract distance object ---
  YVAR <- formula[[2]]
  lhs  <- eval(YVAR, environment(formula), globalenv())
  environment(formula) <- environment()
  if (!inherits(lhs, "dist"))
    stop("lhs of formula must be an adonis2-compatible 'dist' object")
  
  D <- lhs
  
  # --- 2. Make sure distance matrix matches data rownames ---
  if (!all(rownames(as.matrix(D)) == rownames(data)))
    stop("Row names of distance matrix must match row names of data")
  
  # --- 3. Subset to blocking variable + RHS terms ---
  rhs_terms <- labels(terms(formula))
  data_sub  <- data[, c(blocking_variable, rhs_terms), drop = FALSE]
  
  if (any(is.na(data_sub)) && !na.rm)
    stop("data must not have NA values (or set na.rm = TRUE)")
  
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
  
  # Optional: handle edge cases where one of them is empty
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
  
  # --- 8. Prepare metadata_order as in your core function ---
  metadata_order <- c(colnames(permute_within), colnames(block_data))
  
  # --- 9. Call your original (core) engine ---
  res <- PERMANOVA_repeat_measures_core(
    D              = D,
    permute_within = permute_within,
    blocks         = blocks,
    block_data     = block_data,
    permutations   = permutations,
    metadata_order = metadata_order,
    na.rm          = na.rm
  )
  
  # Optional: heading
  heading <- sprintf(paste("Permutation test for adonis under reduced model",
                           "Terms added sequentially (first to last)",
                           "Permutation: blocked by %s",
                           "Number of permutations: %d",
                           sep = "\n"),
                     blocking_variable, permutations)
  attr(res, "heading") <- paste0(heading, "\n", paste0(deparse(sys.call()), collapse = "\n"))
  
  res
}

# ---------- Hájek (subject-respecting split; no GUI) ----------
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
  library(dplyr)
  mtdat <- cbind(permute_within, block_data[blocks,,drop=FALSE])
  
  df <- mtdat %>%
    rownames_to_column("id") %>%
    full_join(bugs %>% rownames_to_column("id"), by = "id") %>%
    dplyr::select(-matches("\\.x$"), -matches("\\.y$"))
  
  mtdat[[covariate_name]] <- factor(mtdat[[covariate_name]])
  cov_levels <- levels(mtdat[[covariate_name]])
  stopifnot(length(cov_levels) == 2)
  
  if (is.null(blocks)) stop("Need 'blocks' (subject ids per sample).")
  subject_ids <- unique(blocks)
  
  train_contains_both_levels <- function(subject_train) {
    train_rows <- which(blocks %in% subject_train)
    lv <- unique(mtdat[[covariate_name]][train_rows])
    length(lv) == 2
  }
  
###
  # --- Balanced subject-level split requiring BOTH levels in train AND test ---
  
  max_tries <- 2000L
  subject_ids <- unique(blocks)
  n_train_subjects <- max(1L, floor(length(subject_ids) * split_ratio))
  
  found_split <- FALSE
  
  for (try in seq_len(max_tries)) {
    # 1) Sample train subjects; derive row indices
    train_subjects <- sample(subject_ids, n_train_subjects)
    train_indices  <- which(blocks %in% train_subjects)
    test_indices   <- setdiff(seq_len(nrow(df)), train_indices)
    
    # 2) Basic size sanity checks
    if (length(train_indices) < 2L || length(test_indices) < 2L) next
    
    # 3) Compute treatment indicators from metadata for train/test
    T_train <- ifelse(mtdat[[covariate_name]][train_indices] == cov_levels[1], 0L, 1L)
    T_test  <- ifelse(mtdat[[covariate_name]][test_indices]  == cov_levels[1], 0L, 1L)
    
    # 4) Require BOTH levels in BOTH splits
    if (any(T_train == 0L) && any(T_train == 1L) &&
        any(T_test  == 0L) && any(T_test  == 1L)) {
      
      # 5) Build train/test data frames
      df_train <- df[train_indices, , drop = FALSE]
      df_test  <- df[test_indices,  , drop = FALSE]
      
      # 6) Group-1 centroid from TRAIN only
      group1_centroid <- df_train %>%
        dplyr::filter(.data[[covariate_name]] == cov_levels[1]) %>%
        dplyr::summarise(dplyr::across(colnames(bugs), mean)) %>%
        as.numeric()
      
      # 7) Outcome Y on TEST
      if (method == "euclidean") {
        df_test$Y <- apply(df_test[, colnames(bugs), drop = FALSE], 1, function(row) {
          sum((row - group1_centroid)^2)
        })
      } else if (method == "bray") {
        df_test$Y <- apply(df_test[, colnames(bugs), drop = FALSE], 1, function(row) {
          den <- sum(row + group1_centroid); if (den == 0) return(0)
          sum(abs(row - group1_centroid)) / den
        })
      } else {
        stop("Unsupported method.")
      }
      
      # 8) Finalize treatment on TEST (reuse computed T_test)
      df_test$T <- T_test
      
      found_split <- TRUE
      break
    }
  }
  
  if (!found_split) {
    stop("Could not find a subject-level split with both treatment levels in ",
         "both train and test after max_tries. Consider adjusting split_ratio or data.")
  }
###
  
  df_test$e <- mean(df_test$T)
  
  get_tau_hat <- function(df_sub) {
    df_sub %>%
      summarise(
        treat_numerator   = sum(T * Y / e),
        treat_denominator = sum(T / e),
        control_numerator = sum((1 - T) * Y / (1 - e)),
        control_denominator = sum((1 - T) / (1 - e))
      ) %>%
      mutate(tau_hat = (treat_numerator / treat_denominator) -
               (control_numerator / control_denominator)) %>%
      pull(tau_hat)
  }
  
  observed_tau <- get_tau_hat(df_test)
  
  null_taus <- replicate(permutations, {
    within.i <- shuffle(nrow(permute_within), control = how(blocks = blocks))
    block.i  <- sample(seq_len(nrow(block_data)))
    mtdat_perm <- cbind(
      permute_within[within.i, , drop = FALSE],
      block_data[block.i, , drop = FALSE][blocks, , drop = FALSE]
    )
    
    df_perm <- df
    df_perm[[covariate_name]] <- mtdat_perm[[covariate_name]]
    
    df_train_perm <- df_perm[train_indices, , drop = FALSE]
    df_test_perm  <- df_perm[test_indices,  , drop = FALSE]
    
    if (length(unique(df_train_perm[[covariate_name]])) < 2) return(NA_real_)
    
    group1_centroid_perm <- df_train_perm %>%
      filter(get(covariate_name) == cov_levels[1]) %>%
      summarise(across(colnames(bugs), mean)) %>%
      as.numeric()
    
    if (method == "euclidean") {
      df_test_perm$Y <- apply(df_test_perm[, colnames(bugs)], 1, function(row) {
        sum((row - group1_centroid_perm)^2)
      })
    } else {
      df_test_perm$Y <- apply(df_test_perm[, colnames(bugs)], 1, function(row) {
        den <- sum(row + group1_centroid_perm); if (den == 0) return(0)
        sum(abs(row - group1_centroid_perm)) / den
      })
    }
    
    df_test_perm$T <- ifelse(df_test_perm[[covariate_name]] == cov_levels[1], 0, 1)
    if (all(df_test_perm$T == 0) || all(df_test_perm$T == 1)) return(NA_real_)
    df_test_perm$e <- mean(df_test_perm$T)
    
    get_tau_hat(df_test_perm)
  })
  
  null_taus_valid <- null_taus[!is.na(null_taus)]
  mean_null <- mean(null_taus_valid)
  sd_null   <- sd(null_taus_valid)
  ratio     <- abs(observed_tau / mean_null)
  z_score   <- (observed_tau - mean_null) / sd_null
  pval      <- (sum(abs(null_taus_valid) >= abs(observed_tau)) + 1) /
    (length(null_taus_valid) + 1)
  
  list(observed = observed_tau, null_dist = null_taus_valid,
       ratio = ratio, z_score = z_score, pval = pval)
}

# ---------- Hájek (subject-respecting split; logistic e) ----------
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
  library(dplyr)
  library(permute)
  
  # -----------------------------
  # 0. Metadata assembly
  # -----------------------------
  mtdat <- as.data.frame(cbind(permute_within, block_data[blocks, , drop = FALSE]))
  
  # start df as mtdat with an explicit 'id' column
  df <- as.data.frame(mtdat)
  df$id <- rownames(df)
  
  # add bugs
  bugs_df <- as.data.frame(bugs)
  bugs_df$id <- rownames(bugs_df)
  df <- full_join(df, bugs_df, by = "id")
  
  # numerical_metadata should be a character vector of column names
  if (!is.null(numerical_metadata)) {
    # sanity check: those columns must exist in df (they come from mtdat)
    if (!all(numerical_metadata %in% colnames(df))) {
      stop("Some numerical_metadata columns are not present in the data frame 'df'.")
    }
    num_cols <- numerical_metadata   # just store the column *names*
  } else {
    num_cols <- character(0)
  }
  
  mtdat[[covariate_name]] <- factor(mtdat[[covariate_name]])
  cov_levels <- levels(mtdat[[covariate_name]])
  stopifnot(length(cov_levels) == 2)
  
  if (is.null(blocks)) stop("Need 'blocks' (subject ids per sample).")
  subject_ids <- unique(blocks)
  
  bugs_cols <- colnames(bugs)
  
  # -----------------------------
  # 1. Helper: estimate e(T=1|X) via logistic regression
  # -----------------------------
  estimate_e <- function(df_train, df_test, covariate_name, num_cols) {
    if (length(num_cols) == 0) {
      # fallback: constant e
      #e_test <- rep(mean(df_train$T), nrow(df_test)) #old, was using train data
      e_test <- rep(mean(df_test$T), nrow(df_test))
      #stop("1.")
    } else {
      formula_e <- as.formula(
        paste("T ~", paste(num_cols, collapse = " + "))
      )
      fit_e <- glm(formula_e, data = df_train, family = binomial)
      e_test <- predict(fit_e, newdata = df_test, type = "response")
      #stop("2.")
    }
    # clip for numerical stability
    eps <- 1e-3
    e_test[e_test < eps]       <- eps
    e_test[e_test > 1 - eps]   <- 1 - eps
    e_test
  }
  
  # -----------------------------
  # 2. Subject-level train/test split
  # -----------------------------
  train_contains_both_levels <- function(subject_train) {
    train_rows <- which(blocks %in% subject_train)
    lv <- unique(mtdat[[covariate_name]][train_rows])
    length(lv) == 2
  }
  
###
  # --- Balanced subject-level split requiring BOTH levels in train AND test ---
  
  max_tries <- 2000L
  subject_ids <- unique(blocks)
  n_train_subjects <- max(1L, floor(length(subject_ids) * split_ratio))
  
  found_split <- FALSE
  
  for (try in seq_len(max_tries)) {
    # 1) Sample train subjects; derive row indices
    train_subjects <- sample(subject_ids, n_train_subjects)
    train_indices  <- which(blocks %in% train_subjects)
    test_indices   <- setdiff(seq_len(nrow(df)), train_indices)
    
    # 2) Basic size sanity checks
    if (length(train_indices) < 2L || length(test_indices) < 2L) next
    
    # 3) Compute treatment indicators from metadata for train/test
    T_train <- ifelse(mtdat[[covariate_name]][train_indices] == cov_levels[1], 0L, 1L)
    T_test  <- ifelse(mtdat[[covariate_name]][test_indices]  == cov_levels[1], 0L, 1L)
    
    # 4) Require BOTH levels in BOTH splits
    if (any(T_train == 0L) && any(T_train == 1L) &&
        any(T_test  == 0L) && any(T_test  == 1L)) {
      
      # 5) Build train/test data frames
      df_train <- df[train_indices, , drop = FALSE]
      df_test  <- df[test_indices,  , drop = FALSE]
      
      # 6) Group-1 centroid from TRAIN only
      group1_centroid <- df_train %>%
        dplyr::filter(.data[[covariate_name]] == cov_levels[1]) %>%
        dplyr::summarise(dplyr::across(colnames(bugs), mean)) %>%
        as.numeric()
      
      # 7) Outcome Y on TEST
      if (method == "euclidean") {
        df_test$Y <- apply(df_test[, colnames(bugs), drop = FALSE], 1, function(row) {
          sum((row - group1_centroid)^2)
        })
      } else if (method == "bray") {
        df_test$Y <- apply(df_test[, colnames(bugs), drop = FALSE], 1, function(row) {
          den <- sum(row + group1_centroid); if (den == 0) return(0)
          sum(abs(row - group1_centroid)) / den
        })
      } else {
        stop("Unsupported method.")
      }
      
      # 8) Finalize treatment on TEST (reuse computed T_test)
      df_test$T <- T_test
      
      found_split <- TRUE
      break
    }
  }
  
  if (!found_split) {
    stop("Could not find a subject-level split with both treatment levels in ",
         "both train and test after max_tries. Consider adjusting split_ratio or data.")
  }
###
  df_test$e <- estimate_e(df_train, df_test, covariate_name, num_cols)
  
  # -----------------------------
  # 5. Hájek estimator using sample-specific e_i
  # -----------------------------
  get_tau_hat <- function(df_sub) {
    with(df_sub, {
      treat_numerator     <- sum(T       * Y / e)
      treat_denominator   <- sum(T           / e)
      control_numerator   <- sum((1 - T) * Y / (1 - e))
      control_denominator <- sum((1 - T)     / (1 - e))
      (treat_numerator / treat_denominator) -
        (control_numerator / control_denominator)
    })
  }
  
  observed_tau <- get_tau_hat(df_test)
  
  # -----------------------------
  # 6. Permutation null
  # -----------------------------
  null_taus <- replicate(permutations, {
    within.i <- shuffle(nrow(permute_within), control = how(blocks = blocks))
    block.i  <- sample(seq_len(nrow(block_data)))
    
    mtdat_perm <- as.data.frame(cbind(
      permute_within[within.i, , drop = FALSE],
      block_data[block.i, , drop = FALSE][blocks, , drop = FALSE]
    ))
    
    df_perm <- as.data.frame(mtdat_perm)
    df_perm$id <- rownames(df_perm)
    
    bugs_df_p <- as.data.frame(bugs)
    bugs_df_p$id <- rownames(bugs_df_p)
    df_perm <- full_join(df_perm, bugs_df_p, by = "id")
    
    df_train_perm <- df_perm[train_indices, , drop = FALSE]
    df_test_perm  <- df_perm[test_indices,  , drop = FALSE]
    
    # If permuted train has only one level, skip
    if (length(unique(df_train_perm[[covariate_name]])) < 2) return(NA_real_)
    
    # Recompute group-1 centroid under permutation
    group1_centroid_perm <- df_train_perm %>%
      filter(get(covariate_name) == cov_levels[1]) %>%
      summarise(across(bugs_cols, mean)) %>%
      as.numeric()
    
    if (method == "euclidean") {
      df_test_perm$Y <- apply(df_test_perm[, bugs_cols], 1, function(row) {
        sum((row - group1_centroid_perm)^2)
      })
    } else {
      df_test_perm$Y <- apply(df_test_perm[, bugs_cols], 1, function(row) {
        den <- sum(row + group1_centroid_perm); if (den == 0) return(0)
        sum(abs(row - group1_centroid_perm)) / den
      })
    }
    
    df_train_perm$T <- ifelse(df_train_perm[[covariate_name]] == cov_levels[1], 0, 1)
    df_test_perm$T  <- ifelse(df_test_perm[[covariate_name]]  == cov_levels[1], 0, 1)
    if (all(df_test_perm$T == 0) || all(df_test_perm$T == 1)) return(NA_real_)
    
    df_test_perm$e <- estimate_e(df_train_perm, df_test_perm, covariate_name, num_cols)
    
    get_tau_hat(df_test_perm)
  })
  
  # -----------------------------
  # 7. Summary statistics
  # -----------------------------
  null_taus_valid <- null_taus[!is.na(null_taus)]
  mean_null <- mean(null_taus_valid)
  sd_null   <- sd(null_taus_valid)
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

hajek_repeat_measures <- function(formula,
                                                                          data,
                                                                          bugs,
                                                                          blocking_variable = "subject",
                                                                          covariate_name,
                                                                          permutations   = 999,
                                                                          split_ratio    = 0.5,
                                                                          numerical_metadata = NULL,
                                                                          method         = "bray") {
  # --- 0. Checks ---
  if (missing(covariate_name))
    stop("Please provide covariate_name (the binary treatment/covariate).")
  
  # --- 1. Parse formula and (optionally) extract distance object on LHS ---
  YVAR <- formula[[2]]
  lhs  <- eval(YVAR, environment(formula), globalenv())
  environment(formula) <- environment()
  D <- lhs  # not used in the core, but kept for symmetry
  
  # --- 2. Check rownames agreement: data vs bugs (and vs D if dist) ---
  if (!all(rownames(data) == rownames(bugs)))
    stop("Row names of 'data' and 'bugs' must match (same samples, same order).")
  
  if (inherits(D, "dist")) {
    if (!all(rownames(as.matrix(D)) == rownames(data)))
      stop("Row names of 'dist' object and 'data' must match.")
  }
  
  # --- 3. Subset to blocking variable + RHS terms ---
  rhs_terms <- labels(terms(formula))
  data_sub  <- data[, c(blocking_variable, rhs_terms), drop = FALSE]
  
  if (any(is.na(data_sub)))
    stop("data must not have NA values for Hájek (or pre-filter before calling).")
  
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
  
  # Optional: handle edge cases where one of them is empty
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
  
  # --- 9. Call the core Hájek engine ---
  if (is.null(numerical_metadata)) {
    # Use the old, constant‑e engine
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
    # Use logistic-e core only when you actually want it
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
  
  # Optional heading
  heading <- sprintf(
    paste("Hájek estimator with subject-level split",
          "Blocked by %s",
          "Number of permutations: %d\n", sep = "\n"),
    blocking_variable, permutations
  )
  attr(res, "heading") <- paste0(heading, paste0(deparse(sys.call()), collapse = "\n"))
  
  res
}

# ---------- Hájek with Location and Dispersion (subject-respecting split; no GUI) ----------
hajek_repeat_measures_loc_and_disp <- function(D, 
                                               permute_within, 
                                               blocks = NULL,
                                               block_data, 
                                               permutations = 0,         # we only need observed loc/disp
                                               split_ratio = 0.5,
                                               numerical_metadata,
                                               bugs, 
                                               method = "bray",
                                               covariate_name) {
  library(dplyr)
  
  # Combine permute_within and block_data (align on blocks)
  mtdat <- cbind(permute_within, block_data[blocks,,drop=FALSE])
  
  df <- mtdat %>%
    rownames_to_column("id") %>%
    full_join(bugs %>% rownames_to_column("id"), by = "id") %>%
    dplyr::select(-matches("\\.x$"), -matches("\\.y$"))
  
  # Factorize covariate
  mtdat[[covariate_name]] <- factor(mtdat[[covariate_name]])
  cov_levels <- levels(mtdat[[covariate_name]])
  if (length(cov_levels) != 2) stop("Binary covariate required.")
  
  # Subject-level split (train/test by subject; both levels in train)
  if (is.null(blocks)) stop("Provide 'blocks' (subject ids per sample).")
  subject_ids <- unique(blocks)
  
  train_contains_both_levels <- function(subject_train) {
    train_rows <- which(blocks %in% subject_train)
    lv <- unique(mtdat[[covariate_name]][train_rows])
    length(lv) == 2
  }
  
###
  # --- Balanced subject-level split requiring BOTH levels in train AND test ---
  
  max_tries <- 2000L
  subject_ids <- unique(blocks)
  n_train_subjects <- max(1L, floor(length(subject_ids) * split_ratio))
  
  found_split <- FALSE
  
  for (try in seq_len(max_tries)) {
    # 1) Sample train subjects; derive row indices
    train_subjects <- sample(subject_ids, n_train_subjects)
    train_indices  <- which(blocks %in% train_subjects)
    test_indices   <- setdiff(seq_len(nrow(df)), train_indices)
    
    # 2) Basic size sanity checks
    if (length(train_indices) < 2L || length(test_indices) < 2L) next
    
    # 3) Compute treatment indicators from metadata for train/test
    T_train <- ifelse(mtdat[[covariate_name]][train_indices] == cov_levels[1], 0L, 1L)
    T_test  <- ifelse(mtdat[[covariate_name]][test_indices]  == cov_levels[1], 0L, 1L)
    
    # 4) Require BOTH levels in BOTH splits
    if (any(T_train == 0L) && any(T_train == 1L) &&
        any(T_test  == 0L) && any(T_test  == 1L)) {
      
      # 5) Build train/test data frames
      df_train <- df[train_indices, , drop = FALSE]
      df_test  <- df[test_indices,  , drop = FALSE]
      
      # 6) Group-1 centroid from TRAIN only
      group1_centroid <- df_train %>%
        dplyr::filter(.data[[covariate_name]] == cov_levels[1]) %>%
        dplyr::summarise(dplyr::across(colnames(bugs), mean)) %>%
        as.numeric()
      
      # 7) Outcome Y on TEST
      if (method == "euclidean") {
        df_test$Y <- apply(df_test[, colnames(bugs), drop = FALSE], 1, function(row) {
          sum((row - group1_centroid)^2)
        })
      } else if (method == "bray") {
        df_test$Y <- apply(df_test[, colnames(bugs), drop = FALSE], 1, function(row) {
          den <- sum(row + group1_centroid); if (den == 0) return(0)
          sum(abs(row - group1_centroid)) / den
        })
      } else {
        stop("Unsupported method.")
      }
      
      # 8) Finalize treatment on TEST (reuse computed T_test)
      df_test$T <- T_test
      
      found_split <- TRUE
      break
    }
  }
  
  if (!found_split) {
    stop("Could not find a subject-level split with both treatment levels in ",
         "both train and test after max_tries. Consider adjusting split_ratio or data.")
  }
###
  # Constant propensity score in test
  df_test$e <- mean(df_test$T)
  
  # Hájek + location + dispersion
  get_tau_hat <- function(df_sub, c_ref) {
    # Hájek on Y (τ̂)
    tau_hat <- df_sub %>%
      summarise(
        treat_numerator     = sum(T * Y / e),
        treat_denominator   = sum(T / e),
        control_numerator   = sum((1 - T) * Y / (1 - e)),
        control_denominator = sum((1 - T) / (1 - e))
      ) %>%
      mutate(tau_hat = (treat_numerator / treat_denominator) -
               (control_numerator / control_denominator)) %>%
      pull(tau_hat)
    
    # Location & dispersion (feature space)
    bugs_cols <- colnames(bugs)
    x0 <- as.matrix(df_sub[df_sub$T == 0, bugs_cols, drop = FALSE])
    x1 <- as.matrix(df_sub[df_sub$T == 1, bugs_cols, drop = FALSE])
    
    if (nrow(x0) < 2 || nrow(x1) < 2) {
      tau_loc  <- NA_real_
      tau_disp <- NA_real_
    } else {
      mu0 <- colMeans(x0)
      mu1 <- colMeans(x1)
      tr_Sigma0 <- sum(apply(x0, 2, stats::var))
      tr_Sigma1 <- sum(apply(x1, 2, stats::var))
      tau_loc   <- sum((mu1 - c_ref)^2) - sum((mu0 - c_ref)^2)
      tau_disp <- tr_Sigma1 - tr_Sigma0
    }
    list(tau_hat = tau_hat, tau_loc = tau_loc, tau_disp = tau_disp)
  }
  
  obs <- get_tau_hat(df_test, c_ref = group1_centroid)
  
  list(
    observed_tau  = obs$tau_hat,   # <-- FIXED: total τ̂
    observed_loc  = obs$tau_loc,
    observed_disp = obs$tau_disp
  )
}
