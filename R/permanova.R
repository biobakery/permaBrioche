#' PERMANOVA for Repeated Measures
#'
#' A convenience wrapper to run a block-aware permutation test for adonis
#' under a reduced model, preserving subject/cluster structure. Internally uses
#' \code{vegan::adonis2} and \code{permute} to create within-block
#' permutations and between-block shuffles of static (block-level) covariates.
#' All terms on the right-hand side of the formula are tested jointly as a
#' single omnibus test.
#'
#' @param formula A model formula of the form \code{D ~ X1 + X2 + ...}, where
#'   \code{D} is a \code{\link{dist}} object compatible with \code{vegan::adonis2}.
#' @param data A \code{data.frame} with rownames matching the distance object and
#'   columns for the RHS formula terms and the \code{blocking_variable}.
#' @param sample_id Optional character scalar giving the column in \code{data}
#'   that contains sample identifiers matching the distance matrix labels.
#'   If \code{NULL}, rownames(data) are used.
#' @param blocking_variable Character scalar; column in \code{data} defining the
#'   subject/cluster for repeated measures (default: \code{"subject"}).
#' @param permutations Integer; number of permutations for the test (default 999).
#' @param na.rm Logical; if \code{TRUE}, samples with any missing metadata are dropped
#'   in a block-aware manner; otherwise an error is thrown if any NA are present.
#' @param center_R2 Logical; if \code{TRUE}, subtracts the mean of the
#'   permutation null R\eqn{^2} distribution from the observed R\eqn{^2},
#'   returning a \code{R2_centered} column.
#'
#' @return A \code{vegan::adonis2} table with one combined model row (omnibus
#'   test of all RHS terms jointly), followed by \code{Residual} and
#'   \code{Total} rows. The \code{Pr(>F)} column is computed from the
#'   block-aware permutation null rather than \code{adonis2}'s internal
#'   permutation. An optional \code{na.removed} attribute is included if
#'   \code{na.rm = TRUE}.
#'   If \code{center_R2 = TRUE}, an additional \code{R2_centered} column is
#'   included, and the mean null R\eqn{^2} values are stored in the
#'   \code{"null_means_R2"} attribute.
#'
#' @details
#' This function separates metadata into two groups:
#' variables that vary within blocks (subjects) and variables that are static
#' across all samples of the same block. During permutation, the within-block
#' variables are permuted respecting the block structure, and the static
#' block-level variables are shuffled across block identities and then reassigned
#' to samples by their block membership. This yields a permutation that respects
#' the repeated-measures structure.
#'
#' All right-hand-side terms are tested jointly as a single omnibus null
#' hypothesis (i.e., no covariate is associated with the outcome). This test
#' is exact under between- and within-subject exchangeability; see Appendix
#' for a formal proof. Per-term sequential or marginal tests are not supported
#' by this function.
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
#'
#' # Variant covariates: one value per observation, varying within subjects
#' X_var <- rnorm(N)
#' Z_var <- rnorm(N)
#'
#' # Invariant covariates: one value per subject, constant within subjects
#' X_inv <- rep(rnorm(n_subj), each = reps)
#' Z_inv <- rep(rnorm(n_subj), each = reps)
#'
#' meta <- data.frame(
#'   subject = subject,
#'   X_var   = X_var,
#'   Z_var   = Z_var,
#'   X_inv   = X_inv,
#'   Z_inv   = Z_inv,
#'   row.names = paste0("id", seq_len(N))
#' )
#'
#' mat <- matrix(rnorm(N * 5), nrow = N,
#'               dimnames = list(rownames(meta), NULL))
#'
#' D <- vegan::vegdist(mat, method = "bray")
#'
#' ## omnibus test (single variant covariate)
#' PERMANOVA_repeat_measures(
#'   D ~ X_var,
#'   data = meta,
#'   blocking_variable = "subject",
#'   permutations = 49
#' )
#'
#' ## omnibus test (single invariant covariate)
#' PERMANOVA_repeat_measures(
#'   D ~ X_inv,
#'   data = meta,
#'   blocking_variable = "subject",
#'   permutations = 49
#' )
#'
#' ## omnibus test (multiple variant covariates -- same type, exact test)
#' PERMANOVA_repeat_measures(
#'   D ~ X_var + Z_var,
#'   data = meta,
#'   blocking_variable = "subject",
#'   permutations = 49
#' )
#'
#' ## omnibus test (multiple invariant covariates -- same type, exact test)
#' PERMANOVA_repeat_measures(
#'   D ~ X_inv + Z_inv,
#'   data = meta,
#'   blocking_variable = "subject",
#'   permutations = 49
#' )
#'
#' \donttest{
#' ## ---- larger example ----
#' PERMANOVA_repeat_measures(
#'   D ~ X_var,
#'   data = meta,
#'   blocking_variable = "subject",
#'   permutations = 99
#' )
#' }
#'
#' @seealso \code{\link[vegan]{adonis2}}, \code{\link[permute]{how}}
#'
#' @import vegan
#' @import permute
#' @importFrom stats as.formula aggregate binomial glm predict terms var
#' @importFrom dplyr full_join filter summarise across mutate pull
#' @importFrom tibble rownames_to_column
#' @export
PERMANOVA_repeat_measures <- function(formula,
                                      data,
                                      sample_id = NULL, 
                                      blocking_variable = "subject",
                                      permutations = 999,
                                      na.rm = FALSE,
                                      center_R2 = FALSE) {
  by = NULL
  
  if (!is.null(by) && !by %in% c("terms", "margin"))
    stop('`by` must be NULL, "terms", or "margin"')
  
  data <- as.data.frame(data)
  
  # --- 1. Parse formula and extract distance object ---
  YVAR <- formula[[2]]
  lhs  <- eval(YVAR, environment(formula), globalenv())
  environment(formula) <- environment()
  if (!inherits(lhs, "dist"))
    stop("lhs of formula must be an adonis2-compatible 'dist' object")

  D <- lhs

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

    ## Disallow mixing static and varying covariates on the RHS
  if (any(static_vars) && any(!static_vars)) {
    static_names  <- names(static_vars)[static_vars]
    varying_names <- names(static_vars)[!static_vars]

    stop(
      paste0(
        "All right-hand-side covariates must be of the same type with ",
        "respect to the blocking variable '", blocking_variable, "'.\n",
        "Currently detected:\n",
        "  - Block-invariant (static within blocks): ",
        paste(static_names, collapse = ", "), "\n",
        "  - Within-block varying: ",
        paste(varying_names, collapse = ", "), "\n\n",
        "You are requesting a single omnibus test that combines both types, ",
        "which this implementation does not support. Please fit separate ",
        "models using only static or only within-block varying covariates."
      ),
      call. = FALSE
    )
  }

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

  # --- 8. Prepare metadata_order for the core engine ---
  # Between-block (static) terms are tested before within-block-varying terms.
  metadata_order <- c(colnames(block_data), colnames(permute_within))

  # --- 9. Call the core engine ---
  res <- PERMANOVA_repeat_measures_core(
    D              = D,
    permute_within = permute_within,
    blocks         = blocks,
    block_data     = block_data,
    permutations   = permutations,
    metadata_order = metadata_order,
    na.rm          = na.rm,
    center_R2      = center_R2,
    by             = by
  )
  
  by_label <- switch(
    if (is.null(by)) "NULL" else by,
    "NULL" = "Overall (omnibus) test of all terms jointly",
    "terms"  = "Terms added sequentially (first to last)",
    "margin" = "Marginal effects of terms (each adjusted for all others)"
  )

  heading <- sprintf(paste(
    "Permutation test for adonis under reduced model",
    "%s",
    "Permutation: blocked by %s",
    "Number of permutations: %d",
    sep = "\n"),
    by_label, blocking_variable, permutations
  )
  attr(res, "heading") <- paste0(heading, "\n", paste0(deparse(sys.call()), collapse = "\n"))

  res
}

#' @keywords internal
#' @noRd
#' @importFrom vegan adonis2
#' @importFrom permute shuffle how
PERMANOVA_repeat_measures_core <- function(
    D, permute_within, blocks = NULL, block_data,
    permutations = 999,
    metadata_order = c(names(block_data), names(permute_within)),
    na.rm = FALSE,
    center_R2 = FALSE,
    by = NULL) {

  if (!is.null(by) && !by %in% c("terms", "margin"))
    stop('`by` must be NULL, "terms", or "margin"')

  if (!inherits(D, "dist")) stop("D must be a dist object")

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

  ord <- attr(D, "Labels") #ord <- rownames(as.matrix(D))
  if (length(ord) != nrow(permute_within) || length(blocks) != length(ord))
    stop("blocks, permute_within, and D are not the same size")
  
  if (!is.null(rownames(permute_within))) {
    if (!all(ord %in% rownames(permute_within)))
      stop("Some samples in D are missing from permute_within")
    
    permute_within <- permute_within[ord, , drop = FALSE]
  }
  
  if (any(is.na(blocks))) stop("NAs are not allowed in blocks")

  if (is.factor(blocks)) {
    if (any(!(levels(blocks) %in% rownames(block_data))))
      stop("not all block levels are contained in block_data")
    block_data <- block_data[match(levels(blocks), rownames(block_data)), , drop = FALSE]
    blocks     <- as.numeric(blocks)
  } else if (is.numeric(blocks)) {
    if (any(blocks < 1) || max(blocks) > nrow(block_data))
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
      D      <- stats::as.dist(D[keep, keep])

      if (length(blocks) < ncol(permute_within) + ncol(block_data)) {
        stop(sprintf("After omitting samples, samples (%d) < metadata (%d)",
                     length(blocks), ncol(permute_within) + ncol(block_data)))
      }
      na.removed <- n_prerm - length(blocks)
    } else {
      stop("Some metadata is NA! adonis does not support any NA in the metadata")
    }
  }

  # --- Warn about singleton blocks ---
  # A block with only one observation has only one possible within-block
  # arrangement (itself), so it contributes no permutation variability to
  # within-block (permute_within) terms' null distribution.
  if (ncol(permute_within) > 0L) {
    block_sizes <- table(blocks)
    n_singletons <- sum(block_sizes == 1)
    if (n_singletons / length(block_sizes) > 0.5) {
      warning(sprintf(
        paste0(
          "%d of %d block(s) (%.0f%%) contain only a single observation. ",
          "Singleton blocks contribute no permutation variability to ",
          "within-block term(s) (%s); their within-block covariate values ",
          "are never actually shuffled. With a majority of blocks being ",
          "singletons, null distributions for these terms are driven by a ",
          "small minority of blocks with >= 2 observations, which may make ",
          "p-values conservative or coarse-grained."
        ),
        n_singletons, length(block_sizes), 100 * n_singletons / length(block_sizes),
        paste(colnames(permute_within), collapse = ", ")
      ))
    }
  }

  # choose permutations for the adonis2 fit used to get R2:
  # by = "margin" needs >= 1 to avoid an internal vegan error
  perm_for_fit <- if (identical(by, "margin")) 1L else 0L

  mtdat <- cbind(permute_within, block_data[blocks,,drop=FALSE])
  ad    <- vegan::adonis2(D ~ ., permutations = perm_for_fit, by = by, data = mtdat[, metadata_order, drop=FALSE])
  R2    <- ad$R2; names(R2) <- rownames(ad)

  nullsamples <- matrix(NA_real_, nrow = length(R2), ncol = permutations)
  ctrl <- how(blocks = blocks)
  for (i in seq_len(permutations)) {
    within.i <- shuffle(nrow(permute_within), control = ctrl)
    block.i  <- sample(seq_len(nrow(block_data)))
    mtdat.i  <- cbind(
      permute_within[within.i,,drop=FALSE],
      block_data[block.i,,drop=FALSE][blocks,,drop=FALSE]
    )
    perm.ad <- vegan::adonis2(D ~ ., permutations = perm_for_fit, by = by, data = mtdat.i[, metadata_order, drop=FALSE])
    nullsamples[,i] <- perm.ad$R2
  }

  n <- length(R2)
  stopifnot(identical(rownames(ad)[c(n - 1, n)], c("Residual", "Total")))
  R2[n-1]           <- 1 - R2[n-1]
  nullsamples[n-1,] <- 1 - nullsamples[n-1,]
  null_means <- rowMeans(nullsamples, na.rm = TRUE)
  if (center_R2) {
    R2_centered <- R2 - null_means
    R2_centered[c(n - 1, n)] <- NA_real_
    ad$R2_centered <- R2_centered
  }
  attr(ad, "null_means_R2") <- null_means
  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)
  P[c(n - 1, n)] <- NA_real_
  ad$`Pr(>F)` <- P
  if (na.rm) ad$na.removed <- na.removed
  ad
}
