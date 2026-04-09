#' PERMANOVA for Repeated Measures (Block-Aware, macOS-safe)
#'
#' A convenience wrapper to run a block-aware permutation test for adonis
#' under a reduced model, preserving subject/cluster structure. Internally uses
#' \code{vegan::adonis2} and \code{permute} to create within-block
#' permutations and between-block shuffles of static (block-level) covariates.
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
#'   permutation null R\eqn{^2} distribution from the observed R\eqn{^2}
#'   for each term, returning a \code{R2_centered} column.
#'
#' @return A \code{vegan::adonis} table (from \code{adonis2}) with an added
#'   \code{Pr(>F)} column computed from the custom null and an optional
#'   \code{na.removed} attribute if \code{na.rm=TRUE}. A descriptive \code{"heading"}
#'   attribute is also added.
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
#' to samples by their block membership. This yields a reduced-model style
#' permutation respecting repeated measures.
#'
#' @examples
#' \dontrun{
#'   library(vegan)
#'   set.seed(1)
#'   # Toy example
#'   n_subj <- 10
#'   reps   <- 3
#'   N      <- n_subj * reps
#'   subj   <- factor(rep(paste0("S", 1:n_subj), each = reps))
#'   X      <- rnorm(N)
#'   Z      <- rep(rnorm(n_subj), each = reps) # static per subject
#'   data   <- data.frame(subject = subj, X = X, Z = Z, row.names = paste0("id", 1:N))
#'
#'   mat <- matrix(rnorm(N * 20), nrow = N, dimnames = list(rownames(data), NULL))
#'   D   <- vegdist(mat, method = "bray")
#'
#'   # Test X after accounting for Z, blocking by subject
#'   PERMANOVA_repeat_measures(D ~ X + Z, data = data,
#'                             blocking_variable = "subject", permutations = 99)
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
  metadata_order <- c(colnames(permute_within), colnames(block_data))

  # --- 9. Call the core engine ---
  res <- PERMANOVA_repeat_measures_core(
    D              = D,
    permute_within = permute_within,
    blocks         = blocks,
    block_data     = block_data,
    permutations   = permutations,
    metadata_order = metadata_order,
    na.rm          = na.rm,
    center_R2      = center_R2
  )

  heading <- sprintf(paste(
    "Permutation test for adonis under reduced model",
    "Terms added sequentially (first to last)",
    "Permutation: blocked by %s",
    "Number of permutations: %d",
    sep = "\n"),
    blocking_variable, permutations
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
    metadata_order = c(names(permute_within), names(block_data)),
    na.rm = FALSE,
    center_R2 = FALSE) {

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

  mtdat <- cbind(permute_within, block_data[blocks,,drop=FALSE])
  ad    <- vegan::adonis2(D ~ ., permutations = 0, data = mtdat[, metadata_order, drop=FALSE])
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
    perm.ad <- vegan::adonis2(D ~ ., permutations = 0, data = mtdat.i[, metadata_order, drop=FALSE])
    nullsamples[,i] <- perm.ad$R2
  }

  n <- length(R2)
  R2[n-1]           <- 1 - R2[n-1]
  nullsamples[n-1,] <- 1 - nullsamples[n-1,]

  null_means <- rowMeans(nullsamples, na.rm = TRUE)
  if (center_R2) {
    ad$R2_centered <- R2 - null_means
  }
  attr(ad, "null_means_R2") <- null_means

  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)
  P[n] <- NA_real_
  ad$`Pr(>F)` <- P
  if (na.rm) ad$na.removed <- na.removed
  ad
}
