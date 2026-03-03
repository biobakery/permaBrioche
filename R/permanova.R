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
#' @param blocking_variable Character scalar; column in \code{data} defining the
#'   subject/cluster for repeated measures (default: \code{"subject"}).
#' @param permutations Integer; number of permutations for the test (default 999).
#' @param na.rm Logical; if \code{TRUE}, samples with any missing metadata are dropped
#'   in a block-aware manner; otherwise an error is thrown if any NA are present.
#'
#' @return A \code{vegan::adonis} table (from \code{adonis2}) with an added
#'   \code{Pr(>F)} column computed from the custom null and an optional
#'   \code{na.removed} attribute if \code{na.rm=TRUE}. A descriptive \code{"heading"}
#'   attribute is also added.
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
#' @importFrom stats as.formula
#' @importFrom vegan adonis2
#' @export
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

  # Edge cases
  if (ncol(permute_within) == 0L) {
    permute_within <- data.frame(row.names = rownames(data_sub))
  }
  if (ncol(block_data_full) == 0L) {
    block_data_full <- as.data.frame(matrix(0, nrow = 1, ncol = 0))
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
    na.rm          = na.rm
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
    na.rm = FALSE) {

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
  for (i in seq_len(permutations)) {
    within.i <- permute::shuffle(nrow(permute_within), control = permute::how(blocks=blocks))
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

  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)
  P[n] <- NA_real_
  ad$`Pr(>F)` <- P
  if (na.rm) ad$na.removed <- na.removed
  ad
}
