#' PERMANOVA for repeated measures (core engine)
#'
#' **Internal**. Performs a blocked (by subject) permutation procedure for an
#' adonis2 model where some covariates vary within subject and others are
#' subject-level constants.
#'
#' @param D A `dist` object (required).
#' @param permute_within data.frame of covariates that vary within blocks.
#' @param blocks Factor / numeric / character vector indicating block membership
#'   (e.g., subject ID) for each sample. Must align to `D`.
#' @param block_data data.frame of block-level covariates (one row per block level).
#' @param permutations Number of permutations (default 999).
#' @param metadata_order Character vector giving the column order to fit in adonis2.
#' @param na.rm Logical; if TRUE drop any samples/blocks containing NA metadata.
#'
#' @return A data.frame like `vegan::adonis2` output, with permuted p-values in
#'   `Pr(>F)`. Attribute `"heading"` may be attached by caller.
#' @keywords internal
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

  nullsamples <- matrix(NA, nrow = length(R2), ncol = permutations)
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
  P[n] <- NA
  ad$`Pr(>F)` <- P
  if (na.rm) ad$na.removed <- na.removed
  ad
}

#' PERMANOVA for repeated measures (formula interface)
#'
#' @param formula A formula with a `dist` object on the LHS, e.g. `D ~ cov1 + cov2`.
#' @param data Sample metadata (rows = samples) with rownames matching the
#'   `dist` object.
#' @param blocking_variable Character name of the subject/block variable (default "subject").
#' @param permutations Number of permutations (default 999).
#' @param na.rm If `TRUE`, remove samples/blocks with NA metadata subject-wise.
#'
#' @return An object like `vegan::adonis2` result with permuted p-values.
#' @export
#' @examples
#' ## Pseudo-example (not run):
#' # D <- vegan::vegdist(matrix(rpois(100, 5), nrow=10), method = "bray")
#' # md <- data.frame(subject = rep(letters[1:5], each=2),
#' #                  treatment = rep(c("A","B"), 5),
#' #                  row.names = attr(D, "Labels"))
#' # PERMANOVA_repeat_measures(D ~ treatment, data = md,
#' #                           blocking_variable = "subject", permutations = 99)
PERMANOVA_repeat_measures <- function(formula,
                                      data,
                                      blocking_variable = "subject",
                                      permutations = 999,
                                      na.rm = FALSE) {

  YVAR <- formula[[2]]
  lhs  <- eval(YVAR, environment(formula), globalenv())
  environment(formula) <- environment()
  if (!inherits(lhs, "dist"))
    stop("lhs of formula must be an adonis2-compatible 'dist' object")
  D <- lhs

  if (!all(rownames(as.matrix(D)) == rownames(data)))
    stop("Row names of distance matrix must match row names of data")

  rhs_terms <- labels(terms(formula))
  data_sub  <- data[, c(blocking_variable, rhs_terms), drop = FALSE]

  if (any(is.na(data_sub)) && !na.rm)
    stop("data must not have NA values (or set na.rm = TRUE)")

  blocks <- as.factor(data_sub[[blocking_variable]])
  rhs_only <- data_sub[, rhs_terms, drop = FALSE]

  agg_res <- stats::aggregate(
    rhs_only,
    list(block = data_sub[[blocking_variable]]),
    function(x) length(unique(x)) == 1
  )
  static_vars <- sapply(agg_res[, -1, drop = FALSE], all)

  permute_within  <- rhs_only[, names(static_vars)[!static_vars], drop = FALSE]
  block_data_full <- rhs_only[, names(static_vars)[static_vars],  drop = FALSE]

  if (ncol(permute_within) == 0L) {
    permute_within <- data.frame(row.names = rownames(data_sub))
  }
  if (ncol(block_data_full) == 0L) {
    block_data_full <- as.data.frame(matrix(0, nrow = 1, ncol = 0))
    rownames(block_data_full) <- levels(blocks)
  }

  block_data <- block_data_full[!duplicated(blocks), , drop = FALSE]
  rownames(block_data) <- levels(blocks)
  metadata_order <- c(colnames(permute_within), colnames(block_data))

  res <- PERMANOVA_repeat_measures_core(
    D              = D,
    permute_within = permute_within,
    blocks         = blocks,
    block_data     = block_data,
    permutations   = permutations,
    metadata_order = metadata_order,
    na.rm          = na.rm
  )

  heading <- sprintf(paste("Permutation test for adonis under reduced model",
                           "Terms added sequentially (first to last)",
                           "Permutation: blocked by %s",
                           "Number of permutations: %d",
                           sep = "\n"),
                     blocking_variable, permutations)
  attr(res, "heading") <- paste0(heading, "\n", paste0(deparse(sys.call()), collapse = "\n"))
  res
}
