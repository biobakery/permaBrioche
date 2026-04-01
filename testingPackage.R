rm(list = ls())
#library(permute)
library(permabrioche)

set.seed(1)
n_subj <- 5; reps <- 2; N <- n_subj * reps
subj   <- factor(rep(paste0("S", 1:n_subj), each = reps))
X      <- rnorm(N)
dat    <- data.frame(subject = subj, X = X, row.names = paste0("id", 1:N))
mat    <- matrix(rnorm(N * 10), nrow = N, dimnames = list(rownames(dat), NULL))
D      <- vegdist(mat, method = "bray")

permabrioche::PERMANOVA_repeat_measures(
  D ~ X,
  data = dat,
  blocking_variable = "subject",
  permutations = 9
  )
