
# ======================================================================
# PERMANOVA Simulations (macOS/RStudio-safe: PSOCK clusters)
# Invariant + Variant branches, with G1 and G2 implemented
# ======================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# Libraries
suppressPackageStartupMessages({
  required_packages <- c(
    "vegan", "tidyr", "ggplot2", "reshape2",
    "cowplot", "ggforce", "doParallel", "permute",
    "tibble", "dplyr"
  )
  
  missing_packages <- required_packages[
    !sapply(required_packages, requireNamespace, quietly = TRUE)
  ]
  
  if (length(missing_packages) > 0) {
    install.packages(missing_packages, repos = "https://cloud.r-project.org")
  }
  
  library(vegan)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(cowplot)
  library(ggforce)
  library(doParallel)
  library(permute)
  library(tibble)
  library(dplyr)
})

# Paths
my_dir <- "/Users/sithijamanage/HuttenhowerLabSummer2025/"
setwd(my_dir)

# If your helper functions are in separate files, source them here:
# source(paste0(my_dir,"hmp_omnibus_tests.R"))
source(paste0(my_dir,"permanova_spiked_R_file_iterations/repeat_measures_userfriendly.R"))

# Data
load(paste0(my_dir, "d_hmp.Rda")) # loads d
load(paste0(my_dir, "bugs_hmp.Rda"))      # loads bugs

numerical_metadata <- d[, sapply(d, is.numeric)]
numerical_metadata <- numerical_metadata[, colSums(is.na(numerical_metadata)) == 0]

d_orig    <- d
bugs      <- bugs[,-1]
bugs_orig <- bugs
rm(d, bugs)

# Parameters ------------------------------------------------------------
n_iterations      <- 10
variant_indicator <- 0   # <-- set 1 for VARIANT, 0 for INVARIANT
nCores            <- max(1, parallel::detectCores() - 1)

reorgSamples   <- c(FALSE, rep(TRUE, 4))
meanSamples    <- c(NA, 2, 4, 2, 4)
distribSamples <- c(NA, "even", "even", "uneven", "uneven")
param          <- data.frame(n_iterations, reorgSamples, meanSamples, distribSamples)

seed_multiplier <- 1000
prob            <- c(0, 0.1, 0.25, 0.5)  # used as flip/strength parameter

# Output dirs -----------------------------------------------------------
dir_output <- paste0(my_dir, "results/")
dir.create(dir_output, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir_output, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(dir_output, "data"),  recursive = TRUE, showWarnings = FALSE)

# Theme -----------------------------------------------------------------
my_theme <- theme(
  title            = element_text(family = "Helvetica", face = "bold", size = 12),
  axis.text.x      = element_text(family = "Helvetica", size = 10),
  axis.text.y      = element_text(family = "Helvetica", size = 10),
  legend.text      = element_text(family = "Helvetica", size = 10),
  legend.title     = element_text(family = "Helvetica", size = 12),
  axis.title.x     = element_text(family = "Helvetica", size = 12),
  axis.title.y     = element_text(family = "Helvetica", size = 12),
  strip.text       = element_text(family = "Helvetica", face = "bold", size = 10),
  axis.line        = element_line(colour = 'black', size = 0.5),
  axis.ticks       = element_line(colour = "black", size = 0.5),
  legend.position  = "none",
  strip.background = element_rect(fill = "lightblue", colour = "black", size = 1)
)

remove_tail <- function(x, sep = "_", del) {
  sapply(strsplit(x, split = sep, fixed = TRUE),
         function(i) paste(head(i, -del), collapse = sep))
}

# ---------- Resampling helpers (no View calls) ----------
adjust_sample_number <- function(df_meta, df_bugs, subj, des_num) {
  row.names(df_meta) <- NULL
  df_meta$sample     <- as.character(df_meta$sample)
  subj_rows          <- which(df_meta$subject == subj)
  cur_num            <- length(subj_rows)
  
  if (des_num < cur_num) {
    drop_subj_rows <- sample(x = subj_rows, size = (cur_num - des_num), replace = FALSE)
    df_meta <- df_meta[-drop_subj_rows, ]
    df_bugs <- df_bugs[-drop_subj_rows, ]
  } else if (des_num > cur_num) {
    add_subj_rows <- if (length(subj_rows) == 1) rep(subj_rows, des_num - 1)
    else sample(x = subj_rows, size = (des_num - cur_num), replace = TRUE)
    addition      <- df_meta[add_subj_rows, ]
    addition$time <- seq(20, 20 + nrow(addition) - 1, 1)
    add_names     <- paste0("Pair", addition$subject, "-", addition$time)
    addition$sample <- add_names
    df_meta       <- rbind(df_meta, addition)
    
    addition_bugs <- df_bugs[add_subj_rows, ]
    row.names(addition_bugs) <- add_names
    df_bugs       <- rbind(df_bugs, addition_bugs)
  }
  list(df_meta, df_bugs)
}

resample_data <- function(df_meta, df_bugs, type, mean) {
  subjects      <- unique(df_meta$subject)
  new_meta_list <- list()
  new_bugs_list <- list()
  
  if (type == "even") {
    target_samples <- round(mean)
    for (subj in subjects) {
      result <- adjust_sample_number(df_meta, df_bugs, subj, target_samples)
      new_meta_list[[subj]] <- result[[1]]
      new_bugs_list[[subj]] <- result[[2]]
    }
  } else if (type == "uneven") {
    n_subjects  <- length(subjects)
    n_high      <- ceiling(n_subjects * 0.3)
    n_low       <- n_subjects - n_high
    low_count   <- 1
    high_count  <- round((mean * n_subjects - n_low * low_count) / n_high)
    high_count  <- max(high_count, 2)
    high_subjects <- sample(subjects, n_high, replace = FALSE)
    
    for (subj in subjects) {
      des_num <- if (subj %in% high_subjects) high_count else low_count
      result  <- adjust_sample_number(df_meta, df_bugs, subj, des_num)
      new_meta_list[[subj]] <- result[[1]]
      new_bugs_list[[subj]] <- result[[2]]
    }
  } else {
    stop("Invalid type parameter. Must be 'even' or 'uneven'.")
  }
  
  d_new     <- do.call(rbind, new_meta_list)
  bugs_new  <- do.call(rbind, new_bugs_list)
  
  d_new$sample      <- make.unique(as.character(d_new$sample))
  row.names(d_new)  <- d_new$sample
  row.names(bugs_new) <- d_new$sample
  
  d_new    <- d_new[order(d_new$sample), ]
  bugs_new <- bugs_new[order(row.names(bugs_new)), ]
  
  stopifnot(nrow(d_new) == nrow(bugs_new))
  stopifnot(identical(row.names(d_new), row.names(bugs_new)))
  list(d_new, bugs_new)
}

# ---------- Invariant simulation kernel (with G1/G2) ----------
run_sim_invariable <- function(seed, prob, d, bugs_mat, dist_bugs) {
  set.seed(seed)
  
  # Collapse by subject
  bugs_collapsed <- aggregate(bugs_mat, by = list(subject = d$subject), FUN = mean)
  row.names(bugs_collapsed) <- bugs_collapsed$subject
  bugs_collapsed$subject    <- NULL
  
  subjects <- row.names(bugs_collapsed)
  
  kclust <- kmeans(bugs_collapsed, centers = 2, nstart = 50)
  cluster_assignment <- ifelse(kclust$cluster == which.min(table(kclust$cluster)), 0, 1)
  
  # Flip some subjects at random
  set.seed(seed + 100)
  flip_subjects <- which(runif(length(cluster_assignment)) < prob)
  cluster_assignment[flip_subjects] <- 1 - cluster_assignment[flip_subjects]
  set.seed(seed)
  
  df_subject <- data.frame(subject = subjects,
                           nonsense_invariant = cluster_assignment)
  df_subject$nonsense_invariant_f <- ifelse(df_subject$nonsense_invariant == 0, "a", "b")
  row.names(df_subject) <- df_subject$subject
  
  df_sample <- merge(x = d[, c("sample", "subject")], y = df_subject,
                     by = "subject", all.x = TRUE)
  row.names(df_sample) <- df_sample$sample
  
  # Unadjusted
  perm_unadj <- adonis2(bugs_mat ~ nonsense_invariant_f, data = df_sample,
                        method = "bray", permutations = 999)
  
  # -------- NEW: within‑block PERMANOVA via formula interface --------
  # data frame for PERMANOVA_repeat_measures:
  #   - must contain blocking_variable ("subject")
  #   - and the covariate(s) on the RHS
  data_perm <- df_sample[, c("subject", "nonsense_invariant_f")]
  rownames(data_perm) <- rownames(df_sample)   # ensure match to dist_bugs
  
  #***
  #labs <- attr(dist.bugs, "Labels")
  #if(setequal(rownames(df_sample), labs)!=TRUE)# TRUE on both datasets
  #  stop("setequal(rownames(df_sample), labs)!=TRUE")# TRUE on Pair…, FALSE on HMP
  #if(identical(rownames(df_sample), labs)!=TRUE)
  #  stop("identical(rownames(df_sample), labs)!=TRUE")
  #if(identical(rownames(data_perm), labs) !=TRUE)
  #  stop("identical(rownames(data_perm), labs) !=TRUE")
  #stop("End")
  #***
  
  perm_adj_B <- PERMANOVA_repeat_measures(
    formula           = dist_bugs ~ nonsense_invariant_f,
    data              = data_perm,
    blocking_variable = "subject",
    permutations      = 999,
    na.rm             = FALSE
  )
  
  # Hájek (new formula-based wrapper)
  # Build metadata for Hájek: must contain subject and the binary covariate
  data_hajek <- df_sample[, c("subject", "nonsense_invariant_f")]
  rownames(data_hajek) <- rownames(df_sample)   # must match bugs_mat and dist_bugs
  
  hajek <- hajek_repeat_measures(
    formula           = dist_bugs ~ nonsense_invariant_f,
    data              = data_hajek,
    bugs              = bugs_mat,
    blocking_variable = "subject",
    covariate_name    = "nonsense_invariant_f",
    permutations      = 999,
    split_ratio       = 0.5,
    numerical_metadata = NULL,#numerical_metadata,
    method            = "bray"
  )
  
  # v2
  perm_adj_v2 <- adonis2(bugs_mat ~ nonsense_invariant_f + subject,
                         data = df_sample, method = "bray", permutations = 999)
  
  # G1: blocked permutations for adonis2
  h1 <- with(df_sample, permute::how(nperm = 999, blocks = subject))
  perm_adj_G1 <- adonis2(bugs_mat ~ nonsense_invariant_f,
                         data = df_sample, method = "bray", permutations = h1)
  
  # G2: constrained ordination with Condition(subject)
  h1  <- with(df_sample, permute::how(nperm = 999, blocks = subject))
  ord <- rda(bugs_mat ~ nonsense_invariant_f + Condition(subject), data = df_sample)
  perm_adj_G2 <- anova(ord, permutations = h1)  # anova.cca
  
  # Extract G2 stats: first row corresponds to 'nonsense_invariant_f'
  G2_p <- NA_real_; G2_F <- NA_real_
  if (!is.null(perm_adj_G2) && nrow(perm_adj_G2) >= 1) {
    G2_p <- perm_adj_G2$`Pr(>F)`[1]
    G2_F <- perm_adj_G2$F[1]
  }
  
  list(
    unadj_p      = perm_unadj$`Pr(>F)`[1],
    unadj_effect = perm_unadj$F[1],
    adj_B_p      = perm_adj_B$`Pr(>F)`[1],
    adj_B_effect = perm_adj_B$F[1],
    hajek_p      = hajek$pval,
    hajek_effect = hajek$observed,
    adj_v1_p     = 1,  adj_v1_effect = 1,          # not used
    adj_v2_p     = perm_adj_v2$`Pr(>F)`[1],
    adj_v2_effect= perm_adj_v2$F[1],
    adj_G1_p     = perm_adj_G1$`Pr(>F)`[1],
    adj_G1_effect= perm_adj_G1$F[1],
    adj_G2_p     = G2_p,
    adj_G2_effect= G2_F
  )
}

# ---------- Variant simulation kernel (matches your design + G1/G2) ----------
run_sim_variable <- function(seed, prob, d, df_subject, bugs_mat) {
  set.seed(seed)
  
  # -----------------------------
  # 0. Setup
  # -----------------------------
  subject_list <- unique(d$subject)
  df_sample    <- d[, c("sample", "subject")]
  
  # -----------------------------
  # 1. Assign nonsense_variant per subject (biased coin)
  # -----------------------------
  type <- rbinom(length(subject_list), 1, 0.5)  # 0 = type a, 1 = type b
  type <- ifelse(type == 0, 0.8, 1 - 0.8)
  df_sample$nonsense_variant <- NA
  for (i in seq_along(subject_list)) {
    subject_rows <- which(df_sample$subject == subject_list[i])
    df_sample$nonsense_variant[subject_rows] <- rbinom(length(subject_rows), 1, type[i])
  }
  df_sample$nonsense_variant_f <- ifelse(df_sample$nonsense_variant == 0, "a", "b")
  treatment <- df_sample$nonsense_variant
  
  # -----------------------------
  # 2. Spike in microbiome effect (scaled by prob)
  # -----------------------------
  highest_effect <- 50
  effect_size    <- highest_effect - (highest_effect * 2) * prob
  
  bugs_spiked <- bugs_mat
  n_taxa      <- ncol(bugs_spiked)
  
  n_spike_taxa <- max(1L, round(0.1 * n_taxa))
  spike_taxa   <- sample(colnames(bugs_spiked), n_spike_taxa)
  
  for (taxon in spike_taxa) {
    bugs_spiked[treatment == 1, taxon] <-
      bugs_spiked[treatment == 1, taxon] * (1 + effect_size)
  }
  
  # -----------------------------
  # 3. Subject-level confounding
  # -----------------------------
  for (i in seq_along(subject_list)) {
    rows  <- which(d$subject == subject_list[i])
    shift <- matrix(runif(n_taxa, min = 0, max = 0.1),
                    nrow = length(rows), ncol = n_taxa, byrow = TRUE)
    bugs_spiked[rows, ] <- bugs_spiked[rows, ] + shift
  }
  
  # Renormalize
  bugs_spiked <- sweep(bugs_spiked, 1, rowSums(bugs_spiked), FUN = "/")
  bugs_spiked[bugs_spiked < 0] <- 0
  
  # -----------------------------
  # 4. Distances
  # -----------------------------
  dist_bugs_spiked <- vegdist(bugs_spiked, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE)
  
  # -----------------------------
  # 5. Tests
  # -----------------------------
  set.seed(seed)
  perm_unadj <- adonis2(bugs_spiked ~ nonsense_variant_f, data = df_sample,
                        method = "bray", permutations = 999)
  
  # -------- NEW: within‑block PERMANOVA via formula interface --------
  # merge in subject-level my_block_var so it appears as a static covariate
  df_block <- df_subject[, c("subject", "my_block_var")]
  df_data  <- dplyr::left_join(df_sample, df_block, by = "subject")
  rownames(df_data) <- df_data$sample
  
  # data for PERMANOVA_repeat_measures: subject + both covariates
  data_perm <- df_data[, c("subject", "nonsense_variant_f", "my_block_var")]
  rownames(data_perm) <- rownames(df_data)   # should match dist_bugs_spiked
  
  set.seed(seed)
  perm_adj_B <- PERMANOVA_repeat_measures(
    formula           = dist_bugs_spiked ~ nonsense_variant_f + my_block_var,
    data              = data_perm,
    blocking_variable = "subject",
    permutations      = 999,
    na.rm             = FALSE
  )
  
  # Hájek (new formula-based wrapper)
  set.seed(seed)
  
  # df_data already has subject + nonsense_variant_f + my_block_var
  data_hajek <- df_data[, c("subject", "nonsense_variant_f", "my_block_var")]
  rownames(data_hajek) <- rownames(df_data)   # must match bugs_spiked and dist_bugs_spiked
  
  hajek <- hajek_repeat_measures(
    formula           = dist_bugs_spiked ~ nonsense_variant_f + my_block_var,
    data              = data_hajek,
    bugs              = bugs_spiked,
    blocking_variable = "subject",
    covariate_name    = "nonsense_variant_f",
    permutations      = 999,
    split_ratio       = 0.5,
    numerical_metadata = NULL,#numerical_metadata,
    method            = "bray"
  )
  
  set.seed(seed)
  perm_adj_v1 <- adonis2(bugs_spiked ~ subject + nonsense_variant_f, data = df_sample,
                         method = "bray", permutations = 999)
  
  set.seed(seed)
  perm_adj_v2 <- adonis2(bugs_spiked ~ nonsense_variant_f + subject, data = df_sample,
                         method = "bray", permutations = 999)
  
  set.seed(seed)
  h1 <- with(df_sample, permute::how(nperm = 999, blocks = subject))
  perm_adj_G1 <- adonis2(bugs_spiked ~ nonsense_variant_f, data = df_sample,
                         method = "bray", permutations = h1)
  
  set.seed(seed)
  h1  <- with(df_sample, permute::how(nperm = 999, blocks = subject))
  ord <- rda(bugs_spiked ~ nonsense_variant_f + Condition(subject), data = df_sample)
  perm_adj_G2 <- anova(ord, permutations = h1)
  
  G2_p <- NA_real_; G2_F <- NA_real_
  if (!is.null(perm_adj_G2) && nrow(perm_adj_G2) >= 1) {
    G2_p <- perm_adj_G2$`Pr(>F)`[1]
    G2_F <- perm_adj_G2$F[1]
  }
  
  list(
    unadj_p       = perm_unadj$`Pr(>F)`[1],
    unadj_effect  = perm_unadj$F[1],
    adj_B_p       = perm_adj_B$`Pr(>F)`[1],
    adj_B_effect  = perm_adj_B$F[1],
    hajek_p       = hajek$pval,
    hajek_effect  = hajek$observed,
    adj_v1_p      = perm_adj_v1$`Pr(>F)`[2],
    adj_v1_effect = perm_adj_v1$F[2],
    adj_v2_p      = perm_adj_v2$`Pr(>F)`[1],
    adj_v2_effect = perm_adj_v2$F[1],
    adj_G1_p      = perm_adj_G1$`Pr(>F)`[1],
    adj_G1_effect = perm_adj_G1$F[1],
    adj_G2_p      = G2_p,
    adj_G2_effect = G2_F
  )
}

# ---------- Main loop ----------
cat("Checking param object:\n"); print(param); cat("nrow(param):", nrow(param), "\n")

for (j in 1) {  # change to 1:nrow(param) to iterate all parameter sets
  cat("j = ", j, "\n")
  
  if (param$reorgSamples[j] == TRUE) {
    cat("Resampling dataset...\n")
    resampled <- resample_data(df_meta = d_orig, df_bugs = bugs_orig,
                               type = param$distribSamples[j],
                               mean = param$meanSamples[j])
    d    <- resampled[[1]]
    bugs <- resampled[[2]]
  } else {
    cat("Using original dataset.\n")
    d    <- d_orig
    bugs <- bugs_orig
  }
  
  dist.bugs <- vegdist(bugs, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE)
  d$subject <- factor(d$subject)
  
  iterations <- seq_len(n_iterations)
  sim <- tidyr::crossing(iterations, prob)
  sim$seed  <- seq_len(nrow(sim)) * seed_multiplier
  sim$index <- seq_len(nrow(sim))  # explicit index for merges
  
  if (variant_indicator == 0) {
    # ------------------------- INVARIANT -------------------------
    cat("Running Simulation 1 (Invariant)...\n")
    cl <- parallel::makeCluster(nCores, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    
    sim_res <- foreach(i = 1:nrow(sim),
                       .export   = c("run_sim_invariable",
                                     "PERMANOVA_repeat_measures",
                                     "hajek_repeat_measures",
                                     "numerical_metadata",
                                     "bugs", "dist.bugs", "d"),
                       .packages = c("vegan","permute","dplyr","tibble","tidyr")) %dopar% {
                         res <- run_sim_invariable(seed = sim$seed[i], prob = sim$prob[i], d = d,
                                                   bugs_mat = bugs, dist_bugs = dist.bugs)
                         res$index <- sim$index[i]
                         unlist(res)
                       }
    
    parallel::stopCluster(cl)
    cat("Simulation 1 Done. j = ", j, "\n")
    
    
    # --------------------- Invariant: combine & style ---------------------
    res_mat <- do.call(rbind, sim_res)
    df_rs   <- as.data.frame(res_mat, stringsAsFactors = FALSE)
    
    stopifnot("index" %in% colnames(df_rs))
    metrics_names <- setdiff(colnames(df_rs), "index")
    
    df_rs$index <- as.integer(df_rs$index)
    df_rs <- merge(df_rs, sim, by = "index")
    
    df_rs <- reshape2::melt(df_rs, measure.vars = metrics_names)
    df_rs$value <- suppressWarnings(as.numeric(df_rs$value))
    
    # Classify metrics
    effect_types.vc <- metrics_names[grep("_effect", metrics_names)]
    pval_types.vc   <- metrics_names[grep("_p",      metrics_names)]
    
    df_rs$val_type   <- ifelse(df_rs$variable %in% effect_types.vc, "effect", "pval")
    df_rs$perm_type  <- remove_tail(as.character(df_rs$variable), sep = "_", del = 1)
    df_rs$sign_color <- ifelse(df_rs$value < 0.05 & df_rs$val_type == "pval", "red", "blue")
    
    # Method labels in desired order (same as your context)
    df_rs$perm_type_lab <- factor(df_rs$perm_type,
                                  levels = c("unadj","adj_v2","adj_v1","adj_G1","adj_G2","adj_B","hajek"),
                                  labels = c("Unadjusted","~Var+Subj","~Subj+Var","G1","G2","Within block","Hajek"))
    
    # Rename "Within block" -> "Between block" (avoid forcats dependency)
    levels(df_rs$perm_type_lab)[levels(df_rs$perm_type_lab) == "Within block"] <- "Between block"
    
    # (Optional) Limit to the 4 methods you want to show
    df_rs <- df_rs[df_rs$perm_type_lab %in% c("Unadjusted","~Var+Subj","Between block","Hajek"), ]
    
    # Build "Corruption Prob. = x" label from "Flip Prob: x" and reverse row order
    cor_vals <- as.numeric(gsub("Flip Prob: ", "", as.character(df_rs$prob_lab)))
    cor_labs <- paste0("Corruption Prob. = ", cor_vals)
    df_rs$prob_lab <- factor(cor_labs,
                             levels = rev(sort(unique(cor_labs))))
    
    # Housekeeping
    df_rs$line_height <- ifelse(df_rs$val_type == "pval", n_iterations * 0.05, 0)
    
    # Save the styled data frame (invariant)
    save(df_rs, file = paste0(dir_output, "data/permanova_sim_invariant.Rda"))
    
    # Split into pval/effect
    df_pval   <- df_rs[df_rs$val_type == "pval", ]
    df_effect <- df_rs[df_rs$val_type == "effect", ]
    
    # ---------- P-VALUE HISTOGRAM (facet_grid: rows=Corruption Prob., cols=method) ----------
    hist_breaks <- seq(0, 1, by = 0.05)
    
    # Invariant theme (lightgray strips as in your context)
    my_theme_inv <- my_theme + theme(strip.background = element_rect(fill = "lightgray",
                                                                     colour = "black", size = 1))
    
    nr <- nlevels(df_pval$prob_lab)
    nc <- nlevels(df_pval$perm_type_lab)
    
    ggp_pval_inv <- ggplot(df_pval, aes(x = value, fill = sign_color)) +
      geom_histogram(breaks = hist_breaks, color = "black", closed = "left") +
      facet_grid(rows = vars(prob_lab), cols = vars(perm_type_lab), scales = "free_y") +
      geom_hline(aes(yintercept = line_height), colour = "black") +
      my_theme_inv +
      scale_fill_manual(values = c("red" = "#3ded97", "blue" = "#F8766D")) +
      labs(title = "Invariant Metadata", x = "P-value", y = "Count")
    
    # Save a single grid PNG for invariant p-values
    save_plot(paste0(dir_output, "plots/invariant_pval_grid.png"),
              ggp_pval_inv, dpi = 300, base_height = (nr*2), base_width = (nc*2))
    
    # ---------- VOLCANO PLOT (effect vs -log10(p), free scales, per method × corruption prob.) ----------
    # Build join keys to pair effect & pval from same (iteration, prob, seed, method)
    df_effect_volcano <- df_effect %>% mutate(ipsv = paste(iterations, prob, seed, perm_type, sep = "_"))
    df_pval_volcano   <- df_pval   %>% mutate(ipsv = paste(iterations, prob, seed, perm_type, sep = "_"))
    
    df_volcano <- inner_join(
      df_effect_volcano %>% dplyr::select(ipsv, prob, prob_lab, perm_type_lab, effect = value),
      df_pval_volcano   %>% dplyr::select(ipsv, pval = value),
      by = "ipsv"
    ) %>%
      mutate(neglog10p = -log10(pval),
             significant = pval < 0.05)
    
    n_methods_inv <- nlevels(df_volcano$perm_type_lab)
    
    g_volcano_inv <- ggplot(df_volcano, aes(x = effect, y = neglog10p, color = significant)) +
      geom_point(alpha = 0.6, size = 1.5) +
      facet_wrap(~ prob_lab + perm_type_lab, scales = "free", ncol = n_methods_inv) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      my_theme_inv +
      scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      labs(title = "Invariant Metadata",
           x = "Effect size",
           y = expression(-log10))
    
    save_plot(paste0(dir_output, "plots/invariant_volcano_free.png"),
              g_volcano_inv, dpi = 300, base_height = (nr*2), base_width = (n_methods_inv*2))
    
    # ---------- False Positive Rate at Corruption Prob. = 0.5 ----------
    fpr_0.5 <- df_pval %>%
      dplyr::filter(prob == 0.5) %>%   # use numeric prob for robustness
      dplyr::group_by(perm_type_lab) %>%
      dplyr::summarise(
        FPR = mean(value < 0.05),
        n   = dplyr::n(),
        .groups = "drop"
      )
    
    write.csv(fpr_0.5, file = paste0(dir_output, "data/invariant_fpr_prob0.5.csv"), row.names = FALSE)
    print(fpr_0.5)
    
    # ---------- Kolmogorov–Smirnov diagnostics at Corruption Prob. = 0.5 ----------
    # (Robust form without dplyr::pick to avoid version dependency)
    ks_results_inv <- df_pval %>%
      dplyr::filter(prob == 0.5) %>%
      dplyr::group_by(perm_type_lab) %>%
      dplyr::summarise(
        ks_D = {
          vals  <- value
          valsJ <- jitter(vals, amount = 1e-10)
          ks.test(valsJ, "punif", 0, 1)$statistic
        },
        ks_p = {
          vals  <- value
          valsJ <- jitter(vals, amount = 1e-10)
          ks.test(valsJ, "punif", 0, 1)$p.value
        },
        FPR_0.05 = mean(value < 0.05),
        n = dplyr::n(),
        .groups = "drop"
      )
    
    write.csv(ks_results_inv, file = paste0(dir_output, "data/invariant_ks_prob0.5.csv"), row.names = FALSE)
    print(ks_results_inv)
    
    
  } else {
    # ------------------------- VARIANT ---------------------------
    cat("Simulation 2 (Variant) Running...\n")
    
    # Per-subject block metadata once (outside workers)
    df_subject <- data.frame(subject = unique(d$subject),
                             my_block_var = rbinom(length(unique(d$subject)), 1, 0.5))
    rownames(df_subject) <- df_subject$subject
    
    cl <- parallel::makeCluster(nCores, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    
    sim_res <- foreach(i = 1:nrow(sim),
                       .export   = c("run_sim_variable",
                                     "PERMANOVA_repeat_measures",
                                     "hajek_repeat_measures",
                                     "numerical_metadata",
                                     "bugs", "d", "df_subject"),
                       .packages = c("vegan","permute","dplyr","tibble","tidyr")) %dopar% {
                         res <- run_sim_variable(seed = sim$seed[i], prob = sim$prob[i],
                                                 d = d, df_subject = df_subject, bugs_mat = bugs)
                         res$index <- sim$index[i]
                         unlist(res)
                       }
    
    parallel::stopCluster(cl)
    cat("Simulation 2 Done. j = ", j, "\n")
    
    
    # --------------------- Variant: combine & style ---------------------
    res_mat <- do.call(rbind, sim_res)
    df_rs   <- as.data.frame(res_mat, stringsAsFactors = FALSE)
    
    stopifnot("index" %in% colnames(df_rs))
    metrics_names <- setdiff(colnames(df_rs), "index")
    
    df_rs$index <- as.integer(df_rs$index)
    df_rs <- merge(df_rs, sim, by = "index")
    
    df_rs <- reshape2::melt(df_rs, measure.vars = metrics_names)
    df_rs$value <- suppressWarnings(as.numeric(df_rs$value))
    
    # Classify metrics
    effect_types.vc <- metrics_names[grep("_effect", metrics_names)]
    pval_types.vc   <- metrics_names[grep("_p",      metrics_names)]
    
    df_rs$val_type   <- ifelse(df_rs$variable %in% effect_types.vc, "effect", "pval")
    df_rs$perm_type  <- remove_tail(as.character(df_rs$variable), sep = "_", del = 1)
    df_rs$sign_color <- ifelse(df_rs$value < 0.05 & df_rs$val_type == "pval", "red", "blue")
    
    # Method labels in your desired order
    df_rs$perm_type_lab <- factor(df_rs$perm_type,
                                  levels = c("unadj","adj_v2","adj_v1","adj_G1","adj_G2","adj_B","hajek"),
                                  labels = c("Unadjusted","~Var+Subj","~Subj+Var","G1","G2","Within block","Hajek"))
    
    # Build "Effect Size = ..." label from prob and reverse row order
    es_vals <- 50 - (50*2) * df_rs$prob
    es_labs <- paste0("Effect Size = ", es_vals)
    
    df_rs$prob_lab <- factor(es_labs,
                             levels = unique(paste0("Effect Size = ", sort(unique(es_vals), decreasing = FALSE))))
    
    # Housekeeping
    df_rs$line_height <- ifelse(df_rs$val_type == "pval", n_iterations * 0.05, 0)
    
    # Save the styled data frame
    save(df_rs, file = paste0(dir_output, "data/permanova_sim_variant.Rda"))
    
    # Split into pval/effect
    df_pval   <- df_rs[df_rs$val_type == "pval", ]
    df_effect <- df_rs[df_rs$val_type == "effect", ]
    
    # ---------- P-VALUE HISTOGRAM (facet_grid: rows=effect size, cols=method) ----------
    hist_breaks <- seq(0, 1, by = 0.05)
    
    # Variant theme (lightgray strips like your example)
    my_theme_var <- my_theme + theme(strip.background = element_rect(fill = "lightgray",
                                                                     colour = "black", size = 1))
    
    nr <- nlevels(df_pval$prob_lab)
    nc <- nlevels(df_pval$perm_type_lab)
    
    #df_pval$prob_lab <- factor(df_pval$prob_lab, levels = rev(levels(df_pval$prob_lab)))
    #df_effect$prob_lab <- factor(df_pval$prob_lab, levels = rev(levels(df_pval$prob_lab)))
    
    ggp_pval <- ggplot(df_pval, aes(x = value, fill = sign_color)) +
      geom_histogram(breaks = hist_breaks, color = "black", closed = "left") +
      facet_grid(rows = vars(prob_lab), cols = vars(perm_type_lab), scales = "free_y") +
      geom_hline(aes(yintercept = line_height), colour = "black") +
      my_theme_var +
      scale_fill_manual(values = c("red" = "#3ded97", "blue" = "#F8766D")) +
      labs(title = "Variant Metadata", x = "P-value", y = "Count")
    
    # Save a single grid PNG
    save_plot(paste0(dir_output, "plots/variant_pval_grid.png"),
              ggp_pval, dpi = 300, base_height = (nr*2), base_width = (nc*2))
    
    # ---------- EFFECT SIZE HISTOGRAM (facet_wrap_paginate: effect size × method) ----------
    nr_eff <- nr
    nc_eff <- nc
    
    ggp_effect <- ggplot(df_effect, aes(x = value)) +
      geom_histogram(binwidth = 0.05, color = "black", fill = "black") +
      ggforce::facet_wrap_paginate(~ prob_lab + perm_type_lab, nrow = nr_eff, ncol = nc_eff) +
      geom_hline(aes(yintercept = line_height), colour = "black") +
      my_theme_var +
      labs(title = "Variant Metadata", x = "Effect size", y = "Count")
    
    num_pages_effect <- ggforce::n_pages(ggp_effect)
    for (i in 1:num_pages_effect) {
      ggp_i <- ggp_effect + ggforce::facet_wrap_paginate(~ prob_lab + perm_type_lab,
                                                         nrow = nr_eff, ncol = nc_eff,
                                                         scales = "fixed", page = i)
      save_plot(paste0(dir_output, "plots/variant_effect_page_", i, ".png"),
                ggp_i, dpi = 300, base_height = (nr_eff*2), base_width = (nc_eff*2))
    }
    
    # ---------- VOLCANO PLOT (effect vs -log10(p), free scales, per method × effect size) ----------
    # Build join keys to pair effect & pval from same (iteration, prob, seed, method)
    df_effect <- df_effect %>%
      mutate(ipsv = paste(iterations, prob, seed, perm_type, sep = "_"))
    df_pval <- df_pval %>%
      mutate(ipsv = paste(iterations, prob, seed, perm_type, sep = "_"))
    
    df_volcano <- inner_join(
      df_effect %>% dplyr::select(ipsv, prob, prob_lab, perm_type_lab, effect = value),
      df_pval   %>% dplyr::select(ipsv, pval = value),
      by = "ipsv"
    ) %>%
      mutate(neglog10p = -log10(pval),
             significant = pval < 0.05)
    #df_volcano$prob_lab <- factor(df_pval$prob_lab, levels = rev(levels(df_pval$prob_lab)))
    
    
    # Use free scales per facet; one row per (Effect Size × Method)
    n_methods <- nlevels(df_volcano$perm_type_lab)
    g_volcano <- ggplot(df_volcano, aes(x = effect, y = neglog10p, color = significant)) +
      geom_point(alpha = 0.6, size = 1.5) +
      facet_wrap(~ prob_lab + perm_type_lab, scales = "free", ncol = n_methods) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      my_theme_var +
      scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
      labs(title = "Variant Metadata",
           x = "Effect size",
           y = expression(-log10))
    
    save_plot(paste0(dir_output, "plots/variant_volcano_free.png"),
              g_volcano, dpi = 300, base_height = (nr*2), base_width = (n_methods*2))
    
    # ---------- Kolmogorov–Smirnov diagnostics (at prob == 0.5) ----------
    # Compute KS against U(0,1) for p-values at highest effect size (prob == 0.5)
    ks_results <- df_pval %>%
      dplyr::filter(prob == 0.5) %>%                  # robust to label changes
      dplyr::group_by(perm_type_lab) %>%
      dplyr::summarise(
        ks_D = {
          vals  <- dplyr::pick(value) |> dplyr::pull()
          valsJ <- jitter(vals, amount = 1e-10)
          ks.test(valsJ, "punif", 0, 1)$statistic
        },
        ks_p = {
          vals  <- dplyr::pick(value) |> dplyr::pull()
          valsJ <- jitter(vals, amount = 1e-10)
          ks.test(valsJ, "punif", 0, 1)$p.value
        },
        FPR_0.05 = mean(value < 0.05),
        n = dplyr::n(),
        .groups = "drop"
      )
    
    # Save KS table
    write.csv(ks_results, file = paste0(dir_output, "data/variant_ks_prob0.5.csv"), row.names = FALSE)
    print(ks_results)
    
  }#variant
}#for j in 1

cat("All done.\n")
