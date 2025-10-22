library(tidyverse)
library(ggplot2)
library(bayesNMF)
library(ggh4x)

data_dir <- "../../data/study2"
output_dir <- "../../output/study2"
processed_dir <- "../../processed/study2"
figures_dir <- "../../figures/study2"

# SignatureAnalyzer Results
results <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(results) <- c("model", "N", "G", "rep", "learned_rank", "precision", "sensitivity")
dirs <- rev(list.dirs(output_dir, full.names = TRUE, recursive = FALSE))
dirs <- dirs[grepl("signatureanalyzer", dirs)]
reference <- get_cosmic()

for(model_dir in dirs) {
  model_name <- basename(model_dir)
  results_dirs <- list.files(model_dir, full.names = TRUE, recursive = FALSE)

  for(dir in results_dirs) {
    if (file.exists(file.path(dir, "W.csv"))) {
      Phat <- read.csv(file.path(dir, "W.csv")) %>%
        column_to_rownames("X") %>%
        select(starts_with('S'))
      learned_rank = ncol(Phat)
      print(paste(name, learned_rank))
    } else {
      next
    }

    # get true rank N and sample size G and rep from dir names
    name <- basename(dir)
    print(paste(model_name, name))

    name_parts <- str_split(name, "_")[[1]]
    N <- as.numeric(str_replace(name_parts[1], "N", ""))
    G <- as.numeric(str_replace(name_parts[2], "G", ""))
    rep <- as.numeric(str_replace(name_parts[3], "rep", ""))

    data <- readRDS(file.path(data_dir, paste0(name, ".rds")))

    # signature analyzer ignored all 0 rows -- add these rows as all 0s to P
    zero_rownames <- rownames(reference)[rowSums(data$M) == 0]
    setdiff(rownames(reference), rownames(Phat))
    if (length(zero_rownames) > 0) {
      fill_zeros <- matrix(
        rep(0, ncol(Phat) * length(zero_rownames)),
        nrow = length(zero_rownames)
      )
      rownames(fill_zeros) <- zero_rownames
      colnames(fill_zeros) <- colnames(Phat)
      Phat <- rbind(
        Phat,
        fill_zeros
      )
      Phat <- Phat[rownames(reference),]
      if (learned_rank == 1) {  Phat <- matrix(Phat, ncol = 1) }
    }

    # only allow 1-1 assignment
    sim_mat <- hungarian_assignment(Phat, data$P, return_mat = TRUE)
    considered <- min(learned_rank, N)
    sim_mat <- sim_mat[1:considered, 1:considered]
    if (considered == 1) {
      sim_mat <- matrix(sim_mat, nrow = 1, ncol = 1)
    }
    N_match <- sum(rowSums(sim_mat > 0.9))

    # compute precision as propotion of estimated factors with a true match > 0.9 cosine sim 
    precision <- N_match / learned_rank

    # compute sensitivity as propotion of true factors with a learned match > 0.9 cosine sim
    sensitivity <- N_match / N

    results <- rbind(results, data.frame(
      model = model_name,
      N = N,
      G = G,
      rep = rep,
      learned_rank = learned_rank,
      precision = precision,
      sensitivity = sensitivity
    ))
    write.csv(results, file.path(processed_dir, "metrics_signatureanalyzer.csv"), row.names = FALSE)
  }
}
results2 <- results %>%
  mutate(
    likelihood = "Poisson",
    prior = case_when(
      grepl("L2", model) ~ "Truncated Normal",
      grepl("L1", model) ~ "Exponential",
      TRUE ~ model
    ),
    MH = FALSE
  )
results3 <- read.csv(file.path(output_dir, "signatureanalyzer_L2_runtime.csv")) %>% 
  mutate(
    model = "signatureanalyzer_L2"
  ) %>% 
  rbind(
    read.csv(file.path(output_dir, "signatureanalyzer_L1_runtime.csv")) %>%
      mutate(
        model = "signatureanalyzer_L1"
      )
  ) %>%
  merge(results2, by = c("model", "N", "rep"), all = TRUE)

write.csv(results3, file.path(processed_dir, "metrics_signatureanalyzer.csv"), row.names = FALSE)
