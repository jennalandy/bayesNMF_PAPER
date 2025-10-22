library(tidyverse)
library(ggplot2)
library(bayesNMF)
library(ggh4x)

data_dir <- "../../data/study2"
output_dir <- "../../output/study2"
processed_dir <- "../../processed/study2"
figures_dir <- "../../figures/study2"

# just checking what is done so far
counts <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(counts) <- c("model", "N", "G", "rep")
model_dirs <- rev(list.dirs(output_dir, full.names = TRUE, recursive = FALSE))
model_dirs <- model_dirs[grepl("BIC", model_dirs)]
for(model_dir in model_dirs) {
  model_name <- basename(model_dir)
  results_dirs <- list.files(model_dir, full.names = TRUE, recursive = FALSE)
  results_dirs <- results_dirs[!grepl("withsamples", results_dirs)]

  for(dir in results_dirs) {
    name <- basename(dir)
    print(paste(model_name, name))

    name_parts <- str_split(name, "_")[[1]]
    N <- as.numeric(str_replace(name_parts[1], "N", ""))
    G <- as.numeric(str_replace(name_parts[2], "G", ""))
    rep <- as.numeric(str_replace(name_parts[3], "rep", ""))

    counts <- rbind(counts, data.frame(
      model = model_name,
      N = N,
      G = G,
      rep = rep
    ))
  }
}
counts %>%
  group_by(model, N, G) %>%
  tally() %>%
  arrange(n)

counts %>%
  group_by(model, N, G) %>%
  tally() %>%
  group_by(model, G) %>%
  summarize(n = sum(n)) %>%
  arrange(n)

results <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(results) <- c("model", "N", "G", "rep", "learned_rank", "precision", "sensitivity")
model_dirs <- rev(list.dirs(output_dir, full.names = TRUE, recursive = FALSE))
model_dirs <- model_dirs[grepl("BIC", model_dirs)]
for(model_dir in model_dirs) {
  model_name <- basename(model_dir)
  results_dirs <- list.files(model_dir, full.names = TRUE, recursive = FALSE)
  results_dirs <- results_dirs[!grepl("withsamples", results_dirs)]

  for(dir in results_dirs) {
    if (file.exists(file.path(dir, "sampler.rds"))) {
      sampler <- NULL
      tryCatch({
          sampler <- readRDS(file.path(dir, "sampler.rds"))
      }, error = function(e) {})
      if (is.null(sampler)) {
          print(paste(model_name, dir, "sampler.rds couldn't be opened"))
          next
      }
      if (!sampler$state$converged) {
          print(paste(model_name, dir, "sampler.rds not converged"))
          next
      }
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

    # get learned rank
    learned_rank <- sum(sampler$MAP$A)

    # only allow 1-1 assignment
    sim_mat <- hungarian_assignment(sampler$MAP$P, data$P, return_mat = TRUE)
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

    # total time as sum across all ranks considered
    rank_dirs <- list.files(dir, full.names = TRUE, recursive = FALSE)
    rank_dirs <- rank_dirs[grepl("rank", rank_dirs)]
    total_minutes <- sum(sapply(rank_dirs, function(x) {
      as.numeric(readRDS(file.path(x, "sampler.rds"))$time$total)
    }))

    results <- rbind(results, data.frame(
      model = model_name,
      N = N,
      G = G,
      rep = rep,
      learned_rank = learned_rank,
      precision = precision,
      sensitivity = sensitivity,
      total_minutes = total_minutes
    ))
    write.csv(results, file.path(processed_dir, "metrics_heuristic.csv"), row.names = FALSE)
  }
}

results2 <- results %>%
  mutate(
    likelihood = case_when(
      grepl("P_", model) ~ "Poisson",
      grepl("N_", model) ~ "Normal",
      TRUE ~ model
    ),
    prior = case_when(
      grepl("T_", model) ~ "Truncated Normal",
      grepl("E_", model) ~ "Exponential",
      grepl("G_", model) ~ "Gamma",
      TRUE ~ model
    ),
    MH = grepl("MH", model)
  )
  
write.csv(results2, file.path(processed_dir, "metrics_heuristic.csv"), row.names = FALSE)
