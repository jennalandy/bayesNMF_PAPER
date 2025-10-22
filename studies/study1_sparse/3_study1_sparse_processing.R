library(tidyverse)
library(glue)
library(bayesNMF)

data_dir <- "../../data/study1_sparse"
output_dir <- "../../output/study1_sparse"
processed_dir <- "../../processed/study1_sparse"

source("../study1_process.R")

# just checking what is done so far
counts <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(counts) <- c("model", "N", "G", "rep")
for(model_dir in list.dirs(output_dir, full.names = TRUE, recursive = FALSE)) {
  model_name <- basename(model_dir)
  results_dirs <- list.files(model_dir, full.names = TRUE, recursive = FALSE)

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
  arrange(model) %>%
  filter(n < 10) 

########################################################
# Individual model metrics
########################################################
print("Processing individual model metrics")
model_dirs <- list.dirs(output_dir, full.names = TRUE, recursive = FALSE)
model_dirs  <- model_dirs[!grepl("repeated", model_dirs)]
results <- lapply(model_dirs, function(model_dir) {
  print(model_dir)
  model_name <- basename(model_dir)
  model_name_parts <- str_split(model_name, "_")[[1]]
  likelihood <- model_name_parts[1]
  prior <- model_name_parts[2]
  MH <- grepl("MH", model_name)

  results_dirs <- list.dirs(model_dir, full.names = TRUE, recursive = FALSE)
  lapply(results_dirs, function(results_dir) {
    sampler_path <- file.path(results_dir, "sampler.rds")
    if (!file.exists(sampler_path)) {
      return(NULL)
    }
    print(results_dir)
    sampler <- readRDS(sampler_path)
    final_metrics <- sampler$state$MAP_metrics %>% filter(iter == sampler$state$iter)

    results_name <- basename(results_dir)
    results_name_parts <- str_split(results_name, "_")[[1]]
    N <- results_name_parts[1] %>% str_remove("N") %>% as.numeric()
    G <- results_name_parts[2] %>% str_remove("G") %>% as.numeric()
    rep <- results_name_parts[3] %>% str_remove("rep") %>% as.numeric()
    
    data <- readRDS(file.path(data_dir, glue("N{N}_G{G}_rep{rep}.rds")))
    assignments <- hungarian_assignment(sampler$MAP$P, reference_P = data$P)

    metrics <- list(
      N = N,
      G = G,
      rep = rep,
      likelihood = likelihood,
      prior = prior,
      MH = MH,
      iters = sampler$state$iter,
      min_sim = min(assignments$cos_sim),
      RMSE = final_metrics$RMSE,
      KL = final_metrics$KL,
      total_minutes = as.numeric(sampler$time$total),
      warmup_minutes = ifelse(MH, as.numeric(sampler$time$warmup), 0),
      minutes_per_iter = as.numeric(sampler$time$per_iter),
      avg_mut_per_sample = mean(colMeans(data$M))
    )
    return(metrics)
  }) %>%
    do.call(rbind, .)
})

results_df <- results %>%
  do.call(rbind, .) %>%
  data.frame()

# columns are currently lists, so convert to characters and then to numeric
results_df_clean <- dplyr::mutate(results_df, across(everything(), as.character)) %>%
  dplyr::mutate(across(c(N, G, rep, min_sim, RMSE, KL, total_minutes, warmup_minutes, minutes_per_iter), as.numeric))
write_csv(results_df_clean, file.path(processed_dir, "metrics.csv"))

results_df_clean %>%
  group_by(N, G, likelihood, prior, MH) %>%
  tally() %>%
  arrange(n)


########################################################
# Comparing MAPs
########################################################
print("Processing MAP comparisons")

Standard_vs_MH_MAPs_Exponential <- process_compare_MAPs("exponential", c("Poisson", "Poisson_MH"))  
Standard_vs_MH_MAPs <- Standard_vs_MH_MAPs_Exponential %>%
  mutate(prior = 'Exponential')
write_csv(Standard_vs_MH_MAPs, file.path(processed_dir, "Standard_vs_MH_MAPs.csv"))

Standard_vs_Standard_repeated_MAPs_Exponential <- process_compare_MAPs("exponential", c("Poisson", "Poisson_repeated"))
Standard_vs_Standard_repeated_MAPs <- Standard_vs_Standard_repeated_MAPs_Exponential %>%
  mutate(prior = 'Exponential')
write_csv(Standard_vs_Standard_repeated_MAPs, file.path(processed_dir, "Standard_vs_Standard_repeated_MAPs.csv"))


########################################################
# Comparing credible intervals
########################################################

print("Processing credible interval comparisons")

Standard_vs_MH_CIs_Exponential <- process_compare_credible_intervals("exponential", c("Poisson", "Poisson_MH")) %>%
  mutate(prior = 'Exponential')

Standard_vs_MH_CIs_Exponential_repeated <- process_compare_credible_intervals("exponential", c("Poisson", "Poisson_repeated")) %>%
  mutate(prior = 'Exponential')

Standard_vs_MH_CIs <- do.call(rbind, list(
  Standard_vs_MH_CIs_Exponential,
  Standard_vs_MH_CIs_Exponential_repeated
))
write_csv(Standard_vs_MH_CIs, file.path(processed_dir, "Standard_vs_MH_CIs.csv"))