library(bayesNMF)
library(glue)
job_assignments <- read.csv("study2_assignments.csv")
data_dir <- "../../data/study2"

########################################################
############# Model Specifications #####################
########################################################
likelihood <- "poisson"
prior <- "exponential"
MH <- TRUE
rank_method <- "BFI"

L <- toupper(substr(likelihood, 1, 1))
P <- toupper(substr(prior, 1, 1))
output_dir <- glue("../../output/study2/{L}_{P}{ifelse(MH, '_MH', '')}_{rank_method}")
print(output_dir)
########################################################

# identify task id from command line arguments
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
print(task_id)
assignment <- job_assignments[job_assignments$i == task_id, ]

for (rep in 1:10) {
  # Load data
  name <- glue("N{assignment$N}_G{assignment$G}_rep{rep}")
  data <- readRDS(file.path(data_dir, glue("{name}.rds")))

  stopifnot(assignment$N == ncol(data$P))
  if (file.exists(file.path(output_dir, name, "sampler.rds"))) {
    sampler <- readRDS(file.path(output_dir, name, "sampler.rds"))
    if (sampler$state$converged) {
      print(glue("Skipping {name} because it already exists and has converged"))
      next
    }
  }

  print(glue("Running {name}"))

  # Run bayesNMF
  res <- bayesNMF(
    data = as.matrix(data$M),
    rank = 1:20,
    rank_method = rank_method,
    likelihood = likelihood,
    prior = prior,
    MH = MH,
    output_dir = file.path(output_dir, name),
    overwrite = TRUE,
    periodic_save = FALSE,
    save_all_samples = FALSE
  )

  # Clean up
  rm(list  = c('res','data'))
}