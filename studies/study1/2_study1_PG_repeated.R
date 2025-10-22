library(bayesNMF)
library(glue)
job_assignments <- read.csv("study1_assignments.csv")
data_dir <- "../../data/study1"

########################################################
############# Model Specifications #####################
########################################################
likelihood <- "poisson"
prior <- "gamma"
MH <- FALSE

L <- toupper(substr(likelihood, 1, 1))
P <- toupper(substr(prior, 1, 1))
output_dir <- glue("../../output/study1/{L}_{P}{ifelse(MH, '_MH', '')}_repeated")
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
    rank = assignment$N,
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