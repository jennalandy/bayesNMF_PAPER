library(bayesNMF)
library(glue)
job_assignments <- read.csv("study1_sparse_assignments.csv")
data_dir <- "../../data/study1_sparse"

########################################################
############# Model Specifications #####################
########################################################
likelihood <- "poisson"
prior <- "exponential"
MH <- TRUE

L <- toupper(substr(likelihood, 1, 1))
P <- toupper(substr(prior, 1, 1))
output_dir <- glue("../../output/study1_sparse/{L}_{P}{ifelse(MH, '_MH', '')}")
print(output_dir)
########################################################

# identify task id from command line arguments
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
post_fix <- ""
if (length(args) > 1) {
  post_fix <- paste0("_", args[2])
}
print(task_id)
assignment <- job_assignments[job_assignments$i == task_id, ]

for (rep in 1:10) {
  # Load data
  name <- glue("N{assignment$N}_G{assignment$G}_rep{rep}")
  data <- readRDS(file.path(data_dir, glue("{name}.rds")))

  name <- paste0(name, post_fix)

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
    periodic_save = FALSE
  )

  # Clean up
  rm(list  = c('res','data'))
}