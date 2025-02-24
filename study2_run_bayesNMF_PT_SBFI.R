# identify task id from command line arguments
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
print(task_id)

# identify counts matrix from job assignments & task id
job_assignments <- read.csv("study2_assignments.csv")
assignment <- job_assignments[job_assignments$i == task_id,]


library(bayesNMF, lib = "/homes8/jlandy/R/x86_64-pc-linux-gnu-library/4.3")
library(glue)

data_dir <- "data/study2"
output_dir <- "output/study2"

for (rep in 1:10) {
  file_name <- glue("N{assignment$N}_G{assignment$G}_rep{rep}.rds")
  data <- readRDS(file.path(data_dir, file_name))

  stopifnot(assignment$N == ncol(data$P))

  output_file_name <- glue("N{assignment$N}_G{assignment$G}_rep{rep}_PT_SBFI")

  if (file.exists(file.path(output_dir, paste0(output_file_name, '.rds')))) {
    res <- readRDS(file.path(output_dir, paste0(output_file_name, '.rds')))
    if (!is.null(res$converged_at)) {
      next()
    }
  }

  # run bayesNMF
  res <- bayesNMF(
    as.matrix(data$M),
    rank = 1:20,
    learn_rank_method = 'SBFI',
    likelihood = "poisson",
    prior = "truncnormal",
    fast = TRUE,
    file = file.path(output_dir, output_file_name),
    overwrite = TRUE
  )
}