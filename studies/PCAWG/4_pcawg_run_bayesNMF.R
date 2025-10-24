# identify task id from command line arguments
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
print(task_id)

# identify counts matrix from job assignments & task id
job_assignments <- read.csv("pcawg_assignments.csv")
matrix <- job_assignments$matrix[task_id]
name <- substr(matrix, 1, nchar(matrix) - 4) # drops .csv

# load bayesNMF R package
library(bayesNMF)

# load counts data (already nonhyper only)
data <- read.csv(file.path("../../processed", "PCAWG", "matrices_nonhyper", matrix), row.names = 1)
data <- as.matrix(data)

# run bayesNMF
res <- bayesNMF(
  data,
  rank = 1:20,
  likelihood = "poisson",
  prior = "truncnormal",
  MH = TRUE,
  rank_method = "SBFI",
  output_dir = file.path("../../output", "PCAWG", "bayesNMF", name),
  overwrite = TRUE,
  periodic_save = FALSE,
  save_all_samples = FALSE
)
