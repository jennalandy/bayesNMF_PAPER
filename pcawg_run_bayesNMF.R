# identify task id from command line arguments
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
print(task_id)

# identify counts matrix from job assignments & task id
job_assignments <- read.csv("pcawg_assignments.csv")
matrix <- job_assignments$matrix[task_id]
name <- substr(matrix, 1, nchar(matrix) - 4) # drops .csv

# load bayesNMF R package
library(bayesNMF, lib = "/homes8/jlandy/R/x86_64-pc-linux-gnu-library/4.3")

# load counts data (already nonhyper only)
data <- read.csv(file.path("pcawg_counts_matrices_nonhyper", matrix), row.names = 1)
data <- as.matrix(data)

# run bayesNMF
res <- bayesNMF(
  data,
  rank = 1:20,
  likelihood = "poisson",
  prior = "truncnormal",
  fast = TRUE,
  learn_rank_method = "SBFI",
  file = file.path("output", "pcawg_bayesNMF", name),
  overwrite = TRUE
)
