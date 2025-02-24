# identify task id from command line arguments
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
print(task_id)

# identify counts matrix from job assignments & task id
job_assignments <- read.csv("study1_assignments.csv")
assignment <- job_assignments[job_assignments$i == task_id,]


library(bayesNMF, lib = "/homes8/jlandy/R/x86_64-pc-linux-gnu-library/4.3")
library(glue)

data_dir <- "data/study1"
output_dir <- "output/study1"

for (rep in 1:10) {
  name <- glue("N{assignment$N}_G{assignment$G}_rep{rep}")
  filename <- glue("{name}.rds")
  data <- readRDS(file.path(data_dir, file_name))

  stopifnot(assignment$N == ncol(data$P))

  output_file_name <- glue("{name}_NT")

  # run bayesNMF
  res <- bayesNMF(
    as.matrix(data$M),
    rank = assignment$N,
    likelihood = "normal",
    prior = "truncnormal",
    file = file.path(output_dir, output_file_name),
    overwrite = TRUE
  )

  rm(list  = c('res','data'))
}