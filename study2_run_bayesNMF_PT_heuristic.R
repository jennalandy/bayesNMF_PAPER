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

  output_file_name <- glue("N{assignment$N}_G{assignment$G}_rep{rep}_PT_heuristic")

  if (file.exists(file.path(output_dir, paste0(output_file_name, '.rds')))) {
    next()
  }

  # run bayesNMF
  res <- bayesNMF(
    as.matrix(data$M),
    rank = 1:20,
    learn_rank_method = 'heuristic',
    likelihood = "poisson",
    prior = "truncnormal",
    fast = TRUE,
    file = file.path(output_dir, output_file_name),
    overwrite = TRUE,
    store_logs = FALSE
  )

  files <- list.files(output_dir)
  intermediate_files <- files[
    grepl(output_file_name, files) &
    grepl("rank", files)
  ]
  for (f in intermediate_files) {
    file.remove(f)
  }

  rm(list = c('res','data'))
}