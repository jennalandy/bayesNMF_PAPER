library(tidyverse)

calc_normal_BIC <- function(M, Mhat, sigmasq, n) {
  K <- nrow(Mhat)
  G <- ncol(Mhat)
  loglik <- (- K * sum(log(2 * pi * sigmasq)) / 2 -
      sum(colSums((M - Mhat)**2) / (2*sigmasq)))

  n_params = n * (G + K)

  return(n_params * log(G) - 2 * loglik)
}

calc_poisson_BIC <- function(M, Mhat, n) {
  logfac = vector(length = max(M))
  logfac[1] = 0
  for (i in 2:length(logfac)) {
    logfac[i] = log(i) + logfac[i-1]
  }

  K <- nrow(Mhat)
  G <- ncol(Mhat)
  Mhat[Mhat <= 0] <- 1
  M[M <= 0] <- 1
  loglik <- sum(sapply(1:K, function(k) {
      sum(sapply(1:G, function(g) {
          (-Mhat[k,g] +
          M[k,g] * log(Mhat[k,g]) -
          logfac[M[k,g]])
      }))
  }))

  n_params = n * (G + K)

  return(n_params * log(G) - 2 * loglik)
}

# Function to perform hungarian algorithm using RcppHungarian
hungarian_algorithm <- function(matrix, which = "max") {
  if (which == "max") {
    matrix <- -1 * matrix
  }
  hungarian_res <- RcppHungarian::HungarianSolver(matrix)
  hungarian_alignment <- data.frame(hungarian_res$pairs) %>%
    filter(X1 != 0 & X2 != 0)

  aligned_matrix <- matrix[
    hungarian_alignment$X1,
    hungarian_alignment$X2
  ]

  if (which == "max") {
    aligned_matrix <- -1 * aligned_matrix
  }

  if (nrow(matrix) == 1) {
    aligned_matrix = matrix(aligned_matrix)
  }
  return(list(
    "mat" = aligned_matrix,
    "pairs" = hungarian_alignment
  ))
}

KL_div <- function(M, Mhat, fill = 0.001) {
  Mhat[Mhat == 0] = fill
  M[M == 0] = fill
  sum(M * log(M / Mhat) - M + Mhat)
}

get_fixedN_metrics <- function(
  results_dir, data_dir
) {
  metrics <- matrix(nrow = 0, ncol = 20)

  files <- list.files(results_dir)
  files <- files[grepl(".rds", files)]
  for (file in files) {
    print(file)

    file_parts <- str_split(file, "_")[[1]]
    if (grepl("fast", file)) {
      model <- paste(file_parts[1:2], collapse = "_")
      N <- file_parts[3]
      G <- file_parts[4]
      rep <- file_parts[5]
    } else {
      model <- file_parts[1]
      N <- file_parts[2]
      G <- file_parts[3]
      rep <- file_parts[4]
    }
    

    tryCatch({

      res <- readRDS(file.path(results_dir, file))
      if (is.null(res$converged_at)) {
        print(paste(file, "not done"))
        next()
      }
      data <- readRDS(file.path(
        data_dir, paste(N, G, rep, "cosmic.rds", sep = "_")
      ))

      # rescale results
      Phat <- res$MAP$P
      Ehat <- sweep(
        res$MAP$E,
        1,
        colSums(Phat),
        '*'
      )
      Phat <- sweep(
        Phat,
        2,
        colSums(Phat),
        '/'
      )
      Mhat <- Phat %*% Ehat

      # metrics
      P_similarity_matrix <- pairwise_sim(data$P, Phat, which = "cols")
      P_hung <- hungarian_algorithm(
        P_similarity_matrix, which = "max"
      )
      P_aligned_similarity_matrix <- P_hung$mat
      Phat <- Phat[,P_hung$pairs$X2]

      E_similarity_matrix <- pairwise_sim(data$E, Ehat, which = "rows")
      E_hung <- hungarian_algorithm(
        E_similarity_matrix, which = "max"
      )
      E_aligned_similarity_matrix <- E_hung$mat
      Ehat <- Ehat[E_hung$pairs$X2,]

      Mhat_no0 = Mhat
      Mhat_no0[Mhat_no0==0] = 0.01
      KL = sum(data$M * log(data$M / Mhat_no0) - data$M + Mhat_no0)

      this_metrics <- list(
        model = model,
        file = file,
        N = N,
        G = G,
        rep = rep,
        RMSE = sqrt(mean((data$M - Mhat)**2)),
        RMSE_P = sqrt(mean((data$P - Phat)**2)),
        RMSE_E = sqrt(mean((data$E - Ehat)**2)),
        KL = KL,
        P_mean_sim = mean(diag(P_aligned_similarity_matrix)),
        P_median_sim = median(diag(P_aligned_similarity_matrix)),
        P_max_sim = max(diag(P_aligned_similarity_matrix)),
        P_min_sim = min(diag(P_aligned_similarity_matrix)),
        E_mean_sim = mean(diag(E_aligned_similarity_matrix)),
        E_median_sim = median(diag(E_aligned_similarity_matrix)),
        E_max_sim = max(diag(E_aligned_similarity_matrix)),
        E_min_sim = min(diag(E_aligned_similarity_matrix)),
        iters = res$converged_at,
        total_seconds = res$time$total_secs,
        secs_per_iter = res$time$avg_secs_per_iter
      )

      metrics <- rbind(metrics, unlist(this_metrics))
      colnames(metrics) <- names(this_metrics)

    }, error = function(e) {
      print(paste(file, "can't open"))
    })
  }
  return(metrics)
}


get_learnN_metrics <- function(
  results_dir,
  model,
  max_cos = 1,
  data_dir = "../simulation_poisson/data/learnN3/",
  recovery = FALSE,
  G = 64
) {
  metrics <- matrix(nrow = 0, ncol = 27)

  files <- list.files(file.path(results_dir, model, paste0("max", max_cos)))
  files <- files[grepl(".rds", files)]
  files <- files[!grepl("_rank", files)]
  files <- files[grepl(paste0("_", G, "_"), files)]
  if (recovery) {
    files <- files[grepl("recovery", files)]
  } else {
    files <- files[!grepl("recovery", files)]
  }
  for (file in files) {
    print(file)
    file_parts <- str_split(str_replace(file, ".rds", ""), "_")[[1]]
    N <- file_parts[2]
    G <- file_parts[3]
    rep <- file_parts[5]

    res <- readRDS(file.path(results_dir, model, paste0("max", max_cos), file))
    if (is.null(res$converged_at)) {
      next()
    }
    data <- readRDS(file.path(
      data_dir,
      paste0("max", max_cos),
      paste(N, G, rep, "cosmic.rds", sep = "_")
    ))


    # rescale results
    Phat <- res$MAP$P
    if (sum(res$MAP$A) == 1) {
      Ehat <- matrix(res$MAP$E * sum(Phat), nrow = 1)
      Phat <- matrix(Phat / sum(Phat), ncol = 1)
    } else {
      Ehat <- sweep(
        res$MAP$E,
        1,
        colSums(Phat),
        '*'
      )
      Phat <- sweep(
        Phat,
        2,
        colSums(Phat),
        '/'
      )
    }
    Mhat <- Phat %*% Ehat

     if (sum(res$MAP$A) > 0) {
      P_similarity_matrix <- pairwise_sim(data$P, Phat, which = "cols")
      P_hung <- hungarian_algorithm(
        P_similarity_matrix, which = "max"
      )
      P_aligned_similarity_matrix <- P_hung$mat

      E_similarity_matrix <- pairwise_sim(data$E, Ehat, which = "rows")
      E_hung <- hungarian_algorithm(
        E_similarity_matrix, which = "max"
      )
      E_aligned_similarity_matrix <- E_hung$mat

      P_mean_sim = mean(diag(P_aligned_similarity_matrix))
      P_median_sim = median(diag(P_aligned_similarity_matrix))
      P_max_sim = max(diag(P_aligned_similarity_matrix))
      P_min_sim = min(diag(P_aligned_similarity_matrix))
      E_mean_sim = mean(diag(E_aligned_similarity_matrix))
      E_median_sim = median(diag(E_aligned_similarity_matrix))
      E_max_sim = max(diag(E_aligned_similarity_matrix))
      E_min_sim = min(diag(E_aligned_similarity_matrix))
      precision = calc_precision(Phat, data$P)
      sensitivity = calc_sensitivity(Phat, data$P)
    } else {
      P_mean_sim = NA
      P_median_sim = NA
      P_max_sim = NA
      P_min_sim = NA
      E_mean_sim = NA
      E_median_sim = NA
      E_max_sim = NA
      E_min_sim = NA
      precision = NA
      sensitivity = NA
    }

    # metrics
    Mhat[Mhat == 0] <- 1

    this_metrics <- list(
      file = file,
      N = N,
      G = G,
      rep = rep,
      RMSE = sqrt(mean((data$M - Mhat)**2)),
      KL = sum(data$M * log(data$M / Mhat) - data$M + Mhat),
      P_mean_sim = P_mean_sim,
      P_median_sim = P_median_sim,
      P_max_sim = P_max_sim,
      P_min_sim = P_min_sim,
      E_mean_sim = E_mean_sim,
      E_median_sim = E_median_sim,
      E_max_sim = E_max_sim,
      E_min_sim = E_min_sim,
      iters = res$converged_at,
      total_seconds = res$time$total_secs,
      secs_per_iter = res$time$avg_secs_per_iter,
      estimated_N = sum(res$MAP$A),
      precision = precision,
      sensitivity = sensitivity,
      recovery = grepl("recovery", file),
      mean_sigmasq = mean(res$MAP$sigmasq),
      mean_M = mean(data$M),
      normal_BIC = calc_normal_BIC(data$M, Mhat, res$MAP$sigmasq, sum(res$MAP$A)),
      poisson_BIC = calc_poisson_BIC(data$M, Mhat, sum(res$MAP$A)),
      normal_BIC_TRUE = calc_normal_BIC(data$M, data$P%*%data$E, colMeans(data$P%*%data$E), as.numeric(N)),
      poisson_BIC_TRUE = calc_poisson_BIC(data$M, data$P%*%data$E, as.numeric(N))
    )

    metrics <- rbind(metrics, unlist(this_metrics))
    colnames(metrics) <- names(this_metrics)
  }
  return(metrics)
}

get_heuristic_metrics <- function(
  results_dir = "results/learnN",
  max_cos = 0.55,
  data_dir = "../simulation_poisson/data/learnN3/"
) {
  metrics <- matrix(nrow = 0, ncol = 27)

  files <- list.files(file.path(results_dir, paste0("max", max_cos)))
  files <- files[grepl(".rds", files)]
  files <- files[grepl("rank", files)]
  for (file in files) {
    print(file)
    file_parts <- str_split(str_replace(file, ".rds", ""), "_")[[1]]
    N <- file_parts[2]
    G <- file_parts[3]
    model <- file_parts[4]
    rep <- file_parts[5]
    n <- str_replace(file_parts[6], "n", "")

    res <- readRDS(file.path(results_dir, paste0("max", max_cos), file))
    if (is.null(res$converged_at)) {
      next()
    }
    data <- readRDS(file.path(
      data_dir,
      paste0("max", max_cos),
      paste(N, G, rep, "cosmic.rds", sep = "_")
    ))

    # rescale results
    Phat <- res$MAP$P
    if (sum(res$MAP$A) == 1) {
      Ehat <- matrix(res$MAP$E * sum(Phat), nrow = 1)
      Phat <- matrix(Phat / sum(Phat), ncol = 1)
    } else {
      Ehat <- sweep(
        res$MAP$E,
        1,
        colSums(Phat),
        '*'
      )
      Phat <- sweep(
        Phat,
        2,
        colSums(Phat),
        '/'
      )
    }
    Mhat <- Phat %*% Ehat
    if (model %in% c("PG", "PE")) {
      res$MAP$sigmasq <- colMeans(Mhat)
    }

    # metrics
    P_similarity_matrix <- pairwise_sim(data$P, Phat, which = "cols")
    P_hung <- hungarian_algorithm(
      P_similarity_matrix, which = "max"
    )
    P_aligned_similarity_matrix <- P_hung$mat

    E_similarity_matrix <- pairwise_sim(data$E, Ehat, which = "rows")
    E_hung <- hungarian_algorithm(
      E_similarity_matrix, which = "max"
    )
    E_aligned_similarity_matrix <- E_hung$mat

    this_metrics <- list(
      file = file,
      N = N,
      G = G,
      rep = rep,
      n = n,
      model = model,
      RMSE = sqrt(mean((data$M - Mhat)**2)),
      KL = sum(data$M * log(data$M / Mhat) - data$M + Mhat),
      P_mean_sim = mean(diag(P_aligned_similarity_matrix)),
      P_median_sim = median(diag(P_aligned_similarity_matrix)),
      P_max_sim = max(diag(P_aligned_similarity_matrix)),
      P_min_sim = min(diag(P_aligned_similarity_matrix)),
      E_mean_sim = mean(diag(E_aligned_similarity_matrix)),
      E_median_sim = median(diag(E_aligned_similarity_matrix)),
      E_max_sim = max(diag(E_aligned_similarity_matrix)),
      E_min_sim = min(diag(E_aligned_similarity_matrix)),
      iters = res$converged_at,
      total_seconds = res$time$total_secs,
      secs_per_iter = res$time$avg_secs_per_iter,
      estimated_N = sum(res$MAP$A),
      precision = calc_precision(Phat, data$P),
      sensitivity = calc_sensitivity(Phat, data$P),
      recovery = grepl("recovery", file),
      mean_sigmasq = mean(res$MAP$sigmasq),
      mean_M = mean(data$M),
      normal_BIC = calc_normal_BIC(data$M, Mhat, res$MAP$sigmasq, as.numeric(n)),
      poisson_BIC = calc_poisson_BIC(data$M, Mhat, as.numeric(n))
    )

    metrics <- rbind(metrics, unlist(this_metrics))
    colnames(metrics) <- names(this_metrics)
  }
  return(metrics)
}