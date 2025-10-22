
#' Compare credible intervals of two samplers with paired t-test between CI widths
#' separately for P and E matrices
#'
#' @param sampler1 bayesNMF_sampler object
#' @param sampler2 bayesNMF_sampler object
#' @param name1 string, name of the first sampler
#' @param name2 string, name of the second sampler
#'
#' @return data.frame with correlation and paired t-test statistics and p-values
#' @export
compare_CIs <- function(sampler1, sampler2, name1, name2) {
  # align sampler2 to sampler1
  assignments <- hungarian_assignment(sampler2$MAP$P, reference_P = sampler1$MAP$P)

  widths_P_1 <- c(sampler1$credible_intervals$P[[2]] - sampler1$credible_intervals$P[[1]])
  widths_P_2 <- c(sampler2$credible_intervals$P[[2]] - sampler2$credible_intervals$P[[1]])
  widths_P_2 <- widths_P_2[assignments$sig_ref]
  widths_E_1 <- c(sampler1$credible_intervals$E[[2]] - sampler1$credible_intervals$E[[1]])
  widths_E_2 <- c(sampler2$credible_intervals$E[[2]] - sampler2$credible_intervals$E[[1]])
  widths_E_2 <- widths_E_2[assignments$sig_ref]

  P_test = t.test(widths_P_1, widths_P_2, paired = TRUE)
  E_test = t.test(widths_E_1, widths_E_2, paired = TRUE)

  data.frame(
    "type" = c("P", "E"),
    "cor" = c(cor(widths_P_1, widths_P_2), cor(widths_E_1, widths_E_2)),
    "paired_t_test_statistic" = c(P_test$statistic, E_test$statistic),
    "paired_t_test_p_value" = c(P_test$p.value, E_test$p.value)
  )
}

#' Plot a comparison of credible intervals of two samplers
#'
#' @param sampler1 bayesNMF_sampler object
#' @param sampler2 bayesNMF_sampler object
#' @param name1 string, name of the first sampler
#' @param name2 string, name of the second sampler
#'
#' @return ggplot2 object
#' @export
plot_CI_comparison <- function(sampler1, sampler2, name1, name2) {
  # align sampler2 to sampler1
  assignments <- hungarian_assignment(sampler2$MAP$P, reference_P = sampler1$MAP$P)

  widths_P_1 <- c(sampler1$credible_intervals$P[[2]] - sampler1$credible_intervals$P[[1]])
  widths_P_1 <- widths_P_1[assignments$sig_ref]
  widths_P_2 <- c(sampler2$credible_intervals$P[[2]] - sampler2$credible_intervals$P[[1]])
  
  widths_E_1 <- c(sampler1$credible_intervals$E[[2]] - sampler1$credible_intervals$E[[1]])
  widths_E_1 <- widths_E_1[assignments$sig_ref]
  widths_E_2 <- c(sampler2$credible_intervals$E[[2]] - sampler2$credible_intervals$E[[1]])
  
  cor_P = cor(widths_P_1, widths_P_2)
  test_P = t.test(widths_P_1, widths_P_2, paired = TRUE)
  title_P = paste0(
    "P matrix credible intervals:\n", 
    "Correlation = ", round(cor_P, 2), "\n",
    "Paired t-test p-value = ", round(test_P$p.value, 3)
  )

  cor_E = cor(widths_E_1, widths_E_2)
  test_E = t.test(widths_E_1, widths_E_2, paired = TRUE)
  title_E = paste0(
    "E matrix credible intervals:\n", 
    "Correlation = ", round(cor_E, 2), "\n",
    "Paired t-test p-value = ", round(test_E$p.value, 3)
  )

  data.frame(
    widths_1 = c(widths_P_1, widths_E_1),
    widths_2 = c(widths_P_2, widths_E_2),
    type = c(rep("P", length(widths_P_1)), rep("E", length(widths_E_1)))
  ) %>%
  dplyr::mutate(
    type = factor(type, levels = c("P", "E"), labels = c(title_P, title_E))
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = widths_1, y = widths_2)) +
  ggplot2::facet_wrap(ggplot2::vars(type), scales = "free") +
  ggplot2::geom_abline(slope = 1, intercept = 0, col = 'blue') +
  ggplot2::geom_point() +
  ggplot2::labs(
    x = paste0(name1, " width"),
    y = paste0(name2, " width"),
    title = paste0("Credible interval widths for ", name1, " vs ", name2)
  ) +
  ggplot2::theme(
    text = ggplot2::element_text(size = 20)
  )
}

#' Compare MAPs of two samplers
#'
#' @param sampler1 bayesNMF_sampler object
#' @param sampler2 bayesNMF_sampler object
#' @param name1 string, name of the first sampler
#' @param name2 string, name of the second sampler
#'
#' @return data.frame with correlation and paired t-test statistics and p-values
#' @export
compare_MAPs <- function(sampler1, sampler2, name1, name2) {
  # align sampler2 to sampler1
  assignments <- hungarian_assignment(sampler2$MAP$P, reference_P = sampler1$MAP$P)

  MAP_E_1 <- c(sampler1$MAP$E[assignments$sig_ref, ])
  MAP_E_2 <- c(sampler2$MAP$E)
  cor_E = cor(MAP_E_1, MAP_E_2)
  test_E = t.test(MAP_E_1, MAP_E_2, paired = TRUE)
  cosine_E = pairwise_sim(sampler1$MAP$E[assignments$sig_ref, ], sampler2$MAP$E, which = "rows")

  MAP_P_1 <- c(sampler1$MAP$P[, assignments$sig_ref])
  MAP_P_2 <- c(sampler2$MAP$P)
  cor_P = cor(MAP_P_1, MAP_P_2)
  test_P = t.test(MAP_P_1, MAP_P_2, paired = TRUE)
  cosine_P = pairwise_sim(sampler1$MAP$P[, assignments$sig_ref], sampler2$MAP$P, which = "cols")

  RMSE_P = sqrt(mean((MAP_P_1 - MAP_P_2)**2))
  RMSE_E = sqrt(mean((MAP_E_1 - MAP_E_2)**2))
  
  data.frame(
    "type" = c("P", "E"),
    "cor" = c(cor_P, cor_E),
    "paired_t_test_statistic" = c(test_P$statistic, test_E$statistic),
    "paired_t_test_p_value" = c(test_P$p.value, test_E$p.value),
    "RMSE" = c(RMSE_P, RMSE_E),
    "min_cosine" = c(min(diag(cosine_P)), min(diag(cosine_E)))
  )
}

#' Plot a comparison of MAPs of two samplers
#'
#' @param sampler1 bayesNMF_sampler object
#' @param sampler2 bayesNMF_sampler object
#' @param name1 string, name of the first sampler
#' @param name2 string, name of the second sampler
#'
#' @return ggplot2 object
#' @export
plot_MAP_comparison <- function(sampler1, sampler2, name1, name2) {
  # align sampler2 to sampler1
  assignments <- hungarian_assignment(sampler2$MAP$P, reference_P = sampler1$MAP$P)
  
  MAP_E_1 <- c(sampler1$MAP$E)
  MAP_E_2 <- c(sampler2$MAP$E[assignments$sig_ref, ])
  cor_E = cor(MAP_E_1, MAP_E_2)
  test_E = t.test(MAP_E_1, MAP_E_2, paired = TRUE)
  title_E  = paste0(
    "E matrix MAPs:\n", 
    "Correlation = ", round(cor_E, 2), "\n",
    "Paired t-test p-value = ", round(test_E$p.value, 3)
  )

  MAP_P_1 <- c(sampler1$MAP$P)
  MAP_P_2 <- c(sampler2$MAP$P[, assignments$sig_ref])
  cor_P = cor(MAP_P_1, MAP_P_2)
  test_P = t.test(MAP_P_1, MAP_P_2, paired = TRUE)
  title_P = paste0(
    "P matrix MAPs:\n",
    "Correlation = ", round(cor_P, 2), "\n",
    "Paired t-test p-value = ", round(test_P$p.value, 3)
  )

  data.frame(
    MAP_1 = c(MAP_E_1, MAP_P_1),
    MAP_2 = c(MAP_E_2, MAP_P_2),
    type = c(rep("E", length(MAP_E_1)), rep("P", length(MAP_P_1)))
  ) %>%
  dplyr::mutate(
    type = factor(type, levels = c("E", "P"), labels = c(title_E, title_P))
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = MAP_1, y = MAP_2)) +
  ggplot2::facet_wrap(ggplot2::vars(type), scales = "free") +
  ggplot2::geom_abline(slope = 1, intercept = 0, col = 'blue') +
  ggplot2::geom_point() +
  ggplot2::labs(
    x = paste0(name1, " MAP"),
    y = paste0(name2, " MAP"),
    title = paste0("MAPs for ", name1, " vs ", name2)
  ) +
  ggplot2::theme(text = ggplot2::element_text(size = 20))
}

process_compare_MAPs <- function(prior, compare = c("Normal", "Poisson")) {
  prefix <- str_sub(compare, 1, 1)
  subscript <- sapply(str_split(compare, "_"), function(st) {ifelse(length(st) > 1, paste0("_", st[2]), '')})
  compare_short <- paste0(prefix, subscript)
  if (prior == "exponential") {
    dirs = paste(prefix, 'E', sep = '_')
    dirs = paste0(dirs, subscript)
    names(dirs) = compare
  } else {
    dirs = sapply(1:2, function(i) {
      if (startsWith(compare[i], "N") | grepl("MH", compare[i])) {
        return(paste(prefix[i], 'T', sep = '_'))
      } else {
        return(paste(prefix[i], 'G', sep = '_'))
      }
    })
    dirs = paste0(dirs, subscript)
    names(dirs) = compare
  }

  results <- lapply(
    list.dirs(file.path(output_dir, dirs[1]), full.names = FALSE, recursive = FALSE),
    function(dir) {
      sampler_path_1 <- file.path(output_dir, dirs[compare[1]], dir, "sampler.rds")
      sampler_path_2 <- file.path(output_dir, dirs[compare[2]], dir, "sampler.rds")

      if (!file.exists(sampler_path_1) || !file.exists(sampler_path_2)) {
        return(NULL)
      }
      print(dir)

      sampler_1 <- readRDS(sampler_path_1)
      sampler_2 <- readRDS(sampler_path_2)
      MAP_comparison <- compare_MAPs(sampler_1, sampler_2, compare[1], compare[2])

      results_name_parts <- str_split(dir, "_")[[1]]
      N <- results_name_parts[1] %>% str_remove("N") %>% as.numeric()
      G <- results_name_parts[2] %>% str_remove("G") %>% as.numeric()
      rep <- results_name_parts[3] %>% str_remove("rep") %>% as.numeric()

      data <- readRDS(file.path(data_dir, glue("N{N}_G{G}_rep{rep}.rds")))

      list(
        'N' = N,
        'G' = G,
        'rep' = rep,
        'prior' = prior,
        'compare' = paste(compare, collapse = "<vs>"),
        'P_cor' = MAP_comparison %>% filter(type == 'P') %>% pull(cor),
        'E_cor' = MAP_comparison %>% filter(type == 'E') %>% pull(cor),
        'P_RMSE' = MAP_comparison %>% filter(type == 'P') %>% pull(RMSE),
        'E_RMSE' = MAP_comparison %>% filter(type == 'E') %>% pull(RMSE),
        'P_t_test_statistic' = MAP_comparison %>% 
          filter(type == 'P') %>% 
          pull(paired_t_test_statistic),
        'E_t_test_statistic' = MAP_comparison %>% 
          filter(type == 'E') %>% 
          pull(paired_t_test_statistic),
        'P_t_test_p_value' = MAP_comparison %>% 
          filter(type == 'P') %>% 
          pull(paired_t_test_p_value),
        'E_t_test_p_value' = MAP_comparison %>% 
          filter(type == 'E') %>% 
          pull(paired_t_test_p_value),
        'P_min_cosine' = MAP_comparison %>%
          filter(type == 'P') %>% pull(min_cosine),
        'E_min_cosine' = MAP_comparison %>%
          filter(type == 'E') %>% pull(min_cosine),
        "mean_mut_per_sample" = mean(colMeans(data$M))
      )
    }
  )

  results <- results %>%
    do.call(rbind, .) %>%
    data.frame() %>%
    dplyr::mutate(across(everything(), as.character))

  return(results)
}

compare_credible_intervals <- function(sampler_1, sampler_2, data) {
  assignments <- hungarian_assignment(sampler_2$MAP$P, reference_P = sampler_1$MAP$P)

  width_P_1 <- sampler_1$credible_intervals$P[[2]] - sampler_1$credible_intervals$P[[1]]
  width_P_1 <- width_P_1[, assignments$sig_ref] # reorder to match sampler2's order
  width_P_2 <- sampler_2$credible_intervals$P[[2]] - sampler_2$credible_intervals$P[[1]]
  width_ratio_P <- c(width_P_1) / c(width_P_2)
  width_diff_P <- c(width_P_1) - c(width_P_2)
  
  width_E_1 <- sampler_1$credible_intervals$E[[2]] - sampler_1$credible_intervals$E[[1]]
  width_E_1 <- width_E_1[assignments$sig_ref, ] # reorder to match sampler2's order
  width_E_2 <- sampler_2$credible_intervals$E[[2]] - sampler_2$credible_intervals$E[[1]]
  width_ratio_E <- c(width_E_1) / c(width_E_2)
  width_diff_E <- c(width_E_1) - c(width_E_2)

  upper_E_prop_1 <- broadcast(
    colSums(data), across = "rows", 
    of = sampler_1$credible_intervals$E[[2]], with = "/"
  )
  lower_E_prop_1 <- broadcast(
    colSums(data), across = "rows", 
    of = sampler_1$credible_intervals$E[[1]], with = "/"
  )
  width_E_prop_1 <- upper_E_prop_1 - lower_E_prop_1

  upper_E_prop_2 <- broadcast(
    colSums(data), across = "rows", 
    of = sampler_2$credible_intervals$E[[2]], with = "/"
  )
  lower_E_prop_2 <- broadcast(
    colSums(data), across = "rows", 
    of = sampler_2$credible_intervals$E[[1]], with = "/"
  )
  width_E_prop_2 <- upper_E_prop_2 - lower_E_prop_2

  width_ratio_E_prop <- width_E_prop_1 / width_E_prop_2
  width_diff_E_prop <- width_E_prop_1 - width_E_prop_2

  data.frame(
    which = c('P','E','E_prop'),
    mean_width_ratio = c(mean(width_ratio_P), mean(width_ratio_E), mean(width_ratio_E_prop)),
    median_width_ratio = c(median(width_ratio_P), median(width_ratio_E), median(width_ratio_E_prop)),
    q10_width_ratio = c(quantile(width_ratio_P, 0.1), quantile(width_ratio_E, 0.1), quantile(width_ratio_E_prop, 0.1)),
    q90_width_ratio = c(quantile(width_ratio_P, 0.9), quantile(width_ratio_E, 0.9), quantile(width_ratio_E_prop, 0.9)),
    q2.5_width_ratio = c(quantile(width_ratio_P, 0.025), quantile(width_ratio_E, 0.025), quantile(width_ratio_E_prop, 0.025)),
    q97.5_width_ratio = c(quantile(width_ratio_P, 0.975), quantile(width_ratio_E, 0.975), quantile(width_ratio_E_prop, 0.975)),
    min_width_ratio = c(min(width_ratio_P ), min(width_ratio_E), min(width_ratio_E_prop)),
    max_width_ratio = c(max(width_ratio_P), max(width_ratio_E), max(width_ratio_E_prop)),
    mean_width_diff = c(mean(width_diff_P), mean(width_diff_E), mean(width_diff_E_prop)),
    median_width_diff = c(median(width_diff_P), median(width_diff_E), median(width_diff_E_prop)),
    q10_width_diff = c(quantile(width_diff_P, 0.1), quantile(width_diff_E, 0.1), quantile(width_diff_E_prop, 0.1)),
    q90_width_diff = c(quantile(width_diff_P, 0.9), quantile(width_diff_E, 0.9), quantile(width_diff_E_prop, 0.9)),
    q2.5_width_diff = c(quantile(width_diff_P, 0.025), quantile(width_diff_E, 0.025), quantile(width_diff_E_prop, 0.025)),
    q97.5_width_diff = c(quantile(width_diff_P, 0.975), quantile(width_diff_E, 0.975), quantile(width_diff_E_prop, 0.975)),
    min_width_diff = c(min(width_diff_P), min(width_diff_E), min(width_diff_E_prop)),
    max_width_diff = c(max(width_diff_P), max(width_diff_E), max(width_diff_E_prop))
  )
}
process_compare_credible_intervals <- function(prior, compare = c("Poisson", "Poisson_MH")) {
  prefix <- str_sub(compare, 1, 1)
  subscript <- sapply(str_split(compare, "_"), function(st) {ifelse(length(st) > 1, paste0("_", st[2]), '')})
  compare_short <- paste0(prefix, subscript)
  if (prior == "exponential") {
    dirs = paste(prefix, 'E', sep = '_')
    dirs = paste0(dirs, subscript)
    names(dirs) = compare
  } else {
    dirs = sapply(1:2, function(i) {
      if (startsWith(compare[i], "N") | grepl("MH", compare[i])) {
        return(paste(prefix[i], 'T', sep = '_'))
      } else {
        return(paste(prefix[i], 'G', sep = '_'))
      }
    })
    dirs = paste0(dirs, subscript)
    names(dirs) = compare
  }

  results <- lapply(
    list.dirs(file.path(output_dir, dirs[1]), full.names = FALSE, recursive = FALSE),
    function(dir) {
      sampler_path_1 <- file.path(output_dir, dirs[compare[1]], dir, "sampler.rds")
      sampler_path_2 <- file.path(output_dir, dirs[compare[2]], dir, "sampler.rds")
      print(paste(sampler_path_1, sampler_path_2))

      if (!file.exists(sampler_path_1) || !file.exists(sampler_path_2)) {
        return(NULL)
      }
      print(dir)

      data <- readRDS(file.path(data_dir, glue("{dir}.rds")))

      sampler_1 <- readRDS(sampler_path_1)
      sampler_2 <- readRDS(sampler_path_2)
      credible_interval_comparison <- compare_credible_intervals(sampler_1, sampler_2, data$M)

      results_name_parts <- str_split(dir, "_")[[1]]
      N <- results_name_parts[1] %>% str_remove("N") %>% as.numeric()
      G <- results_name_parts[2] %>% str_remove("G") %>% as.numeric()
      rep <- results_name_parts[3] %>% str_remove("rep") %>% as.numeric()


      list(
        'N' = N,
        'G' = G,
        'rep' = rep,
        'prior' = prior,
        'compare' = paste(compare, collapse = "<vs>"),
        'P_mean_width_ratio' = credible_interval_comparison %>% filter(which == 'P') %>% pull(mean_width_ratio),
        'E_mean_width_ratio' = credible_interval_comparison %>% filter(which == 'E') %>% pull(mean_width_ratio),
        'E_mean_width_ratio_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(mean_width_ratio),
        'P_median_width_ratio' = credible_interval_comparison %>% filter(which == 'P') %>% pull(median_width_ratio),
        'E_median_width_ratio' = credible_interval_comparison %>% filter(which == 'E') %>% pull(median_width_ratio),
        'E_median_width_ratio_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(median_width_ratio),
        'P_q10_width_ratio' = credible_interval_comparison %>% filter(which == 'P') %>% pull(q10_width_ratio),
        'E_q10_width_ratio' = credible_interval_comparison %>% filter(which == 'E') %>% pull(q10_width_ratio),
        'E_q10_width_ratio_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(q10_width_ratio),
        'P_q90_width_ratio' = credible_interval_comparison %>% filter(which == 'P') %>% pull(q90_width_ratio),
        'E_q90_width_ratio' = credible_interval_comparison %>% filter(which == 'E') %>% pull(q90_width_ratio),
        'E_q90_width_ratio_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(q90_width_ratio),
        'P_q2.5_width_ratio' = credible_interval_comparison %>% filter(which == 'P') %>% pull(q2.5_width_ratio),
        'E_q2.5_width_ratio' = credible_interval_comparison %>% filter(which == 'E') %>% pull(q2.5_width_ratio),
        'E_q2.5_width_ratio_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(q2.5_width_ratio),
        'P_q97.5_width_ratio' = credible_interval_comparison %>% filter(which == 'P') %>% pull(q97.5_width_ratio),
        'E_q97.5_width_ratio' = credible_interval_comparison %>% filter(which == 'E') %>% pull(q97.5_width_ratio),
        'E_q97.5_width_ratio_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(q97.5_width_ratio),
        'P_min_width_ratio' = credible_interval_comparison %>% filter(which == 'P') %>% pull(min_width_ratio),
        'E_min_width_ratio' = credible_interval_comparison %>% filter(which == 'E') %>% pull(min_width_ratio),
        'E_min_width_ratio_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(min_width_ratio),
        'P_max_width_ratio' = credible_interval_comparison %>% filter(which == 'P') %>% pull(max_width_ratio),
        'E_max_width_ratio' = credible_interval_comparison %>% filter(which == 'E') %>% pull(max_width_ratio),
        'E_max_width_ratio_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(max_width_ratio),
        'P_mean_width_diff' = credible_interval_comparison %>% filter(which == 'P') %>% pull(mean_width_diff),
        'E_mean_width_diff' = credible_interval_comparison %>% filter(which == 'E') %>% pull(mean_width_diff),
        'E_mean_width_diff_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(mean_width_diff),
        'P_median_width_diff' = credible_interval_comparison %>% filter(which == 'P') %>% pull(median_width_diff),
        'E_median_width_diff' = credible_interval_comparison %>% filter(which == 'E') %>% pull(median_width_diff),
        'E_median_width_diff_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(median_width_diff),
        'P_q10_width_diff' = credible_interval_comparison %>% filter(which == 'P') %>% pull(q10_width_diff),
        'E_q10_width_diff' = credible_interval_comparison %>% filter(which == 'E') %>% pull(q10_width_diff),
        'E_q10_width_diff_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(q10_width_diff),
        'P_q90_width_diff' = credible_interval_comparison %>% filter(which == 'P') %>% pull(q90_width_diff),
        'E_q90_width_diff' = credible_interval_comparison %>% filter(which == 'E') %>% pull(q90_width_diff),
        'E_q90_width_diff_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(q90_width_diff),
        'P_q2.5_width_diff' = credible_interval_comparison %>% filter(which == 'P') %>% pull(q2.5_width_diff),
        'E_q2.5_width_diff' = credible_interval_comparison %>% filter(which == 'E') %>% pull(q2.5_width_diff),
        'E_q2.5_width_diff_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(q2.5_width_diff),
        'P_q97.5_width_diff' = credible_interval_comparison %>% filter(which == 'P') %>% pull(q97.5_width_diff),
        'E_q97.5_width_diff' = credible_interval_comparison %>% filter(which == 'E') %>% pull(q97.5_width_diff),
        'E_q97.5_width_diff_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(q97.5_width_diff),
        'P_min_width_diff' = credible_interval_comparison %>% filter(which == 'P') %>% pull(min_width_diff),
        'E_min_width_diff' = credible_interval_comparison %>% filter(which == 'E') %>% pull(min_width_diff),
        'E_min_width_diff_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(min_width_diff),
        'P_max_width_diff' = credible_interval_comparison %>% filter(which == 'P') %>% pull(max_width_diff),
        'E_max_width_diff' = credible_interval_comparison %>% filter(which == 'E') %>% pull(max_width_diff),
        'E_max_width_diff_prop' = credible_interval_comparison %>% filter(which == 'E_prop') %>% pull(max_width_diff),
        "mean_mut_per_sample" = mean(colMeans(data$M))
      )
    }
  )

  results <- results %>%
    do.call(rbind, .) %>%
    data.frame() %>%
    dplyr::mutate(across(everything(), as.character))

  return(results)
}
