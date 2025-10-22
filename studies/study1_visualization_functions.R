pair_results <- function(results, metric) {
  out <- results %>%
    mutate(
      prior = case_when(
        prior == "E" ~ "Exponential",
        prior %in% c("T", "G") ~ "Truncated Normal / Gamma"
      ),
      likelihood = case_when(
        likelihood == "P" ~ "Poisson",
        likelihood == "N" ~ "Normal"
      )
    ) 

    # pair MH and non-MH results with the same Normal results
    out <- rbind(
      out,
      out %>%
        filter(likelihood == "Normal") %>%
        mutate(MH = TRUE)
    )
    
    out <- out%>%
      dplyr::select(
      N, G, rep, likelihood, prior, MH, !!rlang::sym(metric)
      ) %>%
      pivot_wider(names_from = likelihood, values_from = !!rlang::sym(metric))
    return(out)
}

# plot metrics Normal (x-axis) vs Poisson (y-axis), pairing N, G, rep, and prior
plot_normal_v_poisson <- function(results, metric, log_scale = FALSE, diagonal = TRUE) {
  plot <- pair_results(results, metric) %>%
    filter(!is.na(Poisson), !is.na(Normal)) %>%
    ggplot(aes(x = Normal, y = Poisson, color = MH, shape = prior))

  if (diagonal) {
    plot <- plot +
      geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed")
  }

  plot <- plot +
    geom_point(size = 2) +
    theme_minimal() +
    theme(
        text = element_text(size = 16)
    )
  if (log_scale) {
    plot <- plot +
      scale_x_log10() +
      scale_y_log10()
  }
  return(plot)
}

# plot ratio of Normal to Poisson (y-axis) vs G (x-axis), faceted by N, pairing rep and prior
plot_ratio_normal_v_poisson <- function(results, metric, log_scale = FALSE) {
  plot <- pair_results(results, metric) %>%
    mutate(ratio = Poisson / Normal) %>%
    ggplot(aes(x = G, y = ratio, color = MH, shape = prior)) +
    facet_grid(rows = vars(N)) +
    geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
    geom_point(size = 3) +
    theme_minimal()
  if (log_scale) {
    plot <- plot +
      scale_x_log10() +
      scale_y_log10()
  }
  return(plot)
}