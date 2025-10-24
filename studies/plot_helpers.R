get_idx_annotations <- function(sampler, idx) {
  done_temp <- which(sampler$temperature_schedule == 1)[1]
  line_df <- data.frame(
    x = c(
      ifelse(sampler$state$converged, sampler$state$converged_iter, NA_integer_),
      ifelse(sampler$specs$learning_rank & sampler$state$iter >= done_temp, done_temp, NA_integer_)
    ),
    what = c("Convergence", "Tempering done")
  ) %>% 
    tidyr::drop_na()

  rect_df <- data.frame(
    xmin = sampler$state$iter - sampler$specs$convergence_control$MAP_over,
    xmax = sampler$state$iter,
    what = "Inference"
  )
  used_levels <- unique(c(as.character(line_df$what), if (exists("rect_df")) as.character(rect_df$what)))
  used_levels <- used_levels[!is.na(used_levels)]
  colors <- c(
    "Tempering done" = "black",
    "Convergence" = "black",
    "Inference" = "black"
  )
  colors <- colors[used_levels]
  return(list(line_df = line_df, rect_df = rect_df, colors = colors, used_levels = used_levels))
}

add_annotations <- function(p_main, sampler, idx_annotations, idx, size = 10, title = NULL) {
  main_range <- ggplot2::ggplot_build(p_main)$layout$panel_params[[1]]$x.range

  if (sampler$specs$learning_rank) {
    seg <- data.frame(
      xmin = c(
        1,
        idx_annotations$line_df$x[idx_annotations$line_df$what == "Convergence"],
        idx_annotations$rect_df$xmin
      ),
      xmax = c(
        idx_annotations$line_df$x[idx_annotations$line_df$what == "Tempering done"],
        sampler$state$iter,
        idx_annotations$rect_df$xmax
      ),
      y = 1,
      eta = c(0.01, 0.01, 0.03),
      lab = c("Tempering", "MH Samples", "Inference"),
      color = c("Tempering label", "Convergence label", "Inference label"),
      placement = c(0.5, 0.25, 0.5)
    )
  } else {
    seg <- data.frame(
      xmin = c(
        idx_annotations$line_df$x[idx_annotations$line_df$what == "Convergence"],
        idx_annotations$rect_df$xmin
      ),
      xmax = c(
        sampler$state$iter,
        idx_annotations$rect_df$xmax
      ),
      y = 1,
      eta = c(0.01, 0.03),
      lab = c("MH Samples", "Inference"),
      color = c("Convergence label", "Inference label"),
      placement = c(0.25, 0.5)
    )
  }

  # remove MH Samples if !MH
  if (!sampler$specs$MH) {
    seg <- seg %>%
      dplyr::filter(lab != "MH Samples")
  }

  # only plot as endpoints if they are within idx
  seg <- seg %>%
    dplyr::mutate(
      end_min = xmin >= min(idx),
      end_max = xmax <= max(idx),
      xmin = pmax(min(idx), xmin),
      xmax = pmin(max(idx), xmax)
    )
  
  idx_annotations$colors['Inference label'] = 'black'
  idx_annotations$colors['Convergence label'] = 'black'
  idx_annotations$colors['Tempering label'] = 'black'

  p_anno <- ggplot2::ggplot(seg, ggplot2::aes(color = color)) +
    # horizontal
    ggplot2::geom_segment(ggplot2::aes(x = xmin, xend = xmax, y = 1 + eta, yend = 1 + eta)) +
    # tips
    ggplot2::geom_segment(ggplot2::aes(x = xmin, xend = xmin, y = ifelse(end_min, 1, 1 + eta), yend = 1 + eta)) +
    ggplot2::geom_segment(ggplot2::aes(x = xmax, xend = xmax, y = ifelse(end_max, 1, 1 + eta), yend = 1 + eta)) +
    # labels
    ggplot2::geom_text(ggplot2::aes(
      x = xmin + (xmax - xmin) * placement,
      y = 1.003 + eta,
      label = lab
    ), vjust = 0, size = size) +
    # match x-axis range of main plot
    ggplot2::coord_cartesian(xlim = c(main_range[1], main_range[2])) +
    ggplot2::scale_y_continuous(limits = c(1, 1.05), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 5),
      legend.position = 'none',
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
    ) +
    ggplot2::scale_color_manual(
      values = idx_annotations$colors[seg$color],
      breaks = names(idx_annotations$colors[seg$color])
    )

  # align limits and margins
  p_main2 <- p_main +
    ggplot2::scale_x_continuous(limits = main_range, expand = c(0,0)) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5))

  p_anno2 <- p_anno +
    ggplot2::scale_x_continuous(limits = main_range, expand = c(0,0)) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 5, r = 5, b = 0, l = 5))
    if (!is.null(title)) {
    p_anno2 <- p_anno2 +
        ggplot2::ggtitle(title)
    }


  aligned <- cowplot::align_plots(p_anno2, p_main2, align = "v", axis = "lr")
  
  # combine plots
  p <- cowplot::plot_grid(
    aligned[[1]], aligned[[2]], ncol = 1,
    rel_heights = c(0.05, 1)
  )
  return(p)
}