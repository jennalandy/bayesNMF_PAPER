library(MASS)
library(dirmult)
library(cluster)
library(ggplot2)

fit_repeated_NB_mixture <- function(
  x, n_runs = 10, metric = 'silhouette',
  K = 2, max_iter = 1000, tol = 1e-6
) {
  if (metric == 'silhouette' && K == 1) {
    metric <- "bic"
    print("Cant use silhouette for K = 1, using bic")
  }

  first = TRUE
  for (run in 1:n_runs) {
    out <- fit_NB_mixture(
      x, K, max_iter, tol
    )
    this <- out$model[[metric]]
    if (first) {
      first = FALSE
      best <- this
      best_out <- out
    }
    if (metric %in% c('aic','bic')) {
      if (this < best) {
        best <- this
        best_out <- out
      }
    } else {
      if (this > best) {
        best <- this
        best_out <- out
      }
    }
  }
  return(best_out)
}

# mean = mu
# prob = size / (size + mu)
# var = mu + mu^2/size
fit_NB_mixture <- function(
  x, K = 2, max_iter = 1000, tol = 1e-6,
  mu = sample(x, K), size = runif(K, 0, 5), 
  pi = dirmult::rdirichlet(1, rep(1, K))
) {
  n = length(x)
  
  loglik_prev <- -Inf
  logliks <- c()
  mus <- matrix(nrow = 0, ncol = K)
  pis <- matrix(nrow = 0, ncol = K)
  sizes <- matrix(nrow = 0, ncol = K)
  for (iter in 1:max_iter) {
    log_probs <- matrix(nrow = n, ncol = K)
    for (k in 1:K) {
      log_probs[,k] <- log(pi[k]) + dnbinom(
        x, size = size[k], mu = mu[k], log = TRUE
      )
    }
    
    max_log_probs <- apply(log_probs, 1, max)
    probs <- exp(log_probs - max_log_probs)
    probs <- probs / rowSums(probs)
    
    for (k in 1:K) {
      pi[k] <- mean(probs[,k])
      mu[k] <- sum(probs[,k] * x)/sum(probs[,k])
      var_k <- sum(probs[, k] * (x - mu[k])^2) / sum(probs[, k])
      size[k] <- mu[k]^2 / (var_k - mu[k])
      if (size[k] <= 0.01) {size[k] <- 0.01}
    }
    mus <- rbind(mus, mu)
    pis <- rbind(pis, pi)
    sizes <- rbind(sizes, size)
    
    log_probs <- matrix(nrow = n, ncol = K)
    for (k in 1:K) {
      log_probs[,k] <- log(pi[k]) + dnbinom(
        x, size = size[k], mu = mu[k], log = TRUE
      )
    }
    loglik <- sum(log(rowSums(exp(log_probs))))
    logliks <- c(logliks, loglik)
    if (abs(loglik - loglik_prev) < tol) {
      break
    }
    loglik_prev <- loglik
  }
  
  n_params <- K * 2 + (K-1)
  max_log_probs <- apply(log_probs, 1, max)
  probs <- exp(log_probs - max_log_probs)
  probs <- probs / rowSums(probs)
  clusts <- apply(probs, 1, function(row) {which(row == max(row))})
  clusts <- factor(clusts, levels = 1:K)
  if (length(unique(clusts)) > 1) {
    silhouette_scores <- cluster::silhouette(
      as.numeric(clusts), dist(x)
    )
    silhouette <- mean(silhouette_scores[, "sil_width"])
  } else {
    silhouette <- -Inf
  }
  return(list(
    logs = list(
      'mus' = mus,
      'pis' = pis,
      'sizes' = sizes,
      'logliks' = logliks
    ),
    model = list(
      'mu' = mu,
      'pi' = pi,
      'size' = size,
      'loglik' = loglik,
      'aic' = -2 * loglik + 2 * n_params,
      'bic' = -2 * loglik + log(n) * n_params,
      'silhouette' = silhouette,
      'clusters' = clusts,
      'probs' = probs,
      'x' = x
    ),
    summary = rbind(
      c('mu', mu),
      c('size', size),
      c('pi', pi),
      c('n', table(clusts)[1:K])
    )
  ))
}

dens_fn <- function(x, pi, size, mu) {
  pi * dnbinom(
    x, size = size, mu = mu
  )
}

plot_Nb_mixture <- function(out, title, ymax = NULL) {
  x = seq(0, max(out$model$x) + 100, by = 100)
  K = length(out$model$pi)
  
  ys = list()
  for(k in 1:K) {
    dens = function(x) {dens_fn(
      x, out$model$pi[k], out$model$size[k], out$model$mu[k]
    )}
    ys[[k]] <- dens(x)
  }
  
  if (is.null(ymax)) {
    ymax = max(sapply(ys, max))
  }
  
  plot(
    x, ys[[1]], type = 'l', ylim = c(0, ymax),
    col = 1, 
    xlab = '', ylab = '',
    main = glue::glue('Mixture model for {title}\nK = {K} components')
  )
  if (K > 1) {
    for (k in 2:K) {
      lines(x, ys[[k]], col = k)
    }
  }
  legend("topright", legend = 1:K, col = 1:K, lty = 1)
}

choose_K <- function(x, Ks = 1:10, metric = 'silhouette') {
  Ks = Ks[Ks < length(x)]

  K_best_bic = -1
  if (metric == 'silhouette') {
    # first round to optimize BIC
    # only continue if best BIC with K > 1
    logliks <- c()
    bics <- c()
    aics <- c()
    silhouettes <- c()
  
    for (K in Ks) {
      # best of 10 initialization for each K
      out <- fit_repeated_NB_mixture(x, K = K, metric = 'bic')
      # record metrics
      logliks <- c(logliks, out$model$loglik)
      bics <- c(bics, out$model$bic)
      aics <- c(aics, out$model$aic)
      silhouettes <- c(silhouettes, out$model$silhouette)
    }
    K_best_bic = Ks[which(bics == min(bics))]
    K_best_bic = K_best_bic[length(K_best_bic)]
  }

  if (K_best_bic != 1) { # if metric is not silhouette OR silhouette and BIC ok for 2+
    logliks <- c()
    bics <- c()
    aics <- c()
    silhouettes <- c()
    
    for (K in Ks) {
      # best of 10 initialization for each K
      out <- fit_repeated_NB_mixture(x, K = K, metric = metric)
      # record metrics
      logliks <- c(logliks, out$model$loglik)
      bics <- c(bics, out$model$bic)
      aics <- c(aics, out$model$aic)
      silhouettes <- c(silhouettes, out$model$silhouette)
    }
    
    if (metric == 'bic') {
      best = Ks[which(bics == min(bics))]
    } else if (metric == 'aic') {
      best = Ks[which(aics == min(aics))]
    } else if (metric == 'loglik') {
      best = Ks[which(logliks == max(logliks))]
    } else if (metric == 'silhouette') {
      best = Ks[which(silhouettes == max(silhouettes))]
    }
    best = best[1]
  } else {
    best = 1
  }
  

  plot = data.frame(
    K = Ks,
    loglik = logliks,
    bic = bics,
    aic = aics,
    silhouette = silhouettes
  ) %>%
    pivot_longer(2:5, names_to = 'metric', values_to = 'value') %>%
    mutate(value = ifelse(value == Inf | value == -Inf, NA, value)) %>%
    mutate(metric = case_when(
      metric == 'loglik' ~ 'Log Likelihood',
      metric == 'bic' ~ 'BIC',
      metric == 'aic' ~ 'AIC',
      metric == 'silhouette' ~ 'Silhouette Score'
    )) %>%
    ggplot(aes(x = K, y = value)) +
    facet_wrap(vars(metric), scales = 'free') +
    geom_point() +
    geom_vline(xintercept = best, color = 'blue')

  return(list(
    K = best, plot = plot
  ))
}


get_intersections <- function(f1, f2, range) {
  to_zero <- function(x) {
    f1(x) - f2(x)
  }
  x <- range[1]:range[2]
  y <- to_zero(x)
  zeros <- x[which(y[-1] * y[-length(y)] < 0)] + 0.5
  return(zeros)
}


plot_clusters <- function(out, title, bins = 30) {
  ggplot(data.frame(
    x = out$model$x,
    clust = out$model$clusters
  ), aes(x = x, fill = clust)) +
    geom_histogram(bins = bins, position = 'dodge') +
    ggtitle(glue::glue('Mixture model for {title}\nK = {K} components')) +
    labs(
      x = 'Mutation Counts per Sample',
      fill = "Cluster"
    )
}

plot_labels <- function(out, nonhyper, title, bins = 30) {
  data.frame(
    x = out$model$x,
    nonhyper = nonhyper
  ) %>%
    mutate(
      class = case_when(
        nonhyper == FALSE ~ "Hypermutated",
        TRUE ~ "Non-Hypermutated"
      )
    ) %>%
    ggplot(aes(x = x, fill = class, color = class)) +
    geom_histogram(bins = bins, position = 'dodge') +
    ggtitle(glue::glue('Mixture model for {title}\nK = {K} components')) +
    labs(
      x = 'Mutation Counts per Sample',
      fill = ""
    ) +
    geom_rug()
}