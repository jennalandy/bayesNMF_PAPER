
library(tidyverse)

get_cosmic <- function() {
  P <- read.csv(
    "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh37.txt", # nolint: line_length_linter.
    sep = "\t"
  )
  rownames(P) <- P$Type
  P <- P[, 2:ncol(P)]
  P <- as.matrix(P)
  return(P)
}

# use hungarian algorithm to maximize sum of aligned cosine similarities
hungarian_algorithm <- function(matrix, which = "max") {
  if (which == "max") { matrix <- -1 * matrix }
  hungarian_res <- RcppHungarian::HungarianSolver(matrix)
  hungarian_alignment <- data.frame(hungarian_res$pairs) %>%
    filter(X1 != 0 & X2 != 0)

  # X1 (rows) is still in order
  stopifnot(all(hungarian_alignment$X1 == 1:nrow(hungarian_alignment)))

  # matrix with columns re-ordered
  aligned_matrix <- matrix[, hungarian_alignment$X2]
  if (which == "max") { aligned_matrix <- -1 * aligned_matrix }

  # make sure column names are still there if 1 signature
  if (nrow(matrix) == 1) {
    aligned_matrix = matrix(
      aligned_matrix, 
      dimnames = list(NULL, names(aligned_matrix))
    )
  }

  # return aligned matrix
  return(aligned_matrix)
}

signature_assignment_inference <- function(res, reference_P) {
  assignment_votes <- lapply(
    # for each posterior sample
    res$posterior_samples$P,

    # assign signatures with hungarian algorithm and record cosine similarity
    function(P_post) {
      # only assign signatures that are included and not recovery
      sigs_idx = res$MAP$A[1,] == 1 & res$final_Theta$recovery == FALSE
      P_post_discovery = P_post[,sigs_idx]
      
      if (sum(sigs_idx) == 0) {
        # no signatures to assign
        return(NULL)
      } else if (sum(sigs_idx) == 1) {
        # one signature -- make it back into a matrix
        P_post_discovery <- matrix(P_post_discovery, ncol = 1)
      }

      # compute pairwise cosine similarities and assign with hungarian algorithm
      sim <- pairwise_sim(P_post_discovery, reference_P)
      assigned <- hungarian_algorithm(sim)

      # this posterior sample "votes" for the assigned signature with weight = cos sim
      votes <- data.frame(
        cosine_sim = diag(assigned),
        sig = colnames(assigned),
        n = which(res$final_Theta$recovery[res$MAP$A[1,] == 1] == FALSE)
      )
      return(votes)
    }
  )

  if (is.null(assignment_votes[[1]])) {
    # if no signatures to assign
    assignments <- NULL
  } else {
    # total votes per reference signature for each estimated signature
    summarized_votes <- do.call(rbind, assignment_votes) %>%
      group_by(sig, n) %>%
      summarize(votes = sum(cosine_sim), .groups = "keep") %>%
      ungroup()

    assignments <- summarized_votes %>%
      dplyr::group_by(n) %>% # for each estimated signature
      dplyr::summarize(
        max_votes = max(votes), # compute max
        total = sum(votes),     # and total number of reference sig votes
        .groups = "keep"
      ) %>%
      merge(summarized_votes) %>%
      dplyr::mutate(score = votes/total) %>% # compute proportion of votes
      dplyr::filter(max_votes == votes) %>% # identify reference sig with most votes
      dplyr::arrange(n) %>% # re-order by estimated signature
      dplyr::select(sig, score, n)
  }
  
  # only if recovery-discovery -- recovery sigs get score of 1
  if (sum(res$final_Theta$recovery) > 0) {
    # only assign included reference signatures
    ref_idx = res$MAP$A[1,1:ncol(reference_P)]
    if (sum(ref_idx == 1) > 0) {
      if (is.null(assignments)) {
        assignments <- data.frame(
          sig = colnames(reference_P)[ref_idx == 1],
          score = 1,
          n = 1:sum(ref_idx == 1)
        )
      } else {
        assignments <- rbind(
          data.frame(
            sig = colnames(reference_P)[ref_idx == 1],
            score = 1,
            n = 1:sum(ref_idx == 1)
          ),
          assignments
        )
      }
    }
  }

  # compute posterior distribution cosine similarities 
  # between estimated and reference signatures
  # columns = reference sig, row = posterior samples, value = cos sim
  assigned_cos_sims <- lapply(res$posterior_samples$P, function(P_post) {
    P_post = P_post[,res$MAP$A[1,] == 1]
    sim <- pairwise_sim(P_post, reference_P)
    sapply(1:nrow(assignments), function(i) {
      sim[assignments$n[i], assignments$sig[i]]
    })
  }) %>%
    do.call(rbind, .)


  assignment_res <- list(
    'assignment' = assignments,
    'MAP'  = list(
      'cos_sim' = assigned_cos_sims %>% colMeans()
    ),
    'credible_intervals' = list(
      'cos_sim' = list(
        apply(assigned_cos_sims, 2, function(col) {quantile(col, 0.025)}),
        apply(assigned_cos_sims, 2, function(col) {quantile(col, 0.975)})
      )
    )
  )
  return(assignment_res)
}



get_results_matrix <- function(files_list) {
  total_counts <- NULL
  total_props <- NULL
  first <- TRUE
  Gs <- c()
  for (file in files_list) {
    name <- str_split(file, '\\.')[[1]][1]
    print(name)
    res <- readRDS(file.path("results", file))
    if (is.null(res$converged_at)) {
      print(paste('not done:', name))
      next()
    }
    Gs <- c(Gs, res$model$dims$G)
    assignment_res <- signature_asssignment_inference(res, cosmic)
    print(assignment_res$MAP$cos_sim)

    if (sum(res$MAP$A) == 1) {
      res$MAP$P = matrix(res$MAP$P, ncol = 1)
      res$MAP$E = matrix(res$MAP$E, nrow = 1)
    }
    colnames(res$MAP$P) <- assignment_res$assignment$sig
    rownames(res$MAP$E) <- assignment_res$assignment$sig

    res$MAP$E <- sweep(res$MAP$E, 1, colSums(res$MAP$P), '*')
    res$MAP$P <- sweep(res$MAP$P, 2, colSums(res$MAP$P), '/')

    cosmic_sig_counts <- rep(0, ncol(cosmic))
    names(cosmic_sig_counts) <- colnames(cosmic)
    this_counts <- apply(res$MAP$E, 1, median)
    cosmic_sig_counts[names(this_counts)] <- this_counts

    cosmic_cos <- rep(0, ncol(cosmic))
    names(cosmic_cos) <- colnames(cosmic)
    cosmic_cos[names(assignment_res$MAP$cos_sim)] <- assignment_res$MAP$cos_sim

    # proportion of samples with signature
    prop_contrib <- sweep(res$MAP$E, 2, colSums(res$MAP$E), '/')
    prop_with_sig <- apply(prop_contrib, 1, function(row){mean(row > 0.05)})
    names(prop_with_sig) <- rownames(res$MAP$E)
    cosmic_sig_prop <- rep(0, ncol(cosmic))
    names(cosmic_sig_prop) <- colnames(cosmic)
    cosmic_sig_prop[names(prop_with_sig)] <- prop_with_sig

    if (first) {
      total_counts <- c('Histology' = name, cosmic_sig_counts)
      total_prop <- c('Histology' = name, cosmic_sig_prop)
      total_cos <- c('Histology' = name, cosmic_cos)
      first = FALSE
    } else {
      total_counts <- rbind(
        total_counts, c('Histology' = name, cosmic_sig_counts)
      )
      total_prop <- rbind(
        total_prop, c('Histology' = name, cosmic_sig_prop)
      )
      total_cos <- rbind(
        total_cos, c('Histology' = name, cosmic_cos)
      )
    }
  }

  results <- data.frame(total_counts) %>%
    mutate(G = Gs) %>%
    pivot_longer(
      2:ncol(total_counts), 
      names_to = "Signature", 
      values_to = "Med_Contribution"
    ) %>%
    merge(
      data.frame(total_cos) %>%
        pivot_longer(
          2:ncol(total_cos),
          names_to = "Signature",
          values_to = "Cosine_Similarity"
      )
    ) %>%
    mutate(
      Med_Contribution = as.numeric(Med_Contribution),
      Signature = factor(Signature, levels = rev(colnames(cosmic))),
      Histology = str_replace_all(Histology, "_SBFI", "")
    ) %>%
    filter(Med_Contribution > 0)

  return(results)
}