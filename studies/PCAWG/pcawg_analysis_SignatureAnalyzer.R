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
      dimnames = list(rownames(matrix), colnames(matrix)[hungarian_alignment$X2])
    )
  }

  # return aligned matrix
  return(aligned_matrix)
}