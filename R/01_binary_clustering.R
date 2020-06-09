#' Run logistic singular value decomposition on the binary molecular features
#'
#' This function takes in a binary matrix, assuming that the variable names include the fingerprint name such as "Pubchem"
#'
#' @param receptor_type Used for the kmeans plot from the SVD results, this is a characer string such as "Androgen"
#' @param data matrix or dataframe that has all binary data except for possibly Name, Label_names, and Label, which are removed
#' @param rank rank of the SVD that decomposes the binary molecular feature data before using the PCs in kmeans.
#' @return Results from the logistic SVD
#' @return Plot of Kmeans for visual inspection before creating clusters
#' @export
#'
cluster_chem_data <- function(receptor_type, data, rank = 3, filter_type = "Pubchem") {
  data$Label_names <- tolower(data$Label_names)
  cluster_data_actives <- subset(data, Label_names == "active")

  cluster_data_actives_mol_features <- dplyr::select(
    cluster_data_actives,
    contains(filter_type)
  )

  logsvd_model_results <- logisticPCA::logisticSVD(cluster_data_actives_mol_features, k = rank)

  set.seed(123)
  k.max <- 15
  data <- logsvd_model_results$A
  wss <- sapply(
    1:k.max,
    function(k) {
      kmeans(data, k, nstart = 50, iter.max = 15)$tot.withinss
    }
  )

  kmeans_plot <- ggplot(data = as.data.frame(cbind(k = 1:k.max, wss)), aes(x = k, y = wss, group = 1)) +
    geom_line() +
    geom_point() +
    xlab("Number of clusters K") +
    ylab("Total within-clusters sum of squares") +
    labs(title = paste("K-means WSS for ", receptor_type, sep = ""), x = "Number of clusters K", y = "Total within-clusters sum of squares")


  return(list(logsvd_results = logsvd_model_results, kmean_plot = kmeans_plot))
}
