#' Plot the first two PCs from logistic SVD for molecular features and color by k-means results
#'
#' This function takes in a binary matrix, assuming that the variable names include the fingerprint name such as "Pubchem"
#'
#' @param log_SVD_results Results from the SVD to get data from for plotting
#' @param receptor_type type of receptor for labeling plots
#' @param k number of groups from kmeans from after viewing elbow plot or other methods of determining groups
#' @return Plot of SVD results colored by kmeans groups
#' @export
#'
plot_SVD_kmeans <- function(receptor_type, log_SVD_results, binary_data, k = 4, endo_ligand = FALSE, endo_ligand_id, additional_label = TRUE, outcome_label = "ASSAY_OUTCOME") {
  binary_data$Label_names <- tolower(binary_data$Label_names)
  cluster_data_actives <- subset(binary_data, Label_names == "active")

  kmeans_groups <- kmeans(log_SVD_results$A, k, nstart = 50, iter.max = 15)
  groups <- kmeans_groups$cluster

  if (endo_ligand) {
    groups_endo_idx <- groups
    idx <- match(endo_ligand_id, cluster_data_actives$Name)
    groups_endo_idx[idx] <- 20
  }

  plot_w_kmeans <- plot(log_SVD_results, type = "scores") +
    geom_point(aes(colour = as.factor(groups_endo_idx))) +
    ggtitle(paste("Positive Binding for ", receptor_type, sep = ""))

  data_negatives <- binary_data %>% filter(Label_names == "inconclusive" | Label_names == "inactive")

  cluster_data_actives$log_SVM_group <- groups
  cluster_data_actives$PC_1 <- log_SVD_results$A[, 1]
  cluster_data_actives$PC_2 <- log_SVD_results$A[, 2]

  data_negatives$log_SVM_group <- 0
  data_negatives$log_SVM_group <- 0

  data_negatives$PC_1 <- NA
  data_negatives$PC_2 <- NA

  data_with_PCs_kmeans_groups <- rbind(cluster_data_actives, data_negatives)

  if (additional_label) {
    plot_outcome <- plot(log_SVD_results, type = "scores") +
      geom_point(aes(colour = as.factor(cluster_data_actives$ASSAY_OUTCOME))) +
      ggtitle(paste("Positive Binding for ", receptor_type, sep = ""))
  } else {
    plot_outcome <- NA
  }

  if (endo_ligand) {
    return(list(plot_kmeans = plot_w_kmeans, plot_outcome_label = plot_outcome, updated_data = data_with_PCs_kmeans_groups, endo_ligand_grp_included = groups_endo_idx))
  } else {
    list(plot_kmeans = plot_w_kmeans, plot_outcome_label = plot_outcome, updated_data = data_with_PCs_kmeans_groups)
  }
}
