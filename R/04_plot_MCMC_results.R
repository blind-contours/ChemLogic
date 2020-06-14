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
plot_model_size <- function(MCMC_results, receptor) {
  MCMC_model_sizes <- MCMC_results$model_sizes
  model_size_plot <- ggplot(MCMC_model_sizes, aes(x = size, y = number)) +
    geom_bar(
      stat = "identity", aes(fill = receptor),
      position = "dodge"
    ) +
    xlab("MCMC Model Sizes") +
    ylab("Count") +
    ggtitle(paste("MCMC Logic Regression for ", receptor, sep = "")) +
    theme_bw()

  return(model_size_plot)
}


plot_compare_model_size <- function(MCMC_results_1, MCMC_results_2, receptor_1, receptor_2) {
  MCMC_model_1_sizes <- MCMC_results_1$model_sizes
  MCMC_model_2_sizes <- MCMC_results_2$model_sizes

  MCMC_model_1_sizes$receptor <- receptor_1
  MCMC_model_2_sizes$receptor <- receptor_2

  combind_receptor_model_sizes <- rbind(MCMC_model_1_sizes, MCMC_model_2_sizes)

  model_size_plot <- ggplot(combind_receptor_model_sizes, aes(x = size, y = number)) +
    geom_bar(
      stat = "identity", aes(fill = receptor),
      position = "dodge"
    ) +
    xlab("MCMC Model Sizes") +
    ylab("Count") +
    ggtitle(paste("Comparing MCMC Model Sizes for ", receptor_1, " and ", receptor_2, sep = "")) +
    theme_bw()

  return(model_size_plot)
}

plot_variables_in_MCMC <- function(MCMC_results, receptor, top = 20, iters = 200000) {
  vars_in_MCMC_iters <- MCMC_results$vars_in_iterations
  vars_in_MCMC_iters$receptor <- receptor

  var_use_top <- vars_in_MCMC_iters[order(-vars_in_MCMC_iters$ni_in_model_MCMC), ]
  var_use_top <- head(var_use_top, top)

  var_use_plot <- ggplot(var_use_top) +
    geom_point(aes(
      x = fraction, y = reorder(Bit_Substructure, fraction),
      colour = factor(receptor)
    )) +
    xlab(paste("Fraction of", iters)) +
    ylab("Bit Substructure") +
    labs(colour = "Receptor") +
    ggtitle(paste("Top 20 Molecular Structures to Fit the ", receptor, "Receptor", sept = ""))

  return(var_use_plot)
}

plot_compare_vars_in_MCMC <- function(MCMC_results_1, MCMC_results_2, receptor_1, receptor_2, top = 20, iters = 200000) {
  recpt_one_vars_in_MCMC_iters <- MCMC_results_1$vars_in_iterations
  recpt_one_vars_in_MCMC_iters$receptor <- receptor_1

  recpt_two_vars_in_MCMC_iters <- MCMC_results_2$vars_in_iterations
  recpt_two_vars_in_MCMC_iters$receptor <- receptor_2

  recpt_one_var_use_top <- recpt_one_vars_in_MCMC_iters[order(-recpt_one_vars_in_MCMC_iters$ni_in_model_MCMC), ]
  var_use_top_recpt1 <- head(recpt_one_var_use_top, top)

  recpt_two_var_use_top <- recpt_two_vars_in_MCMC_iters[order(-recpt_two_vars_in_MCMC_iters$ni_in_model_MCMC), ]
  var_use_top_recpt2 <- head(recpt_two_var_use_top, top)

  rept_one_top_var_use_against_rept_two_merge <- rbind(
    var_use_top_recpt1,
    recpt_two_vars_in_MCMC_iters[recpt_two_vars_in_MCMC_iters$Bit_Substructure %in% var_use_top_recpt1$Bit_Substructure, ]
  )

  rept_two_top_var_use_against_rept_one_merge <- rbind(
    var_use_top_recpt2,
    recpt_one_vars_in_MCMC_iters[recpt_one_vars_in_MCMC_iters$Bit_Substructure %in% var_use_top_recpt2$Bit_Substructure, ]
  )


  top_rept1_vars_against_rept2_plot <- ggplot(rept_one_top_var_use_against_rept_two_merge) +
    geom_point(aes(
      x = fraction, y = reorder(Bit_Substructure, fraction),
      colour = factor(receptor)
    )) +
    xlab(paste("Fraction of ", iters, " iterations")) +
    ylab("Bit Substructure") +
    labs(colour = "Receptor") +
    ggtitle(paste("Top 20 Molecular Structures to Fit ", receptor_1))

  top_rept2_vars_against_rept1_plot <- ggplot(rept_two_top_var_use_against_rept_one_merge) +
    geom_point(aes(
      x = fraction, y = reorder(Bit_Substructure, fraction),
      colour = factor(receptor)
    )) +
    xlab(paste("Fraction of ", iters, " iterations")) +
    ylab("Bit Substructure") +
    labs(colour = "Receptor") +
    ggtitle(paste("Top 20 Molecular Structures to Fit ", receptor_2))



  return(list(
    top_receptor_1_agst_2 = top_rept1_vars_against_rept2_plot,
    top_receptor_2_agst_1 = top_rept2_vars_against_rept1_plot
  ))
}

plot_compare_vars_subgroups <- function(MCMC_results, receptor, top = 20, iters = 200000) {
  subgroup_var_grabber <- function(x, y) {
    recpt_one_vars_in_MCMC_iters <- x$vars_in_iterations
    recpt_one_vars_in_MCMC_iters$group <- y
    return(recpt_one_vars_in_MCMC_iters)
  }

  subgroup_vars_collapsed <- imap_dfr(.x = MCMC_results, .f = subgroup_var_grabber, .id = NULL)

  get_top <- function(x) {
    group <- unique(x$group)
    target_group_ordered <- x[order(-x$ni_in_model_MCMC), ]
    target_group_top <- head(target_group_ordered, top)
    subgroup_vars_collapsed_filtered <- subgroup_vars_collapsed %>% filter(group != !!(group))

    target_group_top_vs_rest <- rbind(
      target_group_top,
      subgroup_vars_collapsed_filtered[subgroup_vars_collapsed_filtered$Bit_Substructure %in% target_group_top$Bit_Substructure, ]
    )

    target_group_top_vs_rest_plot <- ggplot(target_group_top_vs_rest) +
      geom_point(aes(
        x = fraction, y = reorder(Bit_Substructure, fraction),
        colour = factor(group)
      )) +
      xlab(paste("Fraction of ", iters, " iterations")) +
      ylab("Bit Substructure") +
      labs(colour = "Receptor") +
      ggtitle(paste("Top 20 Molecular Structures to Fit Group", group, "vs Rest"))

    return(target_group_top_vs_rest_plot)
  }

  subgroup_top_targets_vs_all_plots <- by(subgroup_vars_collapsed, INDICES = subgroup_vars_collapsed[, "group"], get_top)

  return(target_subgroup_plot_comparisons = subgroup_top_targets_vs_all_plots)
}

extract_MCMC_couples_triples <- function(MCMC_results, receptor, couple_thresh, triple_thresh) {
  data_ij_interaction <- MCMC_results$ij_interactions
  data_ijk_interaction <- MCMC_results$ijk_interactions

  data_ij_interaction$receptor <- receptor
  data_ijk_interaction$receptor <- receptor

  data_ij_interactions_top <- subset(data_ij_interaction, fraction > couple_thresh)
  data_ijk_interactions_top <- subset(data_ijk_interaction, fraction > triple_thresh)


  return(list(
    MCMC_couples = data_ij_interactions_top,
    MCMC_triples = data_ijk_interactions_top
  ))
}

extract_MCMC_couples_triples_subgroups <- function(MCMC_results,
                                                   couple_thresh,
                                                   triple_thresh) {
  subgroup_couples_grabber <- function(x, y) {
    if (dim(x$ijk_interactions)[1] > 0) {
      data_ij_interaction <- x$ij_interactions
      data_ij_interaction$group <- y

      data_ij_interactions_top <- subset(data_ij_interaction, fraction > couple_thresh)
    } else {
      data_ij_interactions_top <- NULL
    }

    return(couples = data_ij_interactions_top)
  }

  subgroup_triples_grabber <- function(x, y) {
    if (dim(x$ijk_interactions)[1] > 0) {
      data_ijk_interaction <- x$ijk_interactions
      data_ijk_interaction$group <- y

      data_ijk_interactions_top <- subset(data_ijk_interaction, fraction >= triple_thresh)
    } else {
      data_ijk_interactions_top <- NULL
    }

    return(triples = data_ijk_interactions_top)
  }

  subgroup_MCMC_couples_results <- imap_dfr(
    .x = MCMC_results,
    .f = subgroup_couples_grabber
  )

  subgroup_MCMC_triples_results <- imap_dfr(
    .x = MCMC_results,
    .f = subgroup_triples_grabber
  )

  return(list(
    couples_df = subgroup_MCMC_couples_results,
    triples_df = subgroup_MCMC_triples_results
  ))
}
