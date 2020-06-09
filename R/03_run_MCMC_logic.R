#' Wrapper for MCMC logic regression package to run monte carlo fitting of logic trees on each labeled subgroup for those found to be active binders
#'
#' This function takes in a binary matrix, assuming that the variable names include the fingerprint name such as "Pubchem"
#'
#' @param log_SVD_results Results from the SVD to get data from for plotting
#' @param receptor_type type of receptor for labeling plots
#' @param k number of groups from kmeans from after viewing elbow plot or other methods of determining groups
#' @return Plot of SVD results colored by kmeans groups
#' @export
#'


logic_MCMC <- function(outcome,
                       exposure,
                       mx_size,
                       ntrees,
                       n_burn_in,
                       niter,
                       hyperparameters,
                       outcome_permuted,
                       top_fraction)
{

  if(outcome_permuted){
    outcome <- sample(AR_pubchem_data$Label, replace = TRUE)
  }

  fit_MCMC <- logreg(resp = outcome,
                     bin=exposure,
                     type = 3,
                     select = 7,
                     nleaves = mx_size,
                     ntrees = ntrees,
                     mc.control = logreg.mc.control(nburn=n_burn_in,
                                                    niter=niter,
                                                    hyperpars=hyperparameters))

  msz_by_iter_MCMC <- fit_MCMC$size
  ni_in_model_MCMC <- fit_MCMC$single
  nij_in_model_MCMC <- as.matrix(fit_MCMC$double)
  nijk_in_model_MCMC <- as.array(fit_MCMC$triple)

  ni_in_model_MCMC_structures <- cbind(pubchem_fingerprints_cleaned, ni_in_model_MCMC)
  ni_in_model_MCMC_structures$fraction <- ni_in_model_MCMC_structures$ni_in_model_MCMC/niter
  ni_in_model_MCMC_structures_ordered <- arrange(ni_in_model_MCMC_structures, desc(fraction))

  #ij in same trees

  rownames(nij_in_model_MCMC) <- pubchem_fingerprints_cleaned[,2]
  colnames(nij_in_model_MCMC) <- pubchem_fingerprints_cleaned[,2]

  nij_in_model_MCMC_collapsed <- setNames(melt(nij_in_model_MCMC), c('Bit Structure 1',
                                                                     'Bit Structure 2',
                                                                     'count'))

  nij_in_model_MCMC_collapsed$fraction <- nij_in_model_MCMC_collapsed$count/niter
  nij_in_model_MCMC_collapsed_top <- subset(nij_in_model_MCMC_collapsed, fraction > top_fraction)
  #nij_in_model_MCMC_collapsed_top

  chem_obs <- dim(AR_pubchem_data)[1]
  pij <- colSums(exposure)/chem_obs
  names(pij) <- pubchem_fingerprints_cleaned[,2]
  pij <- as.data.frame(pij)

  pij_calcs <- merge(nij_in_model_MCMC_collapsed_top,
                     pij,
                     by.x = "Bit Structure 1",
                     by.y = "row.names")

  pij_calcs <- merge(pij_calcs,
                     pij,
                     by.x = "Bit Structure 2",
                     by.y = "row.names")

  names(pij_calcs)[names(pij_calcs) == 'pij.y'] <- 'pij Bit Structure 1'
  names(pij_calcs)[names(pij_calcs) == 'pij.x'] <- 'pij Bit Structure 2'

  pij_calcs$expected_fraction <- pij_calcs$`pij Bit Structure 1` * pij_calcs$`pij Bit Structure 2`
  pij_calcs$ratio <- pij_calcs$fraction / pij_calcs$expected_fraction

  #ijk combinations

  dimnames(nijk_in_model_MCMC)[[1]] <- pubchem_fingerprints_cleaned[,2]
  dimnames(nijk_in_model_MCMC)[[2]] <- pubchem_fingerprints_cleaned[,2]
  dimnames(nijk_in_model_MCMC)[[3]] <- pubchem_fingerprints_cleaned[,2]

  nijk_in_model_MCMC_collapsed <- setNames(melt(nijk_in_model_MCMC), c('Bit Structure 1',
                                                                       'Bit Structure 2',
                                                                       'Bit Structure 3',
                                                                       'count'))

  nijk_in_model_MCMC_collapsed$fraction <- nijk_in_model_MCMC_collapsed$count/niter
  nijk_in_model_MCMC_collapsed_top <- subset(nijk_in_model_MCMC_collapsed, fraction > top_fraction)
  #nijk_in_model_MCMC_collapsed_top

  pijk_calcs <- merge(nijk_in_model_MCMC_collapsed_top,
                      pij,
                      by.x = "Bit Structure 1",
                      by.y = "row.names")

  names(pijk_calcs)[names(pijk_calcs) == 'pij'] <- 'pijk Bit Structure 1'

  pijk_calcs <- merge(pijk_calcs,
                      pij,
                      by.x = "Bit Structure 2",
                      by.y = "row.names")

  names(pijk_calcs)[names(pijk_calcs) == 'pij'] <- 'pijk Bit Structure 2'

  pijk_calcs <- merge(pijk_calcs,
                      pij,
                      by.x = "Bit Structure 3",
                      by.y = "row.names")

  names(pijk_calcs)[names(pijk_calcs) == 'pij'] <- 'pijk Bit Structure 3'

  pijk_calcs$expected_fraction <- pijk_calcs$`pijk Bit Structure 1` *
    pijk_calcs$`pijk Bit Structure 2`*
    pijk_calcs$`pijk Bit Structure 3`

  pijk_calcs$ratio <- pijk_calcs$fraction / pijk_calcs$expected_fraction

  return(list(model = fit_MCMC,
              model_sizes = msz_by_iter_MCMC,
              vars_in_iterations = ni_in_model_MCMC_structures_ordered,
              ij_interactions = pij_calcs,
              ijk_interactions = pijk_calcs))
}
