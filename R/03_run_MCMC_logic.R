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

run_logic_MCMC <- function(
                           outcome_label,
                           filter_type,
                           binary_data,
                           subgroups,
                           subgroup_label,
                           mx_size,
                           ntrees,
                           n_burn_in,
                           niter,
                           hyperparameters = log(2),
                           outcome_permuted = FALSE,
                           top_fraction = 0.2,
                           find_pop_tree = FALSE) {
  if (subgroups) {
    binary_data_split <- binary_data %>%
      dplyr::filter(!!(outcome_label) == 1) %>%
      dplyr::group_split((!!(subgroup_label)))

    control_data <- binary_data %>%
      dplyr::filter(!!(outcome_label) != 1)

    subgroup_dfs <- purrr::map(.x = binary_data_split, ~ rbind(.x, control_data))

    subgroup_exposures <- purrr::map(.x = subgroup_dfs, ~ dplyr::select(
      .x,
      contains(filter_type)
    ))

    subgroup_outcomes <- purrr::map(.x = subgroup_dfs, ~ dplyr::select(
      .x,
      !!(subgroup_label)
    ))

    subgroup_names <-purrr::map(.x = subgroup_dfs, ~ dplyr::select(
      .x, Name)
    )

    subgroup_outcomes <- purrr::map(.x = subgroup_outcomes, droplevels)

    check_outcome <- function(x) {
      if (is.numeric(x[[1]])) {
        x[[1]] <- ifelse(x[[1]] >= 1, 1, 0)
      }
      if (is.factor(x[[1]])) {
        levels(x[[1]]) <- c(0, 1)
        x[[1]] <- as.numeric(x[[1]])
        x[[1]] <- x[[1]] - 1
      }

      return(x[[1]])
    }

    subgroup_outcomes <- purrr::map(.x = subgroup_outcomes, check_outcome)

    ## we can't use map here because we need to check and rename the mcmc log file between each iteration
    subgroup_logreg_MCMC_results <- list()
    subgroup_pop_mcmc_tree_results <- list()

    if (find_pop_tree == TRUE) {
      for (i in 1:length(subgroup_outcomes)){
        subgroup_exposure <- subgroup_exposures[[i]]
        subgroup_outcome <- subgroup_outcomes[[i]]

        subgroup_outcome <- as.factor(subgroup_outcome)

        logreg_results <- LogicReg::logreg(
          resp = subgroup_outcome,
          bin = subgroup_exposure,
          select = 7,
          type = 3,
          nleaves = mx_size,
          ntrees = 1,
          mc.control = logreg.mc.control(
            nburn = n_burn_in,
            niter = niter,
            output=4
          )
        )
        #browser()
        subgroup_logreg_MCMC_results[[i]] <- logreg_results
        ## check the mcmc tmp file name, change it for investigation later
        pop_mcmc_tree_results <- get_visiting_mcmc_tree(group_idx = i)
        subgroup_pop_mcmc_tree_results[[i]] <- pop_mcmc_tree_results
        #browser()
      }
    } else {
      subgroup_logreg_MCMC_results <- purrr::map2(
        .x = subgroup_exposures,
        .y = subgroup_outcomes,
        .f = ~ LogicReg::logreg(
          resp = .y,
          bin = .x,
          select = 7,
          type = 3,
          nleaves = mx_size,
          ntrees = 1,
          mc.control = logreg.mc.control(
            nburn = n_burn_in,
            niter = niter,
            hyperpars = hyperparameters
          )
        )
      )
    }

  } else {
    subgroup_logreg_MCMC_results <- NULL
  }

  outcomes <- binary_data %>%
    dplyr::select(
      !!(outcome_label)
    )

  mol_features <- dplyr::select(
    binary_data,
    contains(filter_type)
  )

  if (find_pop_tree == TRUE) {

    main_outcome_logreg_MCMC_results <- LogicReg::logreg(
      resp = outcomes[[1]],
      bin = mol_features,
      select = 7,
      type = 3,
      nleaves = mx_size,
      ntrees = 1,
      mc.control = logreg.mc.control(
        nburn = n_burn_in,
        niter = niter,
        output = 4
    )
  )
     main_pop_mcmc_tree_results<- get_visiting_mcmc_tree(group_idx = "main")

  } else {
    main_outcome_logreg_MCMC_results <- LogicReg::logreg(
      resp = outcomes[[1]],
      bin = mol_features,
      select = 7,
      type = 3,
      nleaves = mx_size,
      ntrees = ntrees,
      mc.control = logreg.mc.control(
        nburn = n_burn_in,
        niter = niter,
        hyperpars = hyperparameters
      )
    )

    main_pop_mcmc_tree_results <- NULL
  }

  return(
    list(
    MCMC_subgroups_analysis = subgroup_logreg_MCMC_results,
    MCMC_main_outcomes_analysis = main_outcome_logreg_MCMC_results,
    subgroup_pop_mcmc_tree_results = subgroup_pop_mcmc_tree_results,
    main_pop_mcmc_tree_results = main_pop_mcmc_tree_results,
    subgroup_outcomes = subgroup_outcomes,
    subgroup_exposures = subgroup_exposures,
    subgroup_names = subgroup_names
  )
  )
}

collect_logic_MCMC_results <- function(MCMC_results,
                                       fingerprints,
                                       niter,
                                       binary_data,
                                       filter_type,
                                       top_fraction) {

  #browser()

  binary_exposure_data <- binary_data %>% dplyr::select(
    contains(filter_type)
  )

  msz_by_iter_MCMC <- MCMC_results$size
  ni_in_model_MCMC <- as.matrix(MCMC_results$single)
  nij_in_model_MCMC <- as.matrix(MCMC_results$double)
  nijk_in_model_MCMC <- as.array(MCMC_results$triple)

  ni_in_model_MCMC_structures <- as.data.frame(cbind(fingerprints, ni_in_model_MCMC))
  colnames(ni_in_model_MCMC_structures) <- c("fingerprints", "ni_in_model_MCMC")
  ni_in_model_MCMC_structures$ni_in_model_MCMC <- as.numeric(ni_in_model_MCMC_structures$ni_in_model_MCMC)

  ni_in_model_MCMC_structures$fraction <- ni_in_model_MCMC_structures$ni_in_model_MCMC / niter
  ni_in_model_MCMC_structures_ordered <- arrange(ni_in_model_MCMC_structures, desc(fraction))

  # ij in same trees

  rownames(nij_in_model_MCMC) <- fingerprints
  colnames(nij_in_model_MCMC) <- fingerprints

  nij_in_model_MCMC_collapsed <- setNames(reshape2::melt(nij_in_model_MCMC), c(
    "Bit Structure 1",
    "Bit Structure 2",
    "count"
  ))

  nij_in_model_MCMC_collapsed$fraction <- nij_in_model_MCMC_collapsed$count / niter
  nij_in_model_MCMC_collapsed_top <- subset(nij_in_model_MCMC_collapsed, fraction > top_fraction)
  # nij_in_model_MCMC_collapsed_top

  chem_obs <- dim(binary_exposure_data)[1]
  pij <- colSums(binary_exposure_data) / chem_obs
  names(pij) <- fingerprints
  pij <- as.data.frame(pij)

  pij_calcs <- merge(nij_in_model_MCMC_collapsed_top,
    pij,
    by.x = "Bit Structure 1",
    by.y = "row.names"
  )

  pij_calcs <- merge(pij_calcs,
    pij,
    by.x = "Bit Structure 2",
    by.y = "row.names"
  )

  names(pij_calcs)[names(pij_calcs) == "pij.y"] <- "pij Bit Structure 1"
  names(pij_calcs)[names(pij_calcs) == "pij.x"] <- "pij Bit Structure 2"

  pij_calcs$expected_fraction <- pij_calcs$`pij Bit Structure 1` * pij_calcs$`pij Bit Structure 2`
  pij_calcs$ratio <- pij_calcs$fraction / pij_calcs$expected_fraction

  # ijk combinations

  dimnames(nijk_in_model_MCMC)[[1]] <- fingerprints
  dimnames(nijk_in_model_MCMC)[[2]] <- fingerprints
  dimnames(nijk_in_model_MCMC)[[3]] <- fingerprints

  nijk_in_model_MCMC_collapsed <- setNames(reshape2::melt(nijk_in_model_MCMC), c(
    "Bit Structure 1",
    "Bit Structure 2",
    "Bit Structure 3",
    "count"
  ))

  nijk_in_model_MCMC_collapsed$fraction <- nijk_in_model_MCMC_collapsed$count / niter
  nijk_in_model_MCMC_collapsed_top <- subset(nijk_in_model_MCMC_collapsed, fraction > top_fraction)
  # nijk_in_model_MCMC_collapsed_top

  pijk_calcs <- merge(nijk_in_model_MCMC_collapsed_top,
    pij,
    by.x = "Bit Structure 1",
    by.y = "row.names"
  )

  names(pijk_calcs)[names(pijk_calcs) == "pij"] <- "pijk Bit Structure 1"

  pijk_calcs <- merge(pijk_calcs,
    pij,
    by.x = "Bit Structure 2",
    by.y = "row.names"
  )

  names(pijk_calcs)[names(pijk_calcs) == "pij"] <- "pijk Bit Structure 2"

  pijk_calcs <- merge(pijk_calcs,
    pij,
    by.x = "Bit Structure 3",
    by.y = "row.names"
  )

  names(pijk_calcs)[names(pijk_calcs) == "pij"] <- "pijk Bit Structure 3"

  pijk_calcs$expected_fraction <- pijk_calcs$`pijk Bit Structure 1` *
    pijk_calcs$`pijk Bit Structure 2` *
    pijk_calcs$`pijk Bit Structure 3`

  pijk_calcs$ratio <- pijk_calcs$fraction / pijk_calcs$expected_fraction

  return(
    list(
    model = MCMC_results,
    model_sizes = msz_by_iter_MCMC,
    vars_in_iterations = ni_in_model_MCMC_structures_ordered,
    ij_interactions = pij_calcs,
    ijk_interactions = pijk_calcs
  )
  )
}
