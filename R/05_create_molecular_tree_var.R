#' Load Chemical Fingerprint Data
#'
#' This function loads a file as a matrix. All explanatory data should be binary,
#' indicating presence or absence of a bit structure, the first column
#' should have the chemical name in order to join on additional physical data late.
#'
#' @param path Path to the input file
#' @param file_type type of file (xlsx, csv, or txt)
#' @return A matrix
#' @export
#'
run_tree_construction <- function(resp, bin, mol_tree_vars) {

  MCMC_vars_select_tree_construction <- logreg(
    resp = resp,
    bin = bin,
    type = 1,
    select = 1,
    ntrees = 1,
    kfold = 10,
    nleaves = dim(mol_tree_vars)[2]
  )

  tree_formula <- MCMC_vars_select_tree_construction$model$trees[[1]]
  tree_formula <- capture.output(tree_formula)
  tree_formula <- gsub("and", "&", tree_formula)
  tree_formula <- gsub("or", "|", tree_formula)
  tree_formula <- gsub("not", "!", tree_formula)
  tree_formula <- gsub("X", "PubchemFP", tree_formula)

  mol_tree_interaction_var <- binary_data %>%
    mutate(mol_tree_var = ifelse(eval(parse(text = tree_formula)), 1, 0)) %>%
    select(mol_tree_var)

  outcome <- binary_data %>% select(!!outcome_label)
  outcome <- as.vector(outcome[[1]])
  xtable <- table(mol_tree_interaction_var$mol_tree_var, outcome)

  if (dim(xtable)[1] > 1) {
    accuracy <- sum(diag(xtable))/sum(xtable)
  } else {accuracy <- NULL}

  return(list(
    tree_formula = tree_formula,
    tree_variable = mol_tree_interaction_var,
    accuracy = accuracy
  ))

}


create_molecular_tree_var <- function(MCMC_results,
                                      couple_thresh = 0.7,
                                      triple_thresh = 0.5,
                                      top_MCMC_vars_thresh = 0.8,
                                      binary_data,
                                      make_logic_tree = TRUE,
                                      expand_boolean_grid = FALSE,
                                      filter_type = "Pubchem",
                                      fingerprints,
                                      fingerprints_idx,
                                      outcome_label,
                                      iters) {

  top_individual_vars <- subset(MCMC_results$vars_in_iterations, fraction >= top_MCMC_vars_thresh)
  ij_interactions_top <- subset(MCMC_results$ij_interactions, fraction >= couple_thresh)
  ijk_interactions_top <- subset(MCMC_results$ijk_interactions, fraction >= triple_thresh)

  vars_ind_top <- as.vector(top_individual_vars$Bit_Substructure)
  ij_vars_top <- c(as.vector(ij_interactions_top$`Bit Structure 1`), as.vector(ij_interactions_top$`Bit Structure 2`))
  ijk_vars_top <- c(as.vector(ijk_interactions_top$`Bit Structure 1`), as.vector(ijk_interactions_top$`Bit Structure 2`), as.vector(ijk_interactions_top$`Bit Structure 3`))

  vars_to_use <- unique(c(vars_ind_top, ij_vars_top, ijk_vars_top))

  binary_data_exposures <- binary_data %>% select(contains(filter_type))
  colnames(binary_data_exposures) <- fingerprints[, fingerprints_idx]

  mol_tree_vars <- binary_data_exposures %>% select(vars_to_use)
  mol_tree_outcome <- binary_data %>% select(!!(outcome_label))

  if (make_logic_tree) {

    MC_results <- map(seq_len(iters), ~run_tree_construction(resp = mol_tree_outcome[, 1],
                                                             bin = mol_tree_vars,
                                                             mol_tree_vars = mol_tree_vars))

    max_accuracy_idx <- which.max(unlist(sapply(MC_results, `[[`, "accuracy")))
    best_tree_results <- MC_results[max_accuracy_idx]


  } else {

    ## run expand_boolean_grid if size isn't too large

  }

  return(best_tree_results)

}
