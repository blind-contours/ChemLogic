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

run_tree_construction <- function(resp,
                                  mol_tree_vars,
                                  binary_data) {

  #browser()
  if (dim(mol_tree_vars)[2] == dim(binary_data)[2]) {

    MCMC_vars_select_tree_construction <- logreg(
      resp = resp,
      bin = mol_tree_vars,
      type = 1,
      select = 2,
      ntrees = 1,
      kfold = 5,
      nleaves = c(1,10)
    )


  } else {

    MCMC_vars_select_tree_construction <- logreg(
      resp = resp,
      bin = mol_tree_vars,
      type = 1,
      select = 2,
      ntrees = 1,
      kfold = 5,
      nleaves = c(1,dim(mol_tree_vars)[2])
    )
  }


  #browser()

  #names(mol_tree_vars) <- paste0('X', 1:(ncol(mol_tree_vars)))

  min_model_indx <- which.min(MCMC_vars_select_tree_construction$allscores[,1])
  tree_formula <- MCMC_vars_select_tree_construction$alltrees[min_model_indx][[1]]

  tree_formula <- capture.output(tree_formula)
  tree_formula <- gsub(" \\+1 \\* ", "", tree_formula)

  tree_formula <- gsub("and", "&", tree_formula)
  tree_formula <- gsub("or", "|", tree_formula)
  tree_formula <- gsub("not", "!", tree_formula)
  #tree_formula <- gsub("X", "PubchemFP", tree_formula)

  for (i in seq(dim(mol_tree_vars)[2])) {
    target <- paste('X',i, sep ="")
    target_space <- paste(target, "")
    #target <- paste("^", target, sep = "")
    target_end <- paste(target, ')', sep = "")
    target_start <- paste('\\(', target_space,  sep = "")

    pubchem_name <- colnames(mol_tree_vars)[i]
    pubchem_name_end <-paste(pubchem_name, ')', sep = "")
    pubchem_name_start <-paste('(', pubchem_name, sep = "")


    tree_formula <- gsub(target_start, pubchem_name_start, tree_formula)
    tree_formula <- gsub(target_end, pubchem_name_end, tree_formula)

    #print(i)
    #print(tree_formula)
  }

  binary_data_w_mol_tree_var <- mol_tree_vars %>%
    mutate(mol_tree_var = ifelse(eval(parse(text = tree_formula)), 1, 0)) %>%
    select(mol_tree_var)


  #outcome <- binary_data %>% select(!!outcome_label)
  #outcome <- as.vector(outcome[[1]])
  xtable <- table(binary_data_w_mol_tree_var$mol_tree_var, resp)

  if (dim(xtable)[1] > 1) {
    xtable_cm <- caret::confusionMatrix(xtable)
    sensitivity <- xtable_cm$byClass[1]
    specificity <- xtable_cm$byClass[2]
    balanced_accuracy <- xtable_cm$byClass[11]
  } else {
    balanced_accuracy <- 0
    sensitivity <- 0
    specificity <-0
  }

  return(list(
    tree_formula = tree_formula,
    data_w_tree_variable = binary_data_w_mol_tree_var,
    balanced_accuracy = balanced_accuracy[[1]],
    sensitivity = sensitivity[[1]],
    specificity = specificity[[1]],
    table = xtable
  ))

}



run_mcmc_tree_construction <- function(resp,
                                  mol_tree_vars,
                                  binary_data,
                                  group,
                                  niter,
                                  n_burn_in) {

  resp <- factor(resp, levels = c(0,1))

  # browser()
  if (dim(mol_tree_vars)[2] == dim(binary_data)[2]) {

    MCMC_vars_select_tree_construction <- LogicReg::logreg(
      resp = resp,
      bin = mol_tree_vars,
      select = 7,
      type = 3,
      nleaves = 20,
      ntrees = 1,
      mc.control = logreg.mc.control(
        nburn = n_burn_in,
        niter = niter,
        output = 4
      )
    )
  } else {
    #browser()
    MCMC_vars_select_tree_construction <- LogicReg::logreg(
      resp = resp,
      bin = mol_tree_vars,
      select = 7,
      type = 3,
      nleaves = dim(mol_tree_vars)[2],
      ntrees = 1,
      mc.control = logreg.mc.control(
        nburn = n_burn_in,
        niter = niter,
        output = 4
      )
    )
  }

  mcmc_results <- get_visiting_mcmc_tree(group = group)
  #browser()
  mcmc_accuracies <- check_mcmc_accuracy(path = here(paste('/R/group_', group,"_mcmc_trees_visited.csv", sep = "")),
                      binary_data = mol_tree_vars,
                      resp = resp,
                      group_idx = group)


  # browser()
  # names(mol_tree_vars) <- paste0('X', 1:(ncol(mol_tree_vars)))

  # min_model_indx <- which.min(MCMC_vars_select_tree_construction$allscores[,1])
  # tree_formula <- MCMC_vars_select_tree_construction$alltrees[min_model_indx][[1]]

  # tree_formula <- capture.output(tree_formula)
  # tree_formula <- gsub(" \\+1 \\* ", "", tree_formula)

  # browser()
  tree_formula <- mcmc_results$tree_formula
  tree_formula <- gsub("and", "&", tree_formula)
  tree_formula <- gsub("or", "|", tree_formula)
  tree_formula <- gsub("not", "!", tree_formula)
  # tree_formula <- gsub("X", "PubchemFP", tree_formula)

  for (i in seq(dim(mol_tree_vars)[2])) {
    target <- paste("X", i, sep = "")
    target_space <- paste(target, "")
    # target <- paste("^", target, sep = "")
    target_end <- paste(target, ")", sep = "")
    target_start <- paste("\\(", target_space, sep = "")

    pubchem_name <- colnames(mol_tree_vars)[i]
    pubchem_name_end <- paste(pubchem_name, ")", sep = "")
    pubchem_name_start <- paste("(", pubchem_name, sep = "")


    tree_formula <- gsub(target_start, pubchem_name_start, tree_formula)
    tree_formula <- gsub(target_end, pubchem_name_end, tree_formula)

    # print(i)
    # print(tree_formula)
  }

  binary_data_w_mol_tree_var <- mol_tree_vars %>%
    mutate(mol_tree_var = ifelse(eval(parse(text = tree_formula)), 1, 0)) %>%
    select(mol_tree_var)


  # outcome <- binary_data %>% select(!!outcome_label)
  # outcome <- as.vector(outcome[[1]])
  xtable <- table(binary_data_w_mol_tree_var$mol_tree_var, resp)

  if (dim(xtable)[1] > 1) {
    xtable_cm <- caret::confusionMatrix(xtable)
    sensitivity <- xtable_cm$byClass[1]
    specificity <- xtable_cm$byClass[2]
    balanced_accuracy <- xtable_cm$byClass[11]
  } else {
    balanced_accuracy <- 0
    sensitivity <- 0
    specificity <- 0
  }

  return(list(
    tree_formula = tree_formula,
    data_w_tree_variable = binary_data_w_mol_tree_var,
    balanced_accuracy = balanced_accuracy[[1]],
    sensitivity = sensitivity[[1]],
    specificity = specificity[[1]],
    table = xtable,
    mcmc_trees_w_accuracies = mcmc_accuracies
  )
  )
}

create_molecular_tree_var <- function(.x,
                                      .y,
                                      .z,
                                      .a,
                                      group,
                                      couple_thresh = 0.7,
                                      triple_thresh = 0.5,
                                      top_MCMC_vars_thresh = 0.8,
                                      binary_data,
                                      make_logic_tree = TRUE,
                                      expand_boolean_grid = FALSE,
                                      filter_type = "Pubchem",
                                      fingerprints,
                                      outcome_label,
                                      iters,
                                      n_burn_in,
                                      subgroup = TRUE,
                                      mcmc = FALSE) {
  #browser()
  MCMC_results <- .x
  top_individual_vars <- subset(MCMC_results$vars_in_iterations, fraction >= top_MCMC_vars_thresh)
  ij_interactions_top <- subset(MCMC_results$ij_interactions, fraction >= couple_thresh)
  ijk_interactions_top <- subset(MCMC_results$ijk_interactions, fraction >= triple_thresh)

  vars_ind_top <- as.vector(top_individual_vars$fingerprints)
  ij_vars_top <- c(as.vector(ij_interactions_top$`Bit Structure 1`), as.vector(ij_interactions_top$`Bit Structure 2`))
  ijk_vars_top <- c(
    as.vector(ijk_interactions_top$`Bit Structure 1`), as.vector(ijk_interactions_top$`Bit Structure 2`),
    as.vector(ijk_interactions_top$`Bit Structure 3`)
  )

  # browser()

  vars_to_use <- unique(c(vars_ind_top, ij_vars_top, ijk_vars_top))

  if (subgroup) {
    binary_data_exposures <- .z
    subgroup_names <- .a

    # browser()
    binary_data <- filter(binary_data, Name %in% subgroup_names$Name)
    binary_data <- binary_data %>% select(contains(filter_type))
    colnames(binary_data) <- fingerprints
    pubchem_match_idxs <- match(vars_to_use, colnames(binary_data))

    pubchem_names <- colnames(binary_data_exposures)[pubchem_match_idxs]

    mol_tree_vars <- binary_data_exposures %>% select(pubchem_names)

    # browser()

    if (length(mol_tree_vars) <= 5) {
      mol_tree_vars <- binary_data_exposures
    } else {
      mol_tree_vars <- mol_tree_vars
    }

    mol_tree_outcome <- .y
  } else {
    binary_data_exposures <- binary_data %>% select(contains(filter_type))
    binary_data_exposures_unlabeled <- binary_data_exposures

    colnames(binary_data_exposures) <- fingerprints
    pubchem_match_idxs <- match(vars_to_use, colnames(binary_data_exposures))
    pubchem_names <- colnames(binary_data_exposures_unlabeled)[pubchem_match_idxs]

    mol_tree_vars <- binary_data_exposures_unlabeled %>% select(pubchem_names)

    mol_tree_outcome <- binary_data %>% select(!!(outcome_label))

    mol_tree_outcome <- mol_tree_outcome[[1]]

    binary_data_exposures <- binary_data_exposures_unlabeled
  }

  mol_tree_outcome <- as.vector(mol_tree_outcome)

  if (make_logic_tree) {

    if (mcmc) {
    MC_results <- run_mcmc_tree_construction(
      resp = mol_tree_outcome,
      mol_tree_vars = mol_tree_vars,
      binary_data = binary_data_exposures,
      group,
      n_burn_in = n_burn_in,
      niter = iters
    )

    }

    MC_results <-  purrr::map(seq_len(5), ~run_tree_construction(resp = mol_tree_outcome,
                                                                     mol_tree_vars = mol_tree_vars,
                                                                     binary_data = binary_data_exposures))

    max_accuracy_idx <- which.max(sapply(MC_results, `[[`, "balanced_accuracy"))
    MC_results <- MC_results[max_accuracy_idx]

  } else {

    ## run expand_boolean_grid if size isn't too large
  }

  return(MC_results)
}
