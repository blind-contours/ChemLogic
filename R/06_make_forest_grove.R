#' Take rules for each subgroup found and create one binary rule that represents all boolean logic trees
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
make_groves <- function(subgroup_tree_formula_results,
                        main_tree_formula_results,
                        input_data,
                        filter_type,
                        outcome_label) {
  formulas <- c()
  for (i in 1:length(subgroup_tree_formula_results)) {
    formulas[i] <- subgroup_tree_formula_results[[i]][[1]]$tree_formula
  }

  input_data <- input_data %>%
    mutate(main_outcome_tree = ifelse(eval(parse(text = main_tree_formula_results$MCMC_main_outcomes_analysis[[1]]$tree_formula[1])), 1, 0))

  tree_variables <- c()
  for (i in 1:length(formulas)) {
    tree_variable <- paste("mol_tree_", i, sep = "")

    input_data <- input_data %>%
      mutate(!!as.name(tree_variable) := ifelse(eval(parse(text = formulas[i])), 1, 0))
    tree_variables[i] <- tree_variable
  }

  trees_concat <- paste(tree_variables, collapse = "== 1 | ")
  trees_concat <- paste(trees_concat, "== 1 ")

  input_data <- input_data %>%
    mutate(subgroup_trees = ifelse(eval(parse(text = trees_concat)), 1, 0))

  input_data <- input_data %>%
    mutate(subgroup_tree_w_main = ifelse(subgroup_trees == 1 | main_outcome_tree == 1, 1, 0))

  input_data <- input_data %>%
    mutate(false_negatives = ifelse(subgroup_trees == 0 & Label_names == 1, 1, 0))

  false_neg_controls_subset <- filter(input_data, !!outcome_label == 0) %>%
    dplyr::select(contains(filter_type), !!outcome_label)

  false_neg_features_subset <- filter(input_data, false_negatives == 1) %>%
    dplyr::select(contains(filter_type), !!outcome_label)

  false_neg_data <- rbind(false_neg_controls_subset, false_neg_features_subset)

  false_neg_outcomes <- false_neg_data %>%
    dplyr::select(!!outcome_label)

  false_neg_features <- false_neg_data %>%
    dplyr::select(contains(filter_type))

  MCMC_vars_select_tree_construction <- logreg(
    resp = as.factor(false_neg_outcomes[[1]]),
    bin = false_neg_features,
    type = 1,
    select = 2,
    ntrees = 1,
    nleaves = c(5, 10)
  )

  fp_best_tree_results <- which.min(MCMC_vars_select_tree_construction$allscores[, 1])
  fp_tree_formula <- MCMC_vars_select_tree_construction$alltrees[fp_best_tree_results]

  tree_formula <- capture.output(fp_tree_formula[[1]])
  tree_formula <- gsub(" \\+1 \\* ", "", tree_formula)

  tree_formula <- gsub("and", "&", tree_formula)
  tree_formula <- gsub("or", "|", tree_formula)
  tree_formula <- gsub("not", "!", tree_formula)

  for (i in seq(dim(false_neg_features)[2])) {
    target <- paste("X", i, sep = "")
    target_space <- paste(target, "")
    # target <- paste("^", target, sep = "")
    target_end <- paste(target, ")", sep = "")
    target_start <- paste("\\(", target_space, sep = "")

    pubchem_name <- colnames(false_neg_features)[i]
    pubchem_name_end <- paste(pubchem_name, ")", sep = "")
    pubchem_name_start <- paste("(", pubchem_name, sep = "")


    tree_formula <- gsub(target_start, pubchem_name_start, tree_formula)
    tree_formula <- gsub(target_end, pubchem_name_end, tree_formula)
  }

  #browser()

  input_data <- input_data %>%
    mutate(fp_formula = ifelse(eval(parse(text = tree_formula[1])), 1, 0))

  input_data <- input_data %>%
    mutate(subgroup_trees_w_fp = ifelse(subgroup_trees == 1 | fp_formula == 1, 1, 0))

  input_data <- input_data %>%
    mutate(subgroup_trees_w_fp_main = ifelse(subgroup_trees == 1 | fp_formula == 1 | main_outcome_tree == 1, 1, 0))

  main_tree_results <- caret::confusionMatrix(table(input_data$Label_names, input_data$main_outcome_tree))
  subgroup_tree_results <- caret::confusionMatrix(table(input_data$Label_names, input_data$subgroup_trees))
  subgroup_main_results <- caret::confusionMatrix(table(input_data$Label_names, input_data$subgroup_tree_w_main))
  #subgroup_fp_results <- caret::confusionMatrix(table(input_data$Label_names, input_data$subgroup_trees_w_fp))
  #subgroup_all_results <- caret::confusionMatrix(table(input_data$Label_names, input_data$subgroup_trees_w_fp_main))

  return(list(main = main_tree_results,
              subgroup = subgroup_tree_results,
              subgroup_main = subgroup_main_results,
              updated_data = input_data))
}
