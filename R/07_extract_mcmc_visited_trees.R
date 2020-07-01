#' Get MCMC most visiting tree results
#'
#' This funtion reaches into the tmp file created by MCMC logic classification and identifies
#  the tree structure that is most visiting and outputs these results in an interpretable format.
#'
#' @param path Path to the input file
#' @param file_type type of file (xlsx, csv, or txt)
#' @return A matrix
#' @export
#'
get_visiting_mcmc_tree <- function(path = here('R/slogiclisting.tmp'),
                                   group_idx = group) {


  mcmc_trees <- read.delim(path, sep = " ",
                           header = FALSE)

  write.csv(mcmc_trees, here(paste('/R/group_', group_idx,"_mcmc_trees_visited.csv",
                                   sep = "")))

  #browser()

  high_mcmc_tree_visit <- mcmc_trees[which.max(mcmc_trees$V2),]
  prop_log_prob <- high_mcmc_tree_visit[1]
  log_lik <- high_mcmc_tree_visit[2]
  num_visited <- high_mcmc_tree_visit[3]
  tree <- high_mcmc_tree_visit[4:length(high_mcmc_tree_visit)]

  tree_df <- data.frame(matrix(nrow = length(tree),
                               ncol = 6))

  colnames(tree_df) <- c("raw", "number", "conc", "knot", "neg", "pick")
  row.names(tree) <- NULL
  tree_df$raw <- t(as.matrix(tree))
  tree_df$number <- 1:dim(tree_df)[1]

  tree_df$conc <- ifelse(tree_df$raw == 1000, tree_df$conc <- 1,
                         ifelse(tree_df$raw == 2000, tree_df$conc <- 2,
                         tree_df$conc <- 3))

  tree_df$knot <- ifelse(tree_df$conc == 3, tree_df$knot <- abs(tree_df$raw), 0)
  tree_df$neg <- ifelse(tree_df$raw < 0, tree_df$neg <- 1,0)
  tree_df$pick <- ifelse(tree_df$raw != 0, tree_df$pick <- 1,0)
  tree_df <- tree_df[complete.cases(tree_df),]

  tree_df <- tree_df[,-1]

  jpeg(file=here(paste("figures/top_mcmc_tree_group_", group_idx,".jpeg", sep = "")))
  plot_mcmc_logregtree(tree_df)
  dev.off()

  mcmc_formula <- capture.output(print_mcmc_logregtree(tree_df))

  #browser()
  return(list(tree_df = tree_df,
              tree_formula = mcmc_formula,
              prop_log_prob = prop_log_prob,
              log_lik = log_lik,
              num_visited = num_visited))

}

create_popular_tree_var <- function(pop_tree_mcmc_results,
                                    binary_data,
                                    filter_type) {

  pubchem_binary_data <- binary_data %>% select(contains(filter_type))
  #browser()

  pop_tree_formula <- pop_tree_mcmc_results$tree_formula


  tree_formula <- gsub("and", "&", pop_tree_formula)
  tree_formula <- gsub("or", "|", tree_formula)
  tree_formula <- gsub("not", "!", tree_formula)



  for (i in seq(dim(pubchem_binary_data)[2])) {
    target <- paste('X',i, sep ="")
    target_space <- paste(target, "")
    #target <- paste("^", target, sep = "")
    target_end <- paste(target, ')', sep = "")
    target_start <- paste('\\(', target_space,  sep = "")

    pubchem_name <- colnames(pubchem_binary_data)[i]
    pubchem_name_end <-paste(pubchem_name, ')', sep = "")
    pubchem_name_start <-paste('(', pubchem_name, sep = "")


    tree_formula <- gsub(target_start, pubchem_name_start, tree_formula)
    tree_formula <- gsub(target_end, pubchem_name_end, tree_formula)
  }

  pop_tree_var <- binary_data %>%
    mutate(mcmc_pop_tree_var = ifelse(eval(parse(text = tree_formula)), 1, 0)) %>%
    select(mcmc_pop_tree_var)

  #browser()
  resp <- binary_data %>% select(!!outcome_label)

  xtable <- table(pop_tree_var$mcmc_pop_tree_var,
                  resp[[1]])

  cm_results <- caret::confusionMatrix(xtable)



  return(list(pop_tree_var = pop_tree_var,
              perf_table = xtable,
              conf_table_results = cm_results))

}

check_mcmc_accuracy <- function(path,
                                binary_data,
                                resp = resp,
                                filter_type,
                                group_idx) {


  mcmc_trees <- read.delim(path, sep = ",",
                           header = TRUE)

  mcmc_trees <- mcmc_trees[,colSums(is.na(mcmc_trees))<nrow(mcmc_trees)]


  calc_mcmc_tree_accuracy <- function(tree_row,
                                      binary_data,
                                      resp) {

    #browser()

    tree <- tree_row[5:length(tree_row)]
    tree_df <- data.frame(matrix(nrow = length(tree),
                                 ncol = 6))

    colnames(tree_df) <- c("raw", "number", "conc", "knot", "neg", "pick")
    row.names(tree) <- NULL
    tree_df$raw <- as.vector(tree)
    tree_df$number <- 1:dim(tree_df)[1]

    tree_df$conc <- ifelse(tree_df$raw == 1000, tree_df$conc <- 1,
                           ifelse(tree_df$raw == 2000, tree_df$conc <- 2,
                                  tree_df$conc <- 3))

    tree_df$knot <- ifelse(tree_df$conc == 3, tree_df$knot <- abs(tree_df$raw), 0)
    tree_df$neg <- ifelse(tree_df$raw < 0, tree_df$neg <- 1,0)
    tree_df$pick <- ifelse(tree_df$raw != 0, tree_df$pick <- 1,0)
    tree_df <- tree_df[complete.cases(tree_df),]

    tree_df <- tree_df[,-1]

    tree_formula <- capture.output(print_mcmc_logregtree(tree_df, nms = colnames(binary_data)))

    tree_formula <- gsub("and", "&", tree_formula)
    tree_formula <- gsub("or", "|", tree_formula)
    tree_formula <- gsub("not", "!", tree_formula)

    tree_var <- binary_data %>%
      mutate(mol_tree_var = ifelse(eval(parse(text = tree_formula)), 1, 0)) %>%
      select(mol_tree_var)

    xtable <- table(tree_var$mol_tree_var, resp)

    perf_results <- caret::confusionMatrix(xtable)

    balanced_accuracy <- perf_results$byClass[[11]]

    return(balanced_accuracy)

  }

  accuracies <- apply(mcmc_trees, 1, calc_mcmc_tree_accuracy,
                      binary_data = binary_data,
                      resp = resp)

  mcmc_trees$accuracy <- accuracies

  return(mcmc_trees)
}



















