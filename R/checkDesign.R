#' Check Experimental Design
#'
#' @param data A data.frame in a wide format.
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param trial A character string indicating the column in data that contains
#' trials.
#' @param traits A character vector specifying the traits for which the models
#' should be fitted.
#' @param rep A character string indicating the column in data that contains
#' replicates.
#' @param block A character string indicating the column in data that contains
#' sub blocks.
#' @param row A character string indicating the column in data that contains the
#' row coordinates
#' @param col A character string indicating the column in data that contains the
#' column coordinates.
#'
#' @return An object of class \code{checkAgri}, with a list of:
#' \item{summ_traits}{A data.frame containing a summary of the traits.}
#' \item{exp_design_resum}{A data.frame containing a summary of the experimental
#'  design.}
#' \item{filter}{A list by trait containing the filtered trials.}
#' \item{exp_design_list}{A data.frame containing the experimental design of
#' each trial.}
#' \item{check_connectivity}{A data.frame with the genotype connectivity.}
#' \item{connectivity_matrix}{A matrix with the amount of genotypes shared
#' between each pair of trial.}
#' \item{data_design}{A data frame containing the data used with an additional
#' column of the experimental design.}
#' \item{inputs}{A list containing the character string that indicates the
#' column in data that contains the genotype, trial, traits, rep, block, row
#' and col.}
#' @export
#'
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_met(
#'   data = dat,
#'   genotype = "gen",
#'   trial = "county",
#'   traits = c("yield"),
#'   rep = "rep",
#'   block = "block",
#'   col = "col",
#'   row = "row"
#' )
#' print(results)
#' # Plotting Connectivity Matrix
#' plot(results, type = "connectivity")
#' # Plotting Missing Data
#' plot(results, type = "missing")
#' }
#' @import dplyr
#' @importFrom stats sd median
#' @importFrom utils type.convert
check_design_met <- function(data = NULL,
                             genotype = NULL,
                             trial = NULL,
                             traits = NULL,
                             rep = NULL,
                             block = NULL,
                             row = NULL,
                             col = NULL) {
  if (is.null(data)) {
    stop("Error: data not found")
  }
  if (is.null(genotype)) {
    stop("No 'genotype' name column provided")
  }
  if (!genotype %in% names(data)) {
    stop(paste("No '", genotype, "' found in the data"))
  }
  if (is.null(trial)) {
    message("No 'trial' name column provided")
    trial <- "trial"
    data <- data %>% mutate(trial = NA)
  }
  if (!trial %in% names(data)) {
    stop(paste("No '", trial, "' found in the data"))
  }
  if (is.null(traits)) {
    stop("No 'traits' argument provided")
  }
  if (is.null(rep)) {
    message("No 'rep' name column provided")
    rep <- "rep"
    data <- data %>% mutate(rep = NA)
  }
  if (!rep %in% names(data)) {
    stop(paste("No '", rep, "' found in the data"))
  }
  if (is.null(block)) {
    message("No 'block' name column provided")
    block <- "block"
    data <- data %>% mutate(block = NA)
  }
  if (!block %in% names(data)) {
    stop(paste("No '", block, "' found in the data"))
  }
  if (is.null(row)) {
    message("No 'row' name column provided")
    row <- "row"
    data <- data %>% mutate(row = NA)
  }
  if (!row %in% names(data)) {
    stop(paste("No '", row, "' found in the data"))
  }
  if (is.null(col)) {
    message("No 'col' name column provided")
    col <- "col"
    data <- data %>% mutate(col = NA)
  }
  if (!col %in% names(data)) {
    stop(paste("No '", col, "' found in the data"))
  }
  for (i in traits) {
    if (!i %in% names(data)) {
      stop(paste("No '", i, "' column found"))
    }
    class_trait <- data[[i]] %>% class()
    if (!class_trait %in% c("numeric", "integer")) {
      stop(
        paste0("The class of the trait '", i, "' should be numeric or integer.")
      )
    }
  }
  trait_issue <- traits[make.names(traits) != traits]
  traits <- traits[make.names(traits) == traits]
  if (length(trait_issue) != 0) {
    message(
      "Issues in variables names. The trait '",
      paste(trait_issue, collapse = " - "),
      "' has been removed."
    )
    if (length(traits) == 0) {
      stop("Check the variable names.")
    }
  }

  summ_traits <- data %>%
    dplyr::select(.data[[trial]], all_of(traits)) %>%
    gather(data = ., key = "traits", value = "value", -.data[[trial]]) %>%
    group_by(.data[[trial]], traits) %>%
    summarise(
      Mean = mean(value, na.rm = TRUE),
      Median = median(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      CV = SD / Mean,
      n = n(),
      n_miss = sum(is.na(value)),
      miss_perc = n_miss / n,
      .groups = "drop"
    )

  exp_design_resum <- data %>%
    group_by(.data[[trial]], .data[[genotype]]) %>%
    mutate(gen_reps = n()) %>%
    group_by(.data[[trial]], .data[[rep]]) %>%
    mutate(size_reps = n()) %>%
    group_by(.data[[trial]]) %>%
    summarise(
      n = n(),
      n_gen = n_distinct(.data[[genotype]], na.rm = TRUE),
      n_rep = n_distinct(.data[[rep]], na.rm = TRUE),
      n_block = n_distinct(.data[[block]], na.rm = TRUE),
      n_col = n_distinct(.data[[col]], na.rm = TRUE),
      n_row = n_distinct(.data[[row]], na.rm = TRUE),
      num_of_reps = paste(sort(unique(gen_reps)), collapse = "_"),
      num_of_gen = paste(
        table(gen_reps) / sort(unique(gen_reps)),
        collapse = "_"
      ),
      reps_size = paste(unique(size_reps), collapse = "_"),
      reps_equal = length(unique(size_reps)) == 1,
      unrep = ifelse(n_gen == n, TRUE, FALSE),
      rcbd = ifelse(reps_equal && n_rep > 1, "rcbd", NA),
      alpha_lattice = ifelse(rcbd == "rcbd" & n_block > 1, "alpha", NA),
      prep = ifelse(nchar(num_of_reps) > 1, "prep", NA),
      spatial = ifelse(
        test = n_col > 1 & n_row > 1 & !rcbd %in% "rcbd",
        yes = "row_col",
        no = ifelse(
          test = n_col > 1 & n_row > 1 & rcbd %in% "rcbd",
          yes = "res_row_col",
          no = NA
        )
      )
    ) %>%
    type.convert(as.is = FALSE)

  trials_spatial <- exp_design_resum %>%
    filter(!is.na(spatial)) %>%
    pull(trial) %>%
    as.character()

  if (!is.null(trials_spatial)) {
    duplicates <- data %>%
      filter(.data[[trial]] %in% trials_spatial) %>%
      select(.data[[trial]], .data[[col]], .data[[row]]) %>%
      na.omit() %>%
      split(x = ., .[, trial]) %>%
      lapply(FUN = function(x) ifelse(sum(duplicated(x)) >= 1, TRUE, FALSE)) %>%
      unlist()
    duplicates <- names(duplicates[which(duplicates %in% TRUE)])
  } else {
    duplicates <- ""
  }

  filters_summary <- list()

  for (i in traits) {

    # Filter missing trials
    trials_to_remove_miss <- summ_traits %>%
      filter(traits %in% i) %>%
      filter(miss_perc >= 0.50) %>%
      pull(.data[[trial]]) %>%
      as.character()

    # Filter by non_var_gen
    non_var_gen <- data %>%
      filter(!.data[[trial]] %in% trials_to_remove_miss) %>%
      group_by(.data[[trial]], .data[[genotype]]) %>%
      summarise(
        mean = mean(.data[[i]], na.rm = TRUE),
        sd = sd(.data[[i]], na.rm = TRUE),
        cv = sd / mean,
        .groups = "drop"
      ) %>%
      group_by(.data[[trial]]) %>%
      summarise(var_total = sum(sd, na.rm = TRUE), .groups = "drop") %>%
      arrange(var_total)

    # Filter variation
    trials_to_remove_novar <- non_var_gen %>%
      filter(var_total == 0) %>%
      pull(trial) %>%
      as.character()

    filters_summary[[i]] <- list(
      "missing_50%" = trials_to_remove_miss,
      "no_variation" = trials_to_remove_novar,
      "row_col_dup" = duplicates,
      "trials_to_remove" = union(
        x = union(trials_to_remove_miss, trials_to_remove_novar),
        y = duplicates
      )
    )
  }

  design <- exp_design_resum %>%
    mutate(
      exp_design = ifelse(
        test = spatial %in% "row_col",
        yes = "row_col",
        no = ifelse(
          test = spatial %in% "res_row_col",
          yes = "res_row_col",
          no = ifelse(
            test = alpha_lattice %in% "alpha",
            yes = "alpha_lattice",
            no = ifelse(
              test = rcbd %in% "rcbd",
              yes = "rcbd",
              no = ifelse(
                test = prep %in% "prep",
                yes = "crd",
                no = "crd"
              )
            )
          )
        )
      )
    ) %>%
    select(.data[[trial]], exp_design) %>%
    as.data.frame()

  data_design <- merge(x = data, y = design, all.x = TRUE) %>%
    mutate(id = seq_len(nrow(.))) %>%
    relocate(id)

  # Conectivity
  conn1 <- check_connectivity(
    data = data,
    genotype = genotype,
    trial = trial,
    response = NULL
  )

  conn2 <- check_connectivity(
    data = data,
    genotype = genotype,
    trial = trial,
    response = NULL,
    return_matrix = TRUE
  )

  objt <- list(
    summ_traits = summ_traits,
    exp_design_resum = exp_design_resum,
    filter = filters_summary,
    exp_design_list = design,
    check_connectivity = conn1,
    connectivity_matrix = conn2,
    data_design = data_design,
    inputs = list(
      genotype = genotype,
      trial = trial,
      traits = traits,
      rep = rep,
      block = block,
      row = row,
      col = col
    )
  )
  class(objt) <- "checkAgri"
  return(objt)
}
