#' Check connectivity between trials
#'
#' @param data A data.frame in a wide format.
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param trial A character string indicating the column in data that contains
#' trials.
#' @param response A character string specifying the trait.
#' @param all Whether or not print all the table.
#' @param return_matrix A logical value indicating if the user wants to return
#' a (n_trial x n_trial) matrix with the amount of genotypes shared between each
#' pair of trial. (\code{FALSE} by default)
#'
#' @return A data.frame with the genotype connectivity. If return_matrix is TRUE,
#' it will return a n_trial x n_trial matrix with the amount of genotypes shared
#' between each pair of trial.
#' @export
#'
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' head(
#'   check_connectivity(
#'     data = dat,
#'     genotype = "gen",
#'     trial = "county",
#'     response = "yield",
#'     all = TRUE
#'   )
#' )
#' }
#' @importFrom rlang .data
#' @import dplyr tidyr tibble
check_connectivity <- function(data = NULL,
                               genotype = "line",
                               trial = "Experiment",
                               response = NULL,
                               all = FALSE,
                               return_matrix = FALSE) {
  tmp_data <- data %>%
    {
      if (!is.null(response)) {
        filter(.data = ., !is.na(.data[[response]]))
      } else {
        .
      }
    } %>%
    select(all_of(c(genotype, trial))) %>%
    unique.data.frame() %>%
    mutate(value = 1) %>%
    tidyr::spread(all_of(trial), value = value)

  if (return_matrix) {
    conn <- tmp_data %>%
      tibble::column_to_rownames(genotype) %>%
      as.matrix()
    conn[is.na(conn)] <- 0
    conectivity <- t(conn) %*% conn
    return(conectivity)
  } else {
    connection_table <- tmp_data %>%
      mutate(
        total = rowSums(
          x = select(., -all_of(genotype)),
          na.rm = TRUE
        ),
        n = ncol(.) - 1,
        percent = total / n
      ) %>%
      arrange(desc(total))
    if (all) {
      return(connection_table)
    } else {
      connection_table <- connection_table %>%
        select(all_of(genotype), total, n, percent)
      return(connection_table)
    }
  }
}
