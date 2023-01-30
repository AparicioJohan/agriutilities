#' Check connectivity between trials
#'
#' @param data A data.frame in a wide format.
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param trial A character string indicating the column in data that contains
#' trials.
#' @param response A character string specifying the trait.
#' @param all Whether or not print all the table.
#'
#' @return Table with the genotype connectivity.
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
#' @import dplyr
check_connectivity <- function(data = NULL,
                               genotype = "line",
                               trial = "Experiment",
                               response = NULL,
                               all = FALSE) {
  tmp_data <- data %>%
    {
      if (!is.null(response)) {
        dplyr::filter(.data = ., !is.na(.data[[response]]))
      } else {
        .
      }
    } %>%
    dplyr::select(.data[[genotype]], .data[[trial]]) %>%
    unique.data.frame() %>%
    dplyr::mutate(value = 1) %>%
    tidyr::spread(.data[[trial]], value = value) %>%
    dplyr::mutate(
      total = rowSums(dplyr::select(., -.data[[genotype]]), na.rm = TRUE),
      n = ncol(.) - 1,
      percent = total / n
    ) %>%
    dplyr::arrange(dplyr::desc(total))

  if (all) {
    tmp_data
  } else {
    tmp_data %>%
      dplyr::select(.data[[genotype]], total, n, percent)
  }
}

#' Connectivity Matrix
#' @description Check the amount of genotypes shared between each pair of trial.
#'
#' @param data A data.frame in a wide format.
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param trial A character string indicating the column in data that contains
#' trials.
#' @param response A character string specifying the trait.
#'
#' @return This function generates a (n_trial x n_trial) matrix with the amount
#' of genotypes shared between each pair of trial.
#' @export
#'
#' @examples
#' \donttest{
#' library(agridat)
#' data(besag.met)
#' dat <- besag.met
#' connectivity_matrix(
#'   data = dat,
#'   genotype = "gen",
#'   trial = "county",
#'   response = "yield"
#' )
#' }
#' @importFrom rlang .data
#' @import dplyr tidyr tibble
connectivity_matrix <- function(data = NULL,
                                genotype = "germplasmName",
                                trial = "trial",
                                response = NULL) {
  tmp_data <- data %>%
    {
      if (!is.null(response)) {
        dplyr::filter(.data = ., !is.na(.data[[response]]))
      } else {
        .
      }
    } %>%
    dplyr::select(.data[[genotype]], .data[[trial]]) %>%
    unique.data.frame() %>%
    dplyr::mutate(value = 1) %>%
    tidyr::spread(.data[[trial]], value = value) %>%
    tibble::column_to_rownames(genotype) %>%
    as.matrix()
  tmp_data[is.na(tmp_data)] <- 0
  conectivity <- t(tmp_data) %*% tmp_data
  return(conectivity)
}
