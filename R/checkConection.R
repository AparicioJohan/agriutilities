#' checkconection
#'
#' @param data MET dataset
#' @param genotype string
#' @param trial string
#' @param response string
#' @param all wheater or not print all the table
#'
#' @return table with gen conection
#' @export
#'
#' @examples
#' # library(agridat)
#' # data(besag.met)
#' # dat <- besag.met
#' # checkConection(data = dat, genotype = "gen", trial = "county", response = "yield", all = T)
#' @importFrom rlang .data
#' @import dplyr
checkConection <- function(data,
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
      total = rowSums(dplyr::select(., -.data[[genotype]]), na.rm = T),
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


#' checkConection2
#' @description This function generates a (n_trial x n_trial) matrix  with the
#'  amount of genotypes shared between each pair of trial.
#'
#' @param data MET dataset
#' @param genotype string
#' @param trial string
#' @param response string
#'
#' @return matrix
#' @export
#'
#' @examples
#' # library(agridat)
#' # data(besag.met)
#' # dat <- besag.met
#' # conn <- checkConection2(data = dat, genotype = "gen", trial = "county", response = "yield")
#' # heatmap(conn)
#' # res <- factoextra::hcut(conn, k = 3, stand = FALSE)
#' # factoextra::fviz_dend(
#' #   x = res,
#' #   rect = FALSE,
#' #   cex = 0.5,
#' #   lwd = 0.5,
#' #   main = "Dendrogram",
#' #   horiz = TRUE
#' # )
#' @importFrom rlang .data
#' @import dplyr tidyr tibble
checkConection2 <- function(data,
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
