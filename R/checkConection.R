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
checkConection <- function(data, genotype = "line", trial = "Experiment", response = "YDHA", all = FALSE){

  if(all){
    data %>%
      dplyr::filter(!is.na(.data[[response]])) %>%
      dplyr::select(.data[[genotype]], .data[[trial]]) %>%
      unique.data.frame() %>%
      dplyr::mutate(value = 1) %>%
      tidyr::spread(.data[[trial]], value = value) %>%
      dplyr::mutate(total = rowSums(dplyr::select(., -.data[[genotype]] ), na.rm = T), n = ncol(.)-1, percent = total/n) %>%
      dplyr::arrange(dplyr::desc(total))
  } else{
    data %>%
      dplyr::filter(!is.na(.data[[response]])) %>%
      dplyr::select(.data[[genotype]], .data[[trial]]) %>%
      unique.data.frame() %>%
      dplyr::mutate(value = 1) %>%
      tidyr::spread(.data[[trial]], value = value) %>%
      dplyr::mutate(total = rowSums(dplyr::select(., -.data[[genotype]] ), na.rm = T), n = ncol(.)-1, percent = total/n) %>%
      dplyr::arrange(dplyr::desc(total)) %>%
      dplyr::select(.data[[genotype]], total, n, percent)
  }

}
