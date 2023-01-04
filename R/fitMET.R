#' Stability Coefficients
#'
#' @param predictions A data.frame with one value per GxE combination.
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param trial A character string indicating the column in data that contains
#' trials.
#' @param response A character vector specifying the response variable.
#' @param best A character string specifying how to define the best genotype
#' by numeric value ("min", "max"). "max" by default.
#'
#' @return A data.frame with several stability measures.
#' "superiority" (cultivar-superiority measure), "static" (Shukla's stability
#' variance) and "wricke" (Wricke's ecovalence)
#' @export
#'
#' @examples
#' # in progress
stability <- function(predictions = NULL,
                      genotype = NULL,
                      trial = NULL,
                      response = NULL,
                      best = "max") {
  names_env <- predictions %>%
    type.convert(as.is = FALSE) %>%
    pull(.data[[trial]]) %>%
    levels()
  n_env <- length(names_env)
  E <- mean(predictions$predicted.value, na.rm = TRUE)
  best_var <- ifelse(test = best == "max", yes = "max_trial", no = "min_trial")
  predictions %>%
    group_by(.data[[trial]]) %>%
    mutate(
      max_trial = max(predicted.value, na.rm = TRUE),
      min_trial = min(predicted.value, na.rm = TRUE),
      mean_trial = mean(predicted.value, na.rm = TRUE)
    ) %>%
    group_by(.data[[genotype]]) %>%
    mutate(mean_gen = mean(predicted.value, na.rm = TRUE)) %>%
    summarise(
      superiority = sum(
        (predicted.value - .data[[best_var]])^2
      ) / (2 * n_env) %>% sqrt(),
      static = sum((predicted.value - mean_gen)^2) / (n_env - 1) %>% sqrt(),
      wricke = sum((predicted.value - mean_gen - mean_trial + E)^2) %>% sqrt(),
      predicted.value = mean(predicted.value, na.rm = TRUE)
    ) %>%
    arrange(predicted.value) %>%
    as.data.frame()
}
