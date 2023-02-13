#' Genetic Gain Parameters
#'
#' @param model Linear regression model (\code{lm} object)
#' @param trait A character string indicating the column in data that contains
#' trials.
#'
#' @return A data.frame with some parameters from the linear regression (Slope,
#' se_Slope, Intercept, r2, Pr(>F)) and the percentage of Genetic Gain.
#' @export
#'
#' @examples
#' library(ggplot2)
#' library(agridat)
#' library(magrittr)
#' library(agriutilities)
#'
#' data(baker.barley.uniformity)
#' dat <- baker.barley.uniformity
#' head(dat)
#'
#' model <- lm(yield ~ year, dat)
#'
#' dat %>%
#'   na.omit() %>%
#'   ggplot(
#'     aes(x = year, y = yield)
#'   ) +
#'   geom_point() +
#'   geom_smooth(method = "lm") +
#'   theme_bw()
#'
#' parameters_gg(model = model, trait = "yield")
#'
#' @importFrom stats anova
parameters_gg <- function(model, trait = "trait") {
  summ <- summary(model)
  first_year <- min(model$model[, 2], na.rm = TRUE)
  last_year <- max(model$model[, 2], na.rm = TRUE)
  intercept <- summ$coefficients[1, 1]
  slope <- summ$coefficients[2, 1]
  se_slope <- summ$coefficients[2, 2]
  r2_lm <- summ$r.squared
  p_value <- as.data.frame(anova(model))$`Pr(>F)`[1]
  table_report <- data.frame(
    trait = trait,
    first_year = first_year,
    last_year = last_year,
    Slope = slope,
    se_Slope = se_slope,
    Intercept = intercept,
    r2 = r2_lm,
    `Pr(>F)` = p_value,
    `Genetic_Gain%` = slope / (intercept + first_year * slope) * 100,
    check.names = FALSE,
    row.names = "Regression"
  )
  return(table_report)
}
