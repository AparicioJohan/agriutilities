#' Genetic Gain parameters
#'
#' @param model linear regression model (lm object)
#' @param trait string with trait evaluated to be added in the data.frame
#'
#' @return data.frame with some parameters from the regression
#' (Slope, se_Slope, Intercept, r2, Pr(>F) and Genetic_Gain%)
#' @export
#'
#' @examples
#' # parameters_GG(model)
parameters_GG <- function(model, trait = "trait") {
  summ <- summary(model)
  first_year <- min(model$model[, 2], na.rm = TRUE)
  last_year <- max(model$model[, 2], na.rm = TRUE)
  intercept <- summ$coefficients[1, 1]
  slope <- summ$coefficients[2, 1]
  se_slope <- summ$coefficients[2, 2]
  r2_lm <- summ$r.squared
  p_value <- as.data.frame(anova(model2))$`Pr(>F)`[1]
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
