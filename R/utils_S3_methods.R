#' Print an object of class \code{checkAgri}
#'
#' @description Prints information about \code{check_design_MET()} function.
#'
#' @aliases print.checkAgri
#' @usage \method{print}{checkAgri}(x, ...)
#' @param x An object fitted with the function \code{check_design_MET()}.
#' @param ... Options used by the tibble package to format the output. See
#' `tibble::print()` for more details.
#' @author Johan Aparicio [aut]
#' @method print checkAgri
#' @export
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_MET(
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
#' }
print.checkAgri <- function(x, ...) {
  cat("---------------------------------------------------------------------\n")
  cat("Summary Traits by Trial:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$summ_traits, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("Experimental Design Detected:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$exp_design_list)
  cat("\n---------------------------------------------------------------------\n")
  cat("Summary Experimental Design:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$exp_design_resum[, 1:9], ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("Connectivity Matrix:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$connectivity_matrix)
  cat("\n")
}

#' Print an object of class \code{smaAgri}
#'
#' @description Prints information about \code{single_trial_analysis} function.
#'
#' @aliases print.smaAgri
#' @usage \method{print}{smaAgri}(x, ...)
#' @param x An object fitted with the function \code{single_trial_analysis()}.
#' @param ... Options used by the tibble package to format the output. See
#' `tibble::print()` for more details.
#' @author Johan Aparicio [aut]
#' @method print smaAgri
#' @importFrom utils head
#' @export
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_MET(
#'   data = dat,
#'   genotype = "gen",
#'   trial = "county",
#'   traits = c("yield"),
#'   rep = "rep",
#'   block = "block",
#'   col = "col",
#'   row = "row"
#' )
#' out <- single_trial_analysis(results, progress = FALSE)
#' print(out)
#' }
print.smaAgri <- function(x, ...) {
  cat("---------------------------------------------------------------------\n")
  cat("Summary Fitted Models:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$resum_fitted_model, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("Outliers Removed:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$outliers, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("First Predicted Values and Standard Errors (BLUEs/BLUPs):\n")
  cat("---------------------------------------------------------------------\n")
  print(head(x$blues_blups))
  cat("\n")
}

#' Print an object of class \code{metAgri}
#'
#' @description Prints information about \code{met_analysis()} function.
#'
#' @aliases print.metAgri
#' @usage \method{print}{metAgri}(x, ...)
#' @param x An object fitted with the function \code{met_analysis()}.
#' @param ... Options used by the tibble package to format the output. See
#' `tibble::print()` for more details.
#' @author Johan Aparicio [aut]
#' @method print metAgri
#' @importFrom utils head
#' @export
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_MET(
#'   data = dat,
#'   genotype = "gen",
#'   trial = "county",
#'   traits = c("yield"),
#'   rep = "rep",
#'   block = "block",
#'   col = "col",
#'   row = "row"
#' )
#' out <- single_trial_analysis(results, progress = FALSE)
#' met_results <- met_analysis(out)
#' print(met_results)
#' }
print.metAgri <- function(x, ...) {
  cat("---------------------------------------------------------------------\n")
  cat("Trial Effects (BLUEs):\n")
  cat("---------------------------------------------------------------------\n")
  print(x$trial_effects, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("Heritability:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$heritability, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("First Overall Predicted Values and Standard Errors (BLUPs):\n")
  cat("---------------------------------------------------------------------\n")
  print(head(x$overall_BLUPs))
  cat("\n---------------------------------------------------------------------\n")
  cat("Variance-Covariance Matrix:\n")
  cat("---------------------------------------------------------------------\n")
  lt <- names(x$VCOV)
  for (i in lt) {
    vc_str <- x$VCOV[[i]]$vc.model
    cat("\nCorrelation Matrix ", "('", vc_str,"'): ", i, "\n", sep = "")
    print(round(x$VCOV[[i]]$CORR, 2))
    cat("\nCovariance Matrix ", "('", vc_str,"'): ", i, "\n", sep = "")
    print(round(x$VCOV[[i]]$VCOV, 2))
    cat("\n---------------------------------------------------------------------")
  }
  cat("\nFirst Stability Coefficients:\n")
  cat("---------------------------------------------------------------------\n")
  print(head(x$stability))
  cat("\n")
}
