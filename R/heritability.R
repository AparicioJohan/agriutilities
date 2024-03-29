#' Heritability for Factor Analytic Models in ASReml-R
#'
#' @param model_fa Factor Analytic ASReml model
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param env A character string indicating the column in data that contains
#' environments or trials.
#' @param vc_model A character string indicating the variance-covariance
#' structure. Can be "fa1", "fa2", "fa3", "fa4" or "us".
#' @param diag \code{TRUE} or \code{FALSE} depending on the user if they want to
#'  take the elements on the diagonal of the variance-covariance matrix or the
#'  elements out of the diagonal to estimate the heritability. \code{FALSE} by
#'  default.
#'
#' @return An object with a list of:
#' \item{h2_cullis}{A numerical value of the Cullis heritability estimate.}
#' \item{h2_se}{A numerical value of the Cullis heritability estimate based on
#'  the standard error.}
#' @export
#'
#' @examples
#' \dontrun{
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
#' out <- single_trial_analysis(results, progress = FALSE)
#' met_results <- met_analysis(out, progress = FALSE)
#' model <- met_results$met_models$yield
#' heritability_fa(
#'    model_fa = model,
#'    genotype = "genotype",
#'    env = "trial",
#'    vc_model = "us"
#'  )
#' }
heritability_fa <- function(model_fa = NULL,
                            genotype = "line",
                            env = "loc",
                            vc_model = c("fa2"),
                            diag = FALSE) {
  G <- extract_vcov(
    model = model_fa,
    gen = genotype,
    env = env,
    vc_model = vc_model
  )$VCOV

  if (diag) {
    Gvar <- mean(diag(G))
  } else {
    Gvar <- mean(G[upper.tri(G, diag = FALSE)])
  }
  pr <- suppressWarnings(
    asreml::predict.asreml(model_fa, classify = genotype, sed = TRUE, trace = 0)
  )
  vdBLUP.mat <- pr$sed^2
  vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag = FALSE)])
  h2_cullis <- 1 - (vdBLUP.avg / 2 / Gvar)
  h2_cullis
  avsed <- pr$pvals %>%
    summarise(mean = mean(std.error^2)) %>%
    as.numeric()
  h2_js <- 1 - avsed / Gvar
  return(list(h2_cullis = h2_cullis, h2_se = h2_js))
}
