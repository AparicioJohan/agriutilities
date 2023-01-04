#' Heritability for Factor Analytic Models in ASReml-R
#'
#' @param model_fa Factor Analytic asreml model
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param env A character string indicating the column in data that contains
#' environments or trials.
#' @param vc.model A character string indicating the variance-covariance
#' structure. Can be "fa1", "fa2", "fa3", "fa4" or "us".
#' @param diag TRUE or FALSE depending on the user if they want to take the
#' elements on the diagonal of the variance-covariance matrix or the elements
#' out of the diagonal to estimate the heritability. FALSE by default.
#'
#' @return A list
#' @export
#'
#' @examples
#' # library(tidyverse)
#' # library(asreml)
#' # library(agridat)
#' # library(agriutilities)
#' # data(besag.met)
#' # dat <- besag.met
#' #
#' # dat <- dat %>% arrange(county)
#' # model <- asreml(
#' #   fixed = yield ~ 1 + county,
#' #   random = ~ fa(county, 2):gen + county:rep + diag(county):rep:block,
#' #   residual = ~ dsum(~ units | county),
#' #   data = dat,
#' #   na.action = list(x="include",y="include")
#' # )
#' # heritability_fa(
#' #  model_fa = model,
#' #  genotype = "gen",
#' #  env = "county",
#' #  vc.model = "fa2"
#' # )
heritability_fa <- function(model_fa = NULL,
                            genotype = "line",
                            env = "loc",
                            vc.model = c("fa2"),
                            diag = FALSE) {
  G <- extractG(
    model = model_fa,
    gen = genotype,
    env = env,
    vc.model = vc.model
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
  H2Cullis <- 1 - (vdBLUP.avg / 2 / Gvar)
  H2Cullis
  avsed <- pr$pvals %>%
    summarise(mean = mean(std.error^2)) %>%
    as.numeric()
  h2_js <- 1 - avsed / Gvar
  return(list(H2Cullis = H2Cullis, h2_js = h2_js))
}
