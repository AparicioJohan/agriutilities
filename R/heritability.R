#' Heritability for Factor Analytic Models in ASReml-R
#'
#' @param model_fa Factor Analytic asreml model
#' @param genotype String
#' @param env String
#' @param vc.model variance covariance structure c("fa1", "fa2", "fa3", "fa4",
#' "us")
#'
#' @return data.frame
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
                            vc.model = c("fa2")) {
  G <- extractG(
    model = model_fa,
    gen = genotype,
    env = env,
    vc.model = vc.model
  )$VCOV
  Gvar <- mean(G[upper.tri(G, diag = FALSE)])
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
