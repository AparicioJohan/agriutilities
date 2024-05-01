#' @title Find the AIC and BIC for a set of models fitted using asreml
#' @description The function is a wrapper for 'icREML' function described
#' in Verbyla (2019).
#' @param fm A \code{list} of asreml fitted model objects
#' @param scale A scalar to scale the variance matrix of the estimated
#' fixed effects (to ensure numerical stability of a log-determinant).
#' Default value is 1.
#' @return A data frame.  The data frame has the following components
#' \itemize{
#' \item \code{model} : the names of the models
#' \item \code{loglik} : the full log-likelihood for each model
#' \item \code{p} :  the number of fixed effects parameters for each model
#' \item \code{q} : the number of (non-zero) variance parameters for each model.
#' \item \code{b} : the number of variance parameters that are fixed or on the
#' boundary.  These parameters are not counted in the AIC or BIC.
#' \item \code{AIC} : the AIC for each model
#' \item \code{BIC} : the BIC for each model
#' \item \code{logdet} : the log-determinant used in adjusting the residual
#' log-likelihood for each model
#' }
#' @export
#'
#' @author Ari Verbyla (averbyla at avdataanalytics.com.au)
#' @references
#' Verbyla, A. P. (2019). A note on model selection using information
#' criteria for general linear models estimated using REML. Australian &
#' New Zealand Journal of Statistics, 61(1), 39-50.
ic_reml_asr <- function(fm, scale = 1, logdet = TRUE) {
  if (!is.list(fm)) stop(" Models need to be in a list\n")
  if (is.null(names(fm))) {
    namesfm <- paste("fm", 1:length(fm))
  } else {
    namesfm <- names(fm)
  }
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop("The package asreml is not loaded.")
  }
  asreml::asreml.options(Cfixed = TRUE, gammaPar = FALSE)
  fm <- lapply(fm, function(el) {
    if (is.null(el$Cfixed)) {
      out <- asreml::update.asreml(el, maxit = 1)
    } else {
      out <- el
    }
    out
  })
  logl <- lapply(fm, function(el) el$loglik)
  summ <- lapply(fm, function(el) summary(el, coef = TRUE)$coef.fixed)
  which.X0 <- lapply(summ, function(el) !is.na(el[, "z.ratio"]))
  p.0 <- lapply(which.X0, function(el) sum(el))
  Cfixed <- lapply(fm, function(el) {
    ord <- rownames(summary(el, coef = TRUE)$coef.fixed)
    el$Cfixed[ord, ord, drop = FALSE]
  })
  logdetC <- lapply(
    X = 1:length(fm),
    FUN = function(el, Cfixed, which.X0, scale) {
      mysvdd <- svd(as.matrix(scale * Cfixed[[el]][which.X0[[el]], which.X0[[el]]]))$d
      mysvdd <- mysvdd[mysvdd > 0]
      sum(log(mysvdd))
    }, Cfixed, which.X0, scale
  )
  vparam <- lapply(fm, function(el) summary(el)$varcomp)
  q.0 <- lapply(vparam, function(el) {
    sum(!(el$bound == "F" | el$bound == "B" | el$bound == "C")) +
      sum(el$bound[!is.na(str_extract(dimnames(el)[[1]], "cor"))] == "B")
  })
  b.0 <- lapply(vparam, function(el) {
    sum(el$bound == "F" | el$bound == "B") -
      sum(el$bound[!is.na(str_extract(dimnames(el)[[1]], "cor"))] == "B")
  })
  full.logl <- lapply(1:length(fm), function(el, logl, logdetC, p.0) {
    logl[[el]] - logdetC[[el]] / 2
  }, logl, logdetC, p.0)
  aic <- unlist(lapply(1:length(fm), function(el, full.logl, p.0, q.0) {
    -2 * full.logl[[el]] + 2 * (p.0[[el]] + q.0[[el]])
  }, full.logl, p.0, q.0))
  bic <- unlist(lapply(
    1:length(fm), function(el, full.logl, p.0, q.0, fm) {
      -2 * full.logl[[el]] + log(fm[[el]]$nedf + p.0[[el]]) * (p.0[[el]] + q.0[[el]])
    },
    full.logl, p.0, q.0, fm
  ))
  results <- data.frame(
    model = namesfm,
    res_loglik = unlist(logl),
    full_loglik = unlist(full.logl),
    p = unlist(p.0), q = unlist(q.0), b = unlist(b.0),
    AIC = aic, BIC = bic, logdet = unlist(logdetC)
  )
  row.names(results) <- 1:dim(results)[1]
  return(results)
}


#' @title Find the AIC and BIC for a model fitted using SpATS
#' @description This function calculates the AIC and BIC for a model fitted in
#' SpATS following the methodology proposed by Verbyla (2019).
#' @param model A model fitted using SpATS.
#' @param scale A scalar to scale the variance matrix of the estimated
#' fixed effects (to ensure numerical stability of a log-determinant).
#' Default value is 1.
#' @param k An integer value to round ratios when identifying boundary variance
#' parameters. Default value is 2.
#' @param label A string to label the model. Default value is "spats".
#' @return A data frame.  The data frame has the following components
#' \itemize{
#' \item \code{model} : the name of the models
#' \item \code{loglik} : the full log-likelihood for each model
#' \item \code{p} :  the number of fixed effects parameters for each model
#' \item \code{q} : the number of (non-zero) variance parameters for each model.
#' \item \code{b} : the number of variance parameters that are fixed or on the
#' boundary.  These parameters are not counted in the AIC or BIC.
#' \item \code{AIC} : the AIC for each model
#' \item \code{BIC} : the BIC for each model
#' \item \code{logdet} : the log-determinant used in adjusting the residual
#' log-likelihood for each model
#' }
#' @export
#'
#' @author Johan Aparicio
#' @references
#' Verbyla, A. P. (2019). A note on model selection using information
#' criteria for general linear models estimated using REML. Australian &
#' New Zealand Journal of Statistics, 61(1), 39-50.
ic_reml_spt <- function(model, scale = 1, k = 2, label = "spats") {
  loglik <- model$deviance / -2
  fixed_eff <- model$coeff[!attr(model$coeff, "random")]
  p_0 <- length(fixed_eff)
  nedf <- model$nobs - p_0 # obs - fixed effects
  cmat <- model$vcov$C22_inv[names(fixed_eff), names(fixed_eff)]
  mysvdd <- svd(as.matrix(scale * cmat))$d
  mysvdd <- mysvdd[mysvdd > 0]
  logdetC <- sum(log(mysvdd))
  q_0 <- length(model$var.comp) + 1
  eff_dim <- model$eff.dim
  ids <- names(eff_dim)[!names(eff_dim) %in% names(fixed_eff)]
  b_0 <- sum(round(eff_dim[ids] / model$dim.nom[ids], k) == 0)
  q_0 <- q_0 - b_0
  full_logl <- loglik - logdetC / 2
  aic <- -2 * logl + 2 * (p_0 + q_0)
  bic <- -2 * logl + log(nedf + p_0) * (p_0 + q_0)
  results <- data.frame(
    model = label,
    res_loglik = loglik,
    full_loglik = full_logl,
    p = p_0, q = q_0, b = b_0,
    AIC = aic, BIC = bic,
    logdet = logdetC
  )
  return(results)
}
