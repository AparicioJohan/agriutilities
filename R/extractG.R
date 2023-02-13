#' Extract Variance-Covariance from ASReml-R
#'
#' @param model ASReml object
#' @param gen A character string indicating the column in data that contains
#' genotypes.
#' @param env A character string indicating the column in data that contains
#' environments or trials.
#' @param vc_model A character string indicating the variance-covariance fitted.
#' Can be 'diag', 'corv', 'corh', 'corgv', 'fa1', 'fa2', 'fa3', 'fa4', 'corgh',
#' 'us' or 'rr2'.
#' @return An object with a list of:
#' \item{VCOV}{A matrix of the estimated variance-covariance between trials.}
#' \item{CORR}{A n_trial x n_trial matrix with the correlation between trials.}
#' \item{vc_model}{A character string indicating the variance-covariance fitted.}
#' @export
#'
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#'
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_met(
#'  data = dat,
#'  genotype = "gen",
#'  trial = "county",
#'  traits = c("yield"),
#'  rep = "rep",
#'  block = "block",
#'  col = "col",
#'  row = "row"
#')
#'out <- single_trial_analysis(results, progress = FALSE)
#'met_results <- met_analysis(out, progress = FALSE)
#'extract_vcov(
#'  model = met_results$met_models$yield,
#'  vc_model = "us"
#')
#' }
#' @importFrom stats cov2cor
extract_vcov <- function(model = NULL,
                         gen = "genotype",
                         env = "trial",
                         vc_model = "corv") {
  sites <- data.frame(model$mf)[, env]
  s <- nlevels(sites)

  vc <- summary(model)$varcomp
  VCOV <- matrix(0, ncol = s, nrow = s)
  CORR <- matrix(0, ncol = s, nrow = s)
  diag(CORR) <- rep(1, s)

  gxe <- paste(env, gen, sep = ":")

  if (vc_model == "diag") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    diag(VCOV) <- vc[, 1]
  }
  if (vc_model == "corv") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    CORR <- matrix(1, ncol = s, nrow = s)
    CORR <- vc[1, 1] * CORR
    diag(CORR) <- rep(1, s)
    D <- rep(vc[2, 1], s)
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc_model == "corh") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    CORR <- matrix(1, ncol = s, nrow = s)
    CORR <- vc[1, 1] * CORR
    diag(CORR) <- rep(1, s)
    D <- vc[2:(s + 1), 1]
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc_model == "corgv") {
    vc.corr <- vc[grep(".cor", rownames(vc)), ]
    vc.var <- vc[-grep(".cor", rownames(vc)), ]
    k <- 1
    for (i in 1:s) {
      for (j in 1:i) {
        if (i != j) {
          CORR[i, j] <- vc.corr[k, 1]
          CORR[j, i] <- vc.corr[k, 1]
          k <- k + 1
        }
      }
    }
    D <- rep(vc.var[1, 1], s)
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc_model == "fa1") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    R <- vc.var[, 1]
    L <- vc.fa1[, 1]
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  if (vc_model == "fa2") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    vc.fa2 <- vc[grep("!fa2", rownames(vc)), ]
    R <- vc.var[, 1]
    L1 <- vc.fa1[, 1]
    L2 <- vc.fa2[, 1]
    L <- cbind(L1, L2)
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  if (vc_model == "fa3") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    vc.fa2 <- vc[grep("!fa2", rownames(vc)), ]
    vc.fa3 <- vc[grep("!fa3", rownames(vc)), ]
    R <- vc.var[, 1]
    L1 <- vc.fa1[, 1]
    L2 <- vc.fa2[, 1]
    L3 <- vc.fa3[, 1]
    L <- cbind(L1, L2, L3)
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  if (vc_model == "fa4") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    vc.fa2 <- vc[grep("!fa2", rownames(vc)), ]
    vc.fa3 <- vc[grep("!fa3", rownames(vc)), ]
    vc.fa4 <- vc[grep("!fa4", rownames(vc)), ]
    R <- vc.var[, 1]
    L1 <- vc.fa1[, 1]
    L2 <- vc.fa2[, 1]
    L3 <- vc.fa3[, 1]
    L4 <- vc.fa4[, 1]
    L <- cbind(L1, L2, L3, L4)
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  if (vc_model == "corgh") {
    vc.corr <- vc[grep(".cor", rownames(vc)), ]
    vc.var <- vc[-grep(".cor", rownames(vc)), ]
    k <- 1
    for (i in 1:s) {
      for (j in 1:i) {
        if (i != j) {
          CORR[i, j] <- vc.corr[k, 1]
          CORR[j, i] <- vc.corr[k, 1]
          k <- k + 1
        }
      }
    }
    D <- vc.var[1:s, 1]
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc_model == "us") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    k <- 1
    for (i in 1:s) {
      for (j in 1:i) {
        VCOV[i, j] <- vc[k, 1]
        k <- k + 1
      }
    }
    VCOV[upper.tri(VCOV)] <- t(VCOV)[upper.tri(VCOV)]
    CORR <- cov2cor(VCOV)
  }
  if (vc_model == "rr2") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    vc.fa2 <- vc[grep("!fa2", rownames(vc)), ]
    R <- vc.var[, 1]
    L1 <- vc.fa1[, 1]
    L2 <- vc.fa2[, 1]
    L <- cbind(L1, L2)
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  colnames(VCOV) <- levels(sites)
  colnames(CORR) <- levels(sites)
  rownames(VCOV) <- levels(sites)
  rownames(CORR) <- levels(sites)

  return(list(VCOV = VCOV, CORR = CORR, vc_model = vc_model))
}
