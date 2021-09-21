
#' extractG
#'
#' @param model ASReml object
#' @param gen string with genotypes
#' @param env string with genotypes
#' @param vc.model variance covariance that you fitted
#'
#' @return list VCOV = VCOV , CORR = CORR, vc.model = vc.model
#' @export
#'
#' @examples
#' # extractG(model, gen = "genotype" , env = "trial", vc.model = "corv")
#' @importFrom stats cov2cor
extractG <- function(model, gen = "genotype" , env = "trial", vc.model = "corv"){

  sites <- data.frame(model$mf)[, env]
  s <- nlevels(sites)

  vc <- summary(model)$varcomp
  VCOV <- matrix(0, ncol = s, nrow = s)
  CORR <- matrix(0, ncol = s, nrow = s)
  diag(CORR) <- rep(1, s)

  gxe <- paste(env, gen, sep = ":")

  if (vc.model == "diag") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    diag(VCOV) <- vc[, 1]
  }
  if (vc.model == "corv") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    CORR <- matrix(1, ncol = s, nrow = s)
    CORR <- vc[1, 1] * CORR
    diag(CORR) <- rep(1, s)
    D <- rep(vc[2, 1], s)
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc.model == "corh") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    CORR <- matrix(1, ncol = s, nrow = s)
    CORR <- vc[1, 1] * CORR
    diag(CORR) <- rep(1, s)
    D <- vc[2:(s + 1), 1]
    VCOV <- diag(sqrt(D)) %*% CORR %*% diag(sqrt(D))
  }
  if (vc.model == "corgv") {
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
  if (vc.model == "fa1") {
    vc.var <- vc[grep("!var", rownames(vc)), ]
    vc.fa1 <- vc[grep("!fa1", rownames(vc)), ]
    R <- vc.var[, 1]
    L <- vc.fa1[, 1]
    VCOV <- L %*% t(L) + diag(R)
    CORR <- cov2cor(VCOV)
  }
  if (vc.model == "fa2") {
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
  if (vc.model == "fa3") {
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
  if (vc.model == "fa4") {
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
  if (vc.model == "corgh") {
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
  if (vc.model == "us") {
    vc <- vc[grep(gxe, rownames(vc)), ]
    k <- 1
    for (i in 1:s) {
      for (j in 1:i) {
        VCOV[i,j] <- vc[k,1]
        k <- k+1
      }
    }
    VCOV[upper.tri(VCOV)] = t(VCOV)[upper.tri(VCOV)]
    CORR <- cov2cor(VCOV)
  }
  if (vc.model == "rr2") {
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

  return(list(VCOV = VCOV , CORR = CORR, vc.model = vc.model))

}
