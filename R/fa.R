

#' Factor analytic variance models.
#'
#' @description Summarises an asreml model in which factor analytic (FA) variance structures have been fitted. The function
#' aids in interpreting FA models fitted to genotype by environment effects in the analysis
#' of plant breeding multi-environment trial (MET) data.
#'
#' @param object An asreml object.
#' @param uniplot if TRUE (the default), pairwise plots of FA loadings are produced for k > 1.
#' @param uniplot.tol A numeric value between 0 and 1 that affects the display of pairs of loadings in the uniplots. Those pairs with a (scaled) euclidean distance greater than uniplot.tol are drawn with solid blue lines, otherwise they appear as red dashed lines.
#' @param uniplot.cex Character expansion (relative) size for plotting (the default is 0.5).
#' @param trunc.char A numeric vector of length 2 specifying the truncation of loading labels for use in plotting. The default (c(1,6)) is to use the first 6 characters of the environment factor levels, otherwise if NULL the complete level names are used.
#' @param blups 	If TRUE (the default), summaries of estimates from the asreml object are produced. These include BLUPs of genotype by environment effects and genotype scores, and REML estimates of FA variance parameters.
#' @param regplot If TRUE (the default), regression plots of genotype by environment effects against FA loadings are produced.
#' @param addedplot If TRUE (the default), added variable plots of genotype by environment effects against FA loadings are produced.
#' @param g.list A character vector comprising the subset of genotype levels to be used in the regression and added variable plots. The default (NULL) is to use the complete set of genotype levels.
#' @param heatmap If TRUE (the default), a heatmap of the estimated FA correlation matrix is produced.
#' @param heatmap.ord A character string ("asis", "cluster") or vector of environment levels in the order desired for plotting the heatmap. The default is "asis", where the rows and columns of the correlation matrix are maintained in environment level order. The "cluster" option invokes the agnes clustering algorithm on the FA correlation matrix then orders the rows and columns accordingly.
#' @param agnes.method Clustering strategy; the default is "average".
#'
#' @details The following documentation is taken from the ASExtras4 package. Please refer to this package and associated publications if you are interested in going deeper on this technique. You may be interested in reading and citing this publication if using this methodology
#'
#' @references Cullis BR, Smith AB, Beek C, Cowling WA (2010). “Analysis of Yield and Oil from a Series of Canola Breeding Trials. Part II: Exploring VxE using Factor Analysis.” Genome, 53, 1002-1016. Smith AB, Cullis BR, Thompson R (2001). “Analyzing Variety by Environment Data Using Multiplicative Mixed Models and Adjustments for Spatial Field Trend.” Biometrics, 57, 1138-1147.
#'
#' @return list
#' @export
#'
#' @examples
#' # library(tidyverse)
#' # library(asreml)
#' # library(agridat)
#' # data(besag.met)
#' # dat <- besag.met
#' #
#' # dat <- dat %>% arrange(county)
#' # model <- asreml(fixed = yield ~ 1 + county,
#' #                 random = ~ fa(county, 2):gen + county:rep + diag(county):rep:block,
#' #                 residual = ~ dsum(~ units | county),
#' #                 data = dat,
#' #                 na.action = list(x="include",y="include"))
#' #
#' # ASM <- fa.asreml(model, trunc.char = NULL)
fa.asreml <- function (object, uniplot = F, uniplot.tol = 0.85, uniplot.cex = 0.5,
                       trunc.char = c(1, 6), blups = TRUE, regplot = F, addedplot = F,
                       g.list = NULL, heatmap = F, heatmap.ord = "asis", agnes.method = "average")
{
  fa.var <- function(fat) {
    unlist(lapply(fat, function(x) {
      if (x$Fun == "fa")
        x$FacNam
      else NULL
    }))
  }
  fa.outer <- function(fat) {
    unlist(lapply(fat, function(x) {
      if (x$Fun == "fa")
        x$Obj
      else NULL
    }))
  }
  fa.inner.name <- function(fat) {
    unlist(lapply(fat, function(x) {
      if (x$Fun != "fa")
        x$FacNam
      else NULL
    }))
  }
  fa.inner.var <- function(fat) {
    unlist(lapply(fat, function(x) {
      if (x$Fun != "fa")
        x$Obj
      else NULL
    }))
  }
  fa.inner.fun <- function(fat) {
    unlist(lapply(fat, function(x) {
      if (x$Fun != "fa")
        x$Fun
      else NULL
    }))
  }
  ide.order <- function(fat) {
    ide <- which(unlist(lapply(fat, function(x) {
      as.logical(sum(unlist(lapply(x, function(y) {
        if (y$Fun == "ide") TRUE else FALSE
      }))))
    })))
    if (length(ide) > 0)
      c(seq(along = fat)[-ide], ide)
    else seq(along = fat)
  }
  if (!inherits(object, "asreml"))
    stop("\nObject must be of class 'asreml'\n")
  if (!requireNamespace("asreml", quietly = TRUE))
    stop("Requires package 'asreml'")
  tt <- attr(object$mf, "model.terms")$random$Terms.obj
  if (length(which.fun <- attr(tt, "specials")$fa) == 0)
    stop("No fa() term in model\n")
  ttf <- attr(tt, "factors")
  tt.vars <- dimnames(ttf)[[1]]
  term.labels <- attr(tt, "term.labels")
  which.fa <- list()
  names.fa <- NULL
  for (i in which.fun) {
    x <- which(ttf[i, ] > 0)
    for (j in x) {
      which.fa <- c(which.fa, list(which(ttf[, j] > 0)))
      names.fa <- c(names.fa, term.labels[j])
    }
  }
  nterms <- length(which.fa)
  faterms <- vector(mode = "list", length = nterms)
  for (w in 1:nterms) {
    ww <- which.fa[[w]]
    faterms[[w]] <- lapply(attr(object$mf, "model.terms")$random$Vars[tt.vars[ww]],
                           function(x) {
                             y <- list(x$Fun, x$FacNam, x$Obj)
                             names(y) <- c("Fun", "FacNam", "Obj")
                             y
                           })
  }
  names(faterms) <- names.fa
  idx <- ide.order(faterms)
  nice.sum <- summary(object, vparameters = TRUE)$vparameters
  if (is.null(mf <- object$mf)) {
    if (is.null(object$RDS)) {
      stop(object, "is missing model frame and has no RDS component.")
    }
    else mf <- readRDS(object$RDS)
  }
  if (length(y <- attr(mf, "traits")$lhs) > 1)
    stop("No method for multivariate")
  y <- mf[[y]]
  total.where <- nterms
  gammas.lst <- list()
  fanam <- names(faterms)
  uniplot.lst <- blup.lst <- regplot.lst <- addedplot.lst <- heatmap.lst <- agnes.lst <- NULL
  if (uniplot)
    uniplot.lst <- list()
  if (blups)
    blup.lst <- list()
  if (regplot)
    regplot.lst <- list()
  if (addedplot)
    addedplot.lst <- list()
  if (heatmap) {
    heatmap.lst <- list()
    if (heatmap.ord == "cluster") {
      agnes.lst <- list()
    }
  }
  for (nt in idx) {
    if (length(fa.inner.var(faterms[[nt]])) > 1)
      stop("Only first order interaction with FA term allowed\n")
    Variety <- factor(as.character(mf[[fa.inner.var(faterms[[nt]])]]))
    Site <- mf[[fa.outer(faterms[[nt]])]]
    outer.name <- fa.outer(faterms[[nt]])
    inner.name <- fa.inner.var(faterms[[nt]])
    snam <- levels(Site)
    ll <- min(nchar(snam))
    if (length(trunc.char))
      sn <- substring(snam, min(ll, trunc.char[1]), max(ll,
                                                        trunc.char[2]))
    else sn <- snam
    if (length(unique(sn)) < length(snam))
      stop(paste("Fewer levels in FA term than expected,\n",
                 " set 'trunc.char' to avoid non-unique level names in factor",
                 outer.name))
    vnam <- levels(Variety)
    ns <- length(snam)
    nv <- length(vnam)
    if (ns == 0)
      stop(paste("\n", fa.outer(faterms[[nt]]), " is not a factor\n",
                 sep = ""))
    if (nv == 0)
      stop(paste("\n", fa.inner.var(faterms[[nt]]), " is not a factor\n",
                 sep = ""))
    Variety.fac <- Variety
    if (length(grep(":", unique(Variety.fac))))
      stop(paste("Levels of", inner.name, "cannot contain the ':' character"))
    if (length(grep(":", unique(Site))))
      stop(paste("Levels of", outer.name, "cannot contain the ':' character"))
    gammas <- matrix(nice.sum[[names(faterms)[nt]]], nrow = ns)
    k <- ncol(gammas) - 1
    psi <- gammas[, 1]
    names(psi) <- snam
    Lam <- gammas[, -1, drop = FALSE]
    if (k > 1) {
      ss <- svd(Lam)
      Lam <- -Lam %*% ss$v
    }
    dimnames(Lam) <- list(snam, paste("fac", 1:k, sep = "_"))
    Gmat <- Lam %*% t(Lam) + diag(psi)
    Cmat <- cov2cor(Gmat)
    paf.site <- matrix(0, nrow = ns, ncol = k)
    dimnames(paf.site) <- list(snam, paste("fac", 1:k, sep = "_"))
    for (i in 1:k) {
      paf.site[, i] <- 100 * diag(Lam[, i] %*% t(Lam[,
                                                     i]))/diag(Gmat)
    }
    if (k > 1) {
      all <- 100 * diag(Lam %*% t(Lam))/diag(Gmat)
      paf.site <- cbind(paf.site, all)
    }
    paf.mod <- 100 * sum(diag(Lam %*% t(Lam)))/sum(diag(Gmat))
    dd <- 1/sqrt(diag(Gmat))
    Lamc <- diag(dd) %*% Lam
    dimnames(Lamc) <- dimnames(Lam)
    gammas.lst[[nt]] <- list(Gmat = Gmat, Cmat = Cmat, `site %vaf` = paf.site,
                             `total %vaf` = paf.mod, `rotated loads` = Lam, `specific var` = psi,
                             `rotated loads - c` = Lamc)
    if (k == 1)
      uniplot <- FALSE
    if (uniplot) {
      bp.main <- names(faterms)[nt]
      if (k == 2) {
        uniplot.lst[[nt]] <- xyplot(Lamc[, 1] ~ Lamc[,
                                                     2], xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1),
                                    asp = "s", xlab = "Loading 2", ylab = "Loading 1",
                                    main = bp.main, panel = function(x, y, subscripts,
                                                                     sn, tol, lcex, ...) {
                                      panel.curve(sqrt(1 - x^2), from = -1, to = 1)
                                      panel.curve(-sqrt(1 - x^2), from = -1, to = 1)
                                      ne <- length(sn)
                                      radius <- sqrt(x * x + y * y)
                                      lcol <- rep("blue", ne)
                                      ltyp <- rep(1, ne)
                                      llwd <- rep(1.5, ne)
                                      lcol[radius < tol] <- "red"
                                      ltyp[radius < tol] <- 3
                                      llwd[radius < tol] <- 1
                                      panel.segments(rep(0, ne), rep(0, ne), x,
                                                     y, lty = ltyp, col = lcol, lwd = llwd)
                                      ltext(x, y, sn, cex = lcex, col = lcol,
                                            srt = 45)
                                    }, sn = sn, tol = uniplot.tol, lcex = uniplot.cex)
      }
      else {
        tmp <- vector(mode = "list", length = k * (k -
                                                     1)/2)
        natmp <- character(0)
        for (i in seq(1, k - 1)) {
          for (j in seq(i + 1, k)) {
            natmp <- c(natmp, paste("Ld", i, j, sep = ""))
            tmp[[(i - 1) * (k - 1) + (j - i)]] <- xyplot(Lamc[,
                                                              i] ~ Lamc[, j], xlim = c(-1.1, 1.1), ylim = c(-1.1,
                                                                                                            1.1), asp = "s", xlab = paste("Loading",
                                                                                                                                          j), ylab = paste("Loading", i), main = bp.main,
                                                         panel = function(x, y, subscripts, sn,
                                                                          tol, lcex, ...) {
                                                           panel.curve(sqrt(1 - x^2), from = -1,
                                                                       to = 1)
                                                           panel.curve(-sqrt(1 - x^2), from = -1,
                                                                       to = 1)
                                                           ne <- length(sn)
                                                           radius <- sqrt(x * x + y * y)
                                                           lcol <- rep("blue", ne)
                                                           ltyp <- rep(1, ne)
                                                           llwd <- rep(1.5, ne)
                                                           lcol[radius < tol] <- "red"
                                                           ltyp[radius < tol] <- 3
                                                           llwd[radius < tol] <- 1
                                                           panel.segments(rep(0, ne), rep(0, ne),
                                                                          x, y, lty = ltyp, col = lcol, lwd = llwd)
                                                           ltext(x, y, sn, cex = lcex, col = lcol,
                                                                 srt = 45)
                                                         }, sn = sn, tol = uniplot.tol, lcex = uniplot.cex)
          }
        }
        names(tmp) <- natmp
        uniplot.lst[[nt]] <- tmp
      }
    }
    if (heatmap) {
      dimnames(Cmat) <- list(sn, sn)
      if (heatmap.ord == "asis") {
        sn.ord <- sn
      }
      else if (heatmap.ord == "cluster") {
        dis.mat <- 1 - Cmat
        agnes.lst[[nt]] <- agnes(x = dis.mat, diss = TRUE,
                                 method = agnes.method)
        sn.ord <- agnes.lst[[nt]]$order.lab
      }
      else {
        sn.ord <- heatmap.ord
      }
      Cmat.ord <- Cmat[rev(sn.ord), rev(sn.ord)]
      diag(Cmat.ord) <- NA
      range(Cmat.ord, na.rm = T)
      atss <- seq(-1, 1, 0.1)
      hh <- rev(rainbow(256, start = 0, end = 2/3))
      hm.lab <- fa.outer(faterms[[nt]])
      hm.main <- names(faterms)[nt]
      heatmap.lst[[nt]] <- levelplot(Cmat.ord[sn.ord,
                                              rev(sn.ord)], at = atss, col.regions = hh, scales = list(x = list(rot = 60,
                                                                                                                cex = 0.9), y = list(cex = 0.9)), xlab = hm.lab,
                                     ylab = hm.lab, main = hm.main)
    }
    if (blups) {
      cc <- coef(object, list = TRUE)[[names(faterms)[nt]]]
      blup.df <- data.frame(blup = as.vector(cc), stringsAsFactors = TRUE)
      nn <- dimnames(cc)[[1]]
      temp <- strsplit(nn, split = ":", fixed = TRUE)
      tt.inner <- sapply(temp, function(x) x[2])
      tt.outer <- sapply(temp, function(x) x[1])
      tt.outer <- sapply(strsplit(tt.outer, split = "_",
                                  fixed = T), function(x) paste(x[-1], collapse = "_"))
      tt.inner <- sapply(strsplit(tt.inner, split = "_",
                                  fixed = T), function(x) paste(x[-1], collapse = "_"))
      blup.df[[outer.name]] <- tt.outer
      blup.df[[inner.name]] <- tt.inner
      score.df <- subset(blup.df, is.element(blup.df[[outer.name]],
                                             paste("Comp", 1:k, sep = "")))
      score.mat <- matrix(score.df$blup, ncol = k)
      if (k > 1) {
        score.mat <- -score.mat %*% ss$v
      }
      score.df$blupr <- as.vector(score.mat)
      blup.df <- subset(blup.df, !is.element(blup.df[[outer.name]],
                                             paste("Comp", 1:k, sep = "")))
      blup.df$regblup <- as.vector(score.mat %*% t(Lam))
      blup.inmet.df <- subset(blup.df, is.element(blup.df[[inner.name]],
                                                  unique(Variety.fac)))
      score.inmet.df <- subset(score.df, is.element(score.df[[inner.name]],
                                                    unique(Variety.fac)))
      blup.inmet.df <- blup.inmet.df[order(type.convert(blup.inmet.df[[outer.name]]),
                                           type.convert(blup.inmet.df[[inner.name]])),
                                     ]
      pres <- tapply(y, list(Variety, Site), function(x) length(x[!is.na(x)]))
      blup.inmet.df$pres <- as.vector(pres)
      score.inmet.df <- score.inmet.df[order(type.convert(score.inmet.df[[outer.name]]),
                                             type.convert(score.inmet.df[[inner.name]])),
                                       ]
      blup.lst[[nt]] <- list(blups = blup.df, scores = score.df,
                             blups.inmet = blup.inmet.df, scores.inmet = score.inmet.df)
    }
    if (regplot) {
      if (!blups) {
        stop("'blups' must be set to TRUE for regplots\n")
      }
      score.mat <- matrix(score.inmet.df$blupr, ncol = k)
      dimnames(score.mat) <- list(unique(score.inmet.df[[inner.name]]),
                                  paste("fac", 1:k, sep = "_"))
      if (is.null(g.list)) {
        g.list <- list()
        dd <- dimnames(score.mat)[[1]]
        for (kk in 1:k) {
          g.list[[kk]] <- dd[order(score.mat[, kk])][1]
          g.list[[kk + k]] <- dd[order(-score.mat[,
                                                  kk])][1]
        }
        g.list <- unlist(g.list)
        g.list <- unique(g.list)
        if (length(g.list) != 2 * k) {
          tt <- table(Variety)
          tt <- tt[!is.element(names(tt), g.list)]
          tt <- names(tt)[order(-tt)]
          tt <- tt[1:(2 * k - length(g.list))]
          g.list <- c(g.list, tt)
        }
      }
      score.mat <- score.mat[sort(g.list), , drop = FALSE]
      regplot.df <- subset(blup.inmet.df, is.element(blup.inmet.df[[inner.name]],
                                                     g.list))
      for (kk in 1:k) {
        facnam <- paste("fac", kk, sep = ".")
        regplot.df[[facnam]] <- rep(Lam[, kk], each = length(g.list))
      }
      regplot.df[[inner.name]] <- factor(regplot.df[[inner.name]])
      xform <- paste(paste("fac", 1:k, sep = "."), collapse = "+")
      xform <- formula(paste("blup~", xform, "|", inner.name,
                             sep = ""))
      rp.main <- paste(names(faterms)[nt], "BLUPS")
      regplot.lst[[nt]] <- xyplot(xform, data = regplot.df,
                                  outer = TRUE, as.table = TRUE, par.strip.text = list(cex = 0.6),
                                  main = rp.main, panel = function(x, y, slopes,
                                                                   ...) {
                                    panel.xyplot(x, y)
                                    panel.abline(b = slopes[current.column(),
                                                            current.row()], a = 0)
                                  }, slopes = score.mat)
      if (addedplot) {
        Yadd <- list()
        for (kk in 1:k) {
          Lamk <- Lam
          Lamk[, kk] <- 0
          Yadd[[kk]] <- score.mat %*% t(Lamk)
        }
        xform <- paste(paste("fac", 1:k, sep = "."),
                       collapse = "+")
        xform <- formula(paste("blup~", xform, "|",
                               inner.name, sep = ""))
        ap.main <- paste(names(faterms)[nt], "BLUPS")
        addedplot.lst[[nt]] <- xyplot(xform, data = regplot.df,
                                      outer = T, as.table = T, par.strip.text = list(cex = 0.6),
                                      main = ap.main, panel = function(x, y, subscripts,
                                                                       slopes, yadj, ...) {
                                        yadj <- yadj[[current.row()]]
                                        yadj <- yadj[current.column(), ]
                                        panel.xyplot(x, y - yadj)
                                        panel.abline(a = 0, b = slopes[current.column(),
                                                                       current.row()])
                                      }, slopes = score.mat, yadj = Yadd)
      }
    }
    ide.term <- (fa.inner.fun(faterms[[nt]]) == "ide")
    if (ide.term && nterms > 1) {
      ped.term <- character(0)
      for (ped in names(faterms)[-nt]) {
        if (fa.outer(faterms[[ped]]) == fa.outer(faterms[[nt]]) &&
            fa.inner.var(faterms[[ped]]) == fa.inner.var(faterms[[nt]]) &&
            fa.inner.fun(faterms[[ped]]) == "vm")
          ped.term <- match(ped, names(faterms))
      }
      if (length(ped.term) > 0) {
        fanam <- c(fanam, paste(fa.outer(faterms[[nt]]),
                                ":", fa.inner.var(faterms[[nt]]), "-total",
                                sep = ""))
        total.where <- total.where + 1
        Gmat <- gammas.lst[[nt]]$Gmat + gammas.lst[[ped.term]]$Gmat
        Cmat <- cov2cor(Gmat)
        gammas.lst[[total.where]] <- list(Gmat = Gmat,
                                          Cmat = Cmat)
        if (blups) {
          blup.inmet <- data.frame(Site = blup.lst[[nt]]$blups.inmet[[outer.name]],
                                   Variety = blup.lst[[nt]]$blups.inmet[[inner.name]],
                                   pres = blup.lst[[nt]]$blups.inmet$pres,
                                   blup = blup.lst[[nt]]$blups.inmet$blup +
                                     blup.lst[[ped.term]]$blups.inmet$blup,
                                   regblup = blup.lst[[nt]]$blups.inmet$regblup +
                                     blup.lst[[ped.term]]$blups.inmet$regblup,
                                   stringsAsFactors = TRUE)
          names(blup.inmet) <- c(outer.name, inner.name,
                                 "pres", "blup", "regblup")
          blup.lst[[total.where]] <- list(blups.inmet = blup.inmet)
        }
        if (heatmap) {
          dimnames(Cmat) <- list(sn, sn)
          if (heatmap.ord == "asis") {
            sn.ord <- sn
          }
          else if (heatmap.ord == "cluster") {
            dis.mat <- 1 - Cmat
            agnes.lst[[total.where]] <- agnes(x = dis.mat,
                                              diss = TRUE, method = agnes.method)
            sn.ord <- agnes.lst[[total.where]]$order.lab
          }
          else {
            sn.ord <- heatmap.ord
          }
          Cmat.ord <- Cmat[rev(sn.ord), rev(sn.ord)]
          diag(Cmat.ord) <- NA
          range(Cmat.ord, na.rm = T)
          atss <- seq(-1, 1, 0.1)
          hh <- rev(rainbow(256, start = 0, end = 2/3))
          hm.main <- fanam[length(fanam)]
          heatmap.lst[[total.where]] <- levelplot(Cmat.ord[sn.ord,
                                                           rev(sn.ord)], at = atss, col.regions = hh,
                                                  scales = list(x = list(rot = 60, cex = 0.9),
                                                                y = list(cex = 0.9)), xlab = hm.lab, ylab = hm.lab,
                                                  main = hm.main)
        }
      }
    }
  }
  if (uniplot) {
    names(uniplot.lst) <- names(faterms)
  }
  if (regplot) {
    names(regplot.lst) <- names(faterms)
  }
  if (addedplot) {
    names(addedplot.lst) <- names(faterms)
  }
  names(gammas.lst) <- fanam
  if (heatmap) {
    names(heatmap.lst) <- fanam
    if (heatmap.ord == "cluster") {
      names(agnes.lst) <- fanam
    }
  }
  if (blups)
    names(blup.lst) <- fanam
  return(list(gammas = gammas.lst, blups = blup.lst, uniplots = uniplot.lst,
              regplots = regplot.lst, addedplots = addedplot.lst,
              heatmaps = heatmap.lst, agnes = agnes.lst))
}


#' Variance explained
#'
#' @param model An asreml object with Factor analytic structure (fa2)
#'
#' @return vector
#' @export
#'
#' @examples
#' # in progress
var_fa <- function(model){
  ASM <- fa.asreml( model , trunc.char = NULL)
  L.star = ASM$gammas[[1]]$`rotated loads`
  psi = ASM$gammas[[1]]$`specific var`
  VarTot = sum(diag(L.star %*% t(L.star))) / sum(diag(L.star %*% t(L.star) + diag(psi) ))

  paf.site <- ASM$gammas[[1]]$`site %vaf`

  VarGenEnv <- diag(L.star %*% t(L.star) + diag(psi) )
  TotGenVar <- sum(VarGenEnv)

  VarFA1 <- sum(VarGenEnv*paf.site[,1])/100
  VarFA2 <- sum(VarGenEnv*paf.site[,2])/100

  PerVarFA1 <- round(VarFA1/TotGenVar*100,1)
  PerVarFA2 <- round(VarFA2/TotGenVar*100,1)

  return(perc = c(PerVarFA1,PerVarFA2))
}


"circleFun" <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


#' biplot fa2 model
#'
#' @param model An asreml object with Factor analytic structure (fa2 or more)
#' @param predictions dataframe with predicted values by site
#' @param one c(1,-1)
#' @param second c(1,-1)
#' @param fscore score to visualize genotypes and avoid overlapping
#' @param fmult in order to see the sites and genotypes in the same plot we have to play with this value
#' @param alpha opacity segments
#' @param alpha_site opacity sites
#' @param alpha_ind opacity individuals
#' @param subtitle text
#' @param gen vector with genotypes to show (optional)
#' @param size_ind_biplot size_ind_biplot
#'
#' @return list with 3 ggplots
#' @export
#'
#' @examples
#' # library(tidyverse)
#' # library(asreml)
#' # library(agridat)
#' # data(besag.met)
#' # dat <- besag.met
#' #
#' # dat <- dat %>% arrange(county)
#' # model <- asreml(fixed = yield ~ 1 + county,
#' #                 random = ~ fa(county, 2):gen + county:rep + diag(county):rep:block,
#' #                 residual = ~ dsum(~ units | county),
#' #                 data = dat,
#' #                 na.action = list(x="include",y="include"))
#' #
#' # pp <- predict(model, classify = "county")$pvals
#'
#' # sites
#' # p1 <- biplot_fa2(model,predictions = pp, fscore = 1, one = 1, second = 1, fmult = 5)$g1
#' # # sites scaled
#' # p2 <- biplot_fa2(model,predictions = pp, fscore = 1, one = 1, second = 1, fmult = 5)$g2
#' # # biplot
#' # p3 <- biplot_fa2(model,predictions = pp, fscore = 1, one = 1, second = 1, fmult = 5)$g3
biplot_fa2 <- function(model,
                       predictions,
                       one = -1,
                       second = -1,
                       fscore = 2,
                       fmult = 10,
                       alpha = 0.3,
                       alpha_site = 0.5,
                       alpha_ind = 0.5,
                       subtitle=NULL,
                       gen = NULL,
                       size_ind_biplot = 3){

  ASM <- fa.asreml( model , trunc.char = NULL)
  L.star = ASM$gammas[[1]]$`rotated loads`
  L.star[,1] <- L.star[,1]*one
  L.star[,2] <- L.star[,2]*second
  psi = ASM$gammas[[1]]$`specific var`
  Gvar <- ASM$gammas[[1]]$Gmat
  Cmat <- ASM$gammas[[1]]$Cmat
  Env.means <- predictions
  names(Env.means)[1:2] <- c("site", "BLUE")
  faComp <- data.frame(site = rownames(L.star), fa1 = L.star[,1], fa2 = L.star[,2], psi = psi, Vg = diag(Gvar), BLUE = Env.means$BLUE)

  percentg <- var_fa(model)

  # Without Standardize Loadings
  d=data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=L.star[,1], vy=L.star[,2])
  loadings = ggplot2::ggplot(faComp, ggplot2::aes(x = fa1, y = fa2)) +
    ggplot2::geom_point(ggplot2::aes(colour = Vg, size = BLUE)) +
    ggplot2::scale_colour_gradient(low = "pink", high = "blue") +
    ggrepel::geom_label_repel(ggplot2::aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1, alpha = alpha_site) +
    ggplot2::ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    ggplot2::xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ggplot2::ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    ggplot2::theme_bw(base_size = 15)+
    ggplot2::geom_vline(xintercept = 0,linetype = 2) +
    ggplot2::geom_hline(yintercept = 0,linetype = 2)+
    ggplot2::geom_segment(data=d,
                 mapping=ggplot2::aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha= alpha)

  # Standardize Loadings
  faCompR <- faComp
  faCompR[,2:3] <- diag(1/sqrt(diag(Gvar))) %*% L.star
  d <- data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=faCompR[,2], vy=faCompR[,3])
  circle <- circleFun(c(0,0),2,npoints = 100)
  loading_C <- ggplot2::ggplot(faCompR, ggplot2::aes(x = fa1, y = fa2)) +
    ggplot2::geom_point(ggplot2::aes(colour = Vg, size = BLUE)) +
    ggplot2::scale_colour_gradient(low = "pink", high = "blue") +
    ggrepel::geom_label_repel(ggplot2::aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1, alpha = alpha_site) +
    ggplot2::ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    ggplot2::xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ggplot2::ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    ggplot2::theme_bw(base_size = 15)+
    ggplot2::geom_vline(xintercept = 0,linetype = 2) +
    ggplot2::geom_hline(yintercept = 0,linetype = 2)+
    ggplot2::geom_segment(data=d,
                 mapping=ggplot2::aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha=alpha) +
    ggplot2::geom_path(data=circle, ggplot2::aes(x,y))


  # Biplot
  fa12_scores = ASM$blups[[1]]$scores
  names(fa12_scores)[2:3] <- c("comp", "Genotype")
  fa12_scores %<>% select(-blup) %>%  spread(. ,"comp", value = "blupr")
  names(fa12_scores) = c("Genotype", "fa1", "fa2")
  fa12_scores$fa1 <- fa12_scores$fa1*one
  fa12_scores$fa2 <- fa12_scores$fa2*second
  message("summary fscore")
  print(summary(sqrt(fa12_scores$fa1^2+fa12_scores$fa2^2)))
  fa12_scores$Score <-  ifelse(sqrt(fa12_scores$fa1^2+fa12_scores$fa2^2)>fscore,1,0)

  if(!is.null(gen)){

    fa12_scores[fa12_scores$Genotype %in% gen, "Score" ] <- 1

  }

  d=data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=L.star[,1], vy=L.star[,2])
  biplot = ggplot2::ggplot(faComp, ggplot2::aes(x = fa1, y = fa2)) +
    ggplot2::geom_point(ggplot2::aes(colour = Vg, size = BLUE)) +
    ggplot2::scale_colour_gradient(low = "pink", high = "blue") +
    ggrepel::geom_label_repel(ggplot2::aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1,  alpha = alpha_site) +
    ggplot2::ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    ggplot2::xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ggplot2::ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    ggplot2::theme_bw(base_size = 15)+
    ggplot2::geom_vline(xintercept = 0,linetype = 2) +
    ggplot2::geom_hline(yintercept = 0,linetype = 2)+
    ggplot2::geom_segment(data=d,
                 mapping=ggplot2::aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha= alpha) +
    ggrepel::geom_label_repel(data = subset(fa12_scores, Score==1),
                              ggplot2::aes(label = Genotype, x = fmult*fa1 , y= fmult*fa2),
                     colour = "red",segment.colour = "red" , size=size_ind_biplot, alpha = alpha_ind)


  return(list(g1=loadings, g2= loading_C, g3 = biplot))
}


#' Correlation Covariance heatmap
#'
#' @param matrix matrix
#' @param corr logical (TRUE default, correlation matrix)
#' @param size letter size
#'
#' @return ggplot
#' @export
#'
#' @examples
#' # data(iris)
#' # M = cor(iris[,-5])
#' # covcor_heat(M, corr = T)
covcor_heat <- function(matrix, corr = TRUE, size = 4){

  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }

  cormat <- reorder_cormat(matrix)
  upper_tri <- get_upper_tri(matrix)
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

  u <- -1
  m <- 0
  l <- 1
  main <- "Correlation"
  col_pallete <- c("#db4437","white","#4285f4")
  col_letter <- "black"

  if(isFALSE(corr)){
    u <- min(matrix, na.rm = T)
    l <- max(matrix, na.rm = T)
    m <- u+(l-u)/2
    main <- "Covariance"
    col_pallete <- c("#440154","#21908C","#FDE725")
    col_letter <- "white"
  }

  melted_cormat$Var1 <- as.factor(melted_cormat$Var1)
  melted_cormat$Var2 <- as.factor(melted_cormat$Var2)

  ggheatmap <-
    ggplot2::ggplot(melted_cormat, ggplot2::aes(Var2, Var1, fill = value))+
    ggplot2::geom_tile(color = "white")+
    ggplot2::scale_fill_gradient2(low = col_pallete[1], high = col_pallete[3], mid = col_pallete[2],  # color= c("#440154","#21908C","#FDE725")
                         midpoint = m, limit = c(u,l), space = "Lab",
                         name=main) +
    ggplot2::theme_minimal()+ # minimal theme
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1),
          axis.text.y = ggplot2::element_text(size = 12))
  # coord_fixed()


  plot <- ggheatmap +
    ggplot2::geom_text(ggplot2::aes(Var2, Var1, label = value), color = col_letter, size = size) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))

  return(plot)

}

