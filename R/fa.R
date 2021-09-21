

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
        uniplot.lst[[nt]] <- lattice::xyplot(Lamc[, 1] ~ Lamc[,
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
            tmp[[(i - 1) * (k - 1) + (j - i)]] <- lattice::xyplot(Lamc[,
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
      blup.df <- data.frame(blup = as.vector(cc))
      nn <- dimnames(cc)[[1]]
      temp <- strsplit(nn, split = ":", fixed = TRUE)
      tt.inner <- sapply(temp, function(x) x[2])
      tt.outer <- sapply(temp, function(x) x[1])
      tt.outer <- sapply(strsplit(tt.outer, split = "_",
                                  fixed = T), function(x) paste(x[-1], collapse = "_"))
      tt.inner <- sapply(strsplit(tt.inner, split = "_",
                                  fixed = T), function(x) paste(x[-1], collapse = "_"))
      # tt.inner <- sapply(strsplit(tt.inner, split = "_",
      #                             fixed = T), function(x) paste(x[-c(1,2)], collapse = "_"))
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
      regplot.lst[[nt]] <- lattice::xyplot(xform, data = regplot.df,
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
        addedplot.lst[[nt]] <- lattice::xyplot(xform, data = regplot.df,
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
                                     blup.lst[[ped.term]]$blups.inmet$regblup)
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
  loadings = ggplot(faComp, aes(x = fa1, y = fa2)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_label_repel(aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1, alpha = alpha_site) +
    ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    theme_bw(base_size = 15)+
    geom_vline(xintercept = 0,linetype = 2) + geom_hline(yintercept = 0,linetype = 2)+
    geom_segment(data=d,
                 mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha= alpha)

  # Standardize Loadings
  faCompR <- faComp
  faCompR[,2:3] <- diag(1/sqrt(diag(Gvar))) %*% L.star
  d <- data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=faCompR[,2], vy=faCompR[,3])
  circle <- circleFun(c(0,0),2,npoints = 100)
  loading_C <- ggplot(faCompR, aes(x = fa1, y = fa2)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_label_repel(aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1, alpha = alpha_site) +
    ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    theme_bw(base_size = 15)+
    geom_vline(xintercept = 0,linetype = 2) + geom_hline(yintercept = 0,linetype = 2)+
    geom_segment(data=d,
                 mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha=alpha) +
    geom_path(data=circle, aes(x,y))


  # Biplot
  fa12_scores = ASM$blups[[1]]$scores
  names(fa12_scores)[2:3] <- c("comp", "Genotype")
  fa12_scores %<>% select(-blup) %>%  spread(. ,"comp", value = "blupr")
  names(fa12_scores) = c("Genotype", "fa1", "fa2")
  fa12_scores$fa1 <- fa12_scores$fa1*one
  fa12_scores$fa2 <- fa12_scores$fa2*second
  fa12_scores$Score <-  ifelse(sqrt(fa12_scores$fa1^2+fa12_scores$fa2^2)>fscore,1,0)

  if(!is.null(gen)){

    fa12_scores[fa12_scores$Genotype %in% gen, "Score" ] <- 1

  }

  d=data.frame(x=rep(0, nrow(L.star)), y=rep(0, nrow(L.star)), vx=L.star[,1], vy=L.star[,2])
  biplot = ggplot(faComp, aes(x = fa1, y = fa2)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_label_repel(aes(label = site), nudge_y= 0.05, nudge_x=-0.03, force=1,  alpha = alpha_site) +
    ggtitle(paste0("Environment Factor Loadings ", "(",sum(percentg),"%)"), subtitle = subtitle) +
    xlab(paste0("FA1 loading ", "(",percentg[1],"%)" )) +
    ylab(paste0("FA2 loading ", "(",percentg[2],"%)")) +
    theme_bw(base_size = 15)+
    geom_vline(xintercept = 0,linetype = 2) + geom_hline(yintercept = 0,linetype = 2)+
    geom_segment(data=d,
                 mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy),
                 arrow=arrow(), size=0.5, color="black", alpha= alpha) +
    geom_label_repel(data = subset(fa12_scores, Score==1),
                     aes(label = Genotype, x = fmult*fa1 , y= fmult*fa2),
                     colour = "red",segment.colour = "red" , size=size_ind_biplot, alpha = alpha_ind)


  return(list(g1=loadings, g2= loading_C, g3 = biplot))
}
