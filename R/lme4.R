#' ran
#'
#' @param var factor to be taken as random
#'
#' @return string
#' @noRd
#'
#' @examples
#' # ran("rep")
ran <- function(var) {
  effect <- paste0("(", 1, "|", var, ")")
  return(effect)
}

#' varG
#'
#' @param model lmer model
#' @param comp String random component to be extracted
#'
#' @return variance component
#' @noRd
#'
#' @examples
#' # Varg(model, "gen")
VarG <- function(model, comp) {
  v <- as.data.frame(lme4::VarCorr(model))
  v <- v[v$grp == comp, "vcov"]
  return(v)
}

#' VarE
#'
#' @param model lmer model
#'
#' @return residual variance
#' @noRd
#'
#' @examples
#' # VarE(model)
VarE <- function(model) {
  v <- as.data.frame(lme4::VarCorr(model))
  v <- v[v$grp == "Residual", "vcov"]
  return(v)
}

#' Cullis heritability for lme4 models
#'
#' @param model Object of class \code{lmer}.
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param re_MME A logical value to ask if we want to reconstruct the mixed models
#' equations to estimate the Cullis heritability. (\code{FALSE} by default)
#'
#' @author Paul Schmidt, Johan Aparicio.
#'
#' @return A numerical value of the Cullis heritability estimate. If
#' \code{re_MME} is \code{TRUE}, a list with matrices of the mixed models
#' equations is returned.
#' @export
#'
#' @examples
#' \donttest{
#' library(lme4)
#' library(agridat)
#' library(agriutilities)
#' dat <- john.alpha
#' g.ran <- lmer(
#'   formula = yield ~ rep + (1 | gen) + (1 | rep:block),
#'   data = dat
#' )
#' h_cullis(model = g.ran, genotype = "gen")
#' }
h_cullis <- function(model, genotype, re_MME = FALSE) {
  if (re_MME) {
    Gen_levels <- levels(model@frame[, genotype])

    vc <- as.data.frame(lme4::VarCorr(model))
    # Number of random effects
    n.ran <- nrow(vc)
    # R = varcov-matrix for error term
    n <- length(summary(model)$residuals) # numer of observations
    vc.e <- vc[vc$grp == "Residual", "vcov"] # error vc
    R <- diag(n) * vc.e # R matrix = I * vc.e
    # names of random effects
    nomb <- names(summary(model)$ngrps)
    # Genotype
    n.g <- summary(model)$ngrps[which(nomb == genotype)]
    vc.g <- vc[vc$grp == genotype, "vcov"]

    # G matrix of random effects
    G.tmp <- list()
    if (n.ran >= 2) {
      for (i in 1:(n.ran - 1)) { # remove the residual variance
        n.tmp <- summary(model)$ngrps[which(nomb == nomb[i])]
        vc.tmp <- vc[vc$grp == nomb[i], "vcov"]
        G.tmp[[i]] <- diag(n.tmp) * vc.tmp
      }
    }

    G <- Matrix::bdiag(G.tmp) # G is blockdiagonal with G.g and G.b

    # Design Matrices
    X <- as.matrix(lme4::getME(model, "X")) # Design matrix fixed effects
    Z <- as.matrix(lme4::getME(model, "Z")) # Design matrix random effects

    # Mixed model Equation
    C11 <- t(X) %*% solve(R) %*% X
    C12 <- t(X) %*% solve(R) %*% Z
    C21 <- t(Z) %*% solve(R) %*% X
    C22 <- t(Z) %*% solve(R) %*% Z + solve(G)

    C <- as.matrix(rbind(
      cbind(C11, C12), # Combine components into one matrix C
      cbind(C21, C22)
    ))

    # Mixed model Equation Solutions
    C.inv <- solve(C) # Inverse of C
    C22.g <- C.inv[Gen_levels, Gen_levels] # subset of C.inv that refers to genotypic BLUPs

    # Mean variance of BLUP-difference from C22 matrix of genotypic BLUPs
    vdBLUP.sum <- n.g * sum(diag(C22.g)) - sum(C22.g)
    vdBLUP.avg <- vdBLUP.sum * (2 / (n.g * (n.g - 1))) # mean variance of BLUP-difference = divide sum by number of genotype pairs

    H2Cullis <- 1 - (vdBLUP.avg / 2 / vc.g)

    return(
      list(
        X = X,
        Z = Z,
        G = G,
        R = R,
        C11 = C11,
        C12 = C12,
        C21 = C21,
        C22 = C22,
        C = C,
        C.inv = C.inv,
        C22.g = C22.g,
        vdBLUP_avg = vdBLUP.avg,
        H2Cullis = H2Cullis
      )
    )
  } else {
    aveped <- mean(attr(lme4::ranef(model, drop = TRUE)[[genotype]], "postVar"))
    vc.g <- VarG(model, genotype)
    H2Cullis <- ifelse(vc.g == 0, 0, 1 - aveped / vc.g)
    return(H2Cullis)
  }
}

#' residuals lme4
#'
#' @param model lmer model
#' @param returnN logical value for returning the number of extreme observations
#' @param k  number of standard desviations to define values as extremes
#'
#' @return data.frame of number of extreme observations depending on returnN
#' @noRd
#'
#' @examples
#' # in progress
res_lme4 <- function(model, returnN = FALSE, k = 3) {
  res <- stats::residuals(model, scaled = TRUE)
  data <- model@frame
  data$residual <- res
  data$outlier <- NA
  data$outlier[which(abs(data$res) >= k)] <- TRUE
  data$outlier[which(abs(data$res) < k)] <- FALSE
  ix <- ifelse(length(which(abs(res) > k)) >= 1, length(which(abs(res) > k)), 0)
  if (!returnN) {
    return(data)
  } else {
    return(ix)
  }
}

#' @noRd
#' @keywords internal
mult_models <- function(data = NULL,
                        equation = NULL,
                        by = NULL,
                        mixed_model = TRUE,
                        progress = TRUE) {
  models <- list()
  data[[by]] <- as.factor(data[[by]])
  lvls <- levels(data[, by])
  for (exp in lvls) {
    tmp_dt <- dplyr::filter(data, .data[[by]] %in% exp)

    if (progress) {
      trait <- all.vars(equation)[1]
      cat(paste0(
        "Fitting models for ", paste(trait, collapse = ", "),
        " in ", exp, ".\n"
      ))
    }

    if (mixed_model) {
      model <- try(
        expr = suppressMessages(
          lmerTest::lmer(equation, data = tmp_dt, na.action = na.omit)
        ),
        silent = TRUE
      )
    } else {
      model <- try(
        expr = suppressMessages(
          stats::lm(equation, data = tmp_dt, na.action = na.omit)
        ),
        silent = TRUE
      )
    }
    if (inherits(model, "try-error")) {
      models[[exp]] <- NULL
    } else {
      models[[exp]] <- model
    }
  }
  return(models)
}

#' @noRd
#' @keywords internal
mult_summary <- function(models_fixed = NULL,
                         models_rand = NULL,
                         genotype = "Name",
                         exp_design = "unknown") {
  trials <- names(models_fixed)
  gv <- unlist(lapply(models_rand, VarG, genotype))
  ev <- unlist(lapply(models_rand, VarE))
  he <- unlist(lapply(models_rand, h_cullis, genotype))
  cv <- base::sapply(X = models_fixed, function(mf0) {
    100 * summary(mf0)$sigma / abs(mean(stats::fitted(mf0), na.rm = TRUE))
  })
  summ <- data.frame(
    trial = trials,
    heritability = he,
    CV = cv,
    VarGen = gv,
    VarErr = ev,
    design = exp_design,
    row.names = NULL
  )
  return(summ)
}

#' @noRd
#' @keywords internal
get_blup_blues <- function(models_fixed = NULL,
                           models_rand = NULL,
                           genotype = "gen") {
  blues <- sapply(X = models_fixed, FUN = function(mf0) {
    emmeans::emmeans(mf0, specs = genotype) %>%
      as.data.frame() %>%
      dplyr::select(1:3)
  }, simplify = FALSE)
  blues <- dplyr::bind_rows(blues, .id = "trial")
  names(blues) <- c("trial", "genotype", "BLUEs", "seBLUEs")

  blup_func <- function(model, genotype) {
    coefs <- coef(model)[[genotype]]
    predvals <- data.frame(rownames(coefs), coefs[, "(Intercept)"])
    names(predvals) <- c("genotype", "BLUPs")
    ran_effs <- lme4::ranef(model, condVar = TRUE)[[genotype]]
    pred_err <- data.frame(
      rownames(ran_effs),
      as.vector(sqrt(attr(ran_effs, "postVar")))
    )
    colnames(pred_err) <- c("genotype", "seBLUPs")
    predvals <- merge(predvals, pred_err, by = "genotype", all = TRUE)
    return(predvals)
  }
  blups <- sapply(
    X = models_rand,
    FUN = blup_func,
    genotype = genotype,
    simplify = FALSE
  ) %>%
    dplyr::bind_rows(., .id = "trial")

  blues_blups <- merge(
    x = blues,
    y = blups,
    by = c("genotype", "trial"),
    all = TRUE
  ) %>%
    dplyr::mutate(wt = (1 / seBLUEs)^2)
  return(blues_blups)
}

#' @noRd
#' @keywords internal
get_residuals <- function(models, k = 3) {
  data_out <- lapply(models, res_lme4, FALSE, k = k) %>%
    dplyr::bind_rows(., .id = "trial")
  outliers <- data_out[data_out$outlier %in% TRUE, ]
  return(list(residuals = data_out, outliers = outliers))
}


#' Fit multiples CRD models
#'
#' @param data A data.frame in a wide format.
#' @param trial A character string indicating the column in data that contains
#' trials.
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param response A character vector specifying the traits for which the models
#' should be fitted.
#'
#' @return a list of data.frames.
#' @noRd
#'
#' @examples
#' # library(agridat)
#' # library(agriutilities)
#' # data(besag.met)
#' # dat <- besag.met
#' # res <- fit_crd(
#' #   data = dat,
#' #   trial = "county",
#' #   genotype = "gen",
#' #   response = "yield"
#' # )
fit_crd <- function(data = NULL,
                    trial = NULL,
                    genotype = NULL,
                    response = NULL,
                    progress = TRUE) {
  data <- as.data.frame(data)
  # lme4 equation
  equation_ran <- stats::reformulate(ran(genotype), response = response)
  # lm equation
  equation_fixed <- stats::reformulate(genotype, response = response)
  # fitting models
  models_rand <- mult_models(
    data = data,
    equation = equation_ran,
    by = trial,
    mixed_model = TRUE,
    progress = progress
  )
  models_fixed <- mult_models(
    data = data,
    equation = equation_fixed,
    by = trial,
    mixed_model = FALSE,
    progress = FALSE
  )
  # summary
  mt_summ <- mult_summary(
    models_fixed = models_fixed,
    models_rand = models_rand,
    genotype = genotype,
    exp_design = "CRD"
  )
  # BLUEs and BLUPs
  blues_blups <- get_blup_blues(
    models_fixed = models_fixed,
    models_rand = models_rand,
    genotype = genotype
  )
  # Residuals
  residuals <- get_residuals(models = models_rand)
  outliers <- merge(
    x = data,
    y = residuals$outliers,
    by.x = c(trial, genotype, response),
    by.y = c("trial", genotype, response),
    all.y = TRUE
  ) %>%
    dplyr::select(dplyr::all_of(c(trial, genotype)), id, outlier)
  names(outliers) <- c("trial", "genotype", "id", "outlier")

  residuals <- residuals$residuals %>%
    dplyr::select(trial, all_of(genotype), residual)
  names(residuals) <- c("trial", "genotype", response)
  # results
  output <- list(
    models_rand = models_rand,
    models_fixed = models_fixed,
    resum_fitted_model = mt_summ,
    blues_blups = blues_blups,
    residuals = residuals,
    outliers = outliers
  )
  return(output)
}


#' @noRd
res_qqplot <- function(data_out, title = NULL) {
  q <- dplyr::filter(data_out, !is.na(Classify)) %>%
    ggpubr::ggqqplot(
      x = "Residuals",
      fill = "Classify",
      ylab = "Sample Quantile",
      xlab = "Theoretical Quantile", title = title
    )
  return(q)
}
