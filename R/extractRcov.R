#' Extract Residual Variance-Covariance from ASReml-R
#'
#' This function is specially useful for extracting residual variance covariance
#' matrices from ASReml-R when running repeated measurements analysis.
#'
#' @param model An asreml object.
#' @param time An optional character string indicating the "Time". By default the function
#' identifies this parameter.
#' @param plot An optional character string indicating the "PlotID". By default the function
#' identifies this parameter.
#' @param vc_error An optional character string indicating the variance covariance.
#'  It can be "corv", "corh", "corgh", "us", "expv", "exph", "ar1v", "ar1h" or "ante".
#'  By using NULL the function tries to guess which was the variance-covariance used.
#'
#' @return An object with a list of:
#' \item{corr_mat}{A matrix with the residual correlation between time points.}
#' \item{vcov_mat}{A matrix of the estimated residual variance-covariance between time points.}
#' \item{vc}{A character string indicating the variance-covariance fitted.}
#'
#' @details The expected residual variance covariance structure must be of the form:
#'  `~id(Plot):corv(Time)`, where `Plot` is a unique identifier of each experimental unit,
#'  and `Time`, represents the variable that contains the time when
#'  the experimental units were measured. This form also requires that the levels of
#'  the factor Time are nested in the levels of the factor Plot. If it is not in that form
#'  you can sort the dataset by using the following  command `arrange(grassUV, Plant, Time)`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggpubr)
#' library(agriutilities)
#' library(tidyverse)
#' library(asreml)
#'
#' head(grassUV)
#' str(grassUV)
#'
#' # Exploration -------------------------------------------------------------
#'
#' grassUV %>%
#'   ggplot(
#'     aes(x = Time, y = y, group = Plant, color = Plant)
#'   ) +
#'   geom_point() +
#'   geom_line() +
#'   facet_wrap(~Tmt) +
#'   theme_minimal(base_size = 15)
#'
#' tmp <- grassUV %>%
#'   group_by(Time, Plant) %>%
#'   summarise(mean = mean(y, na.rm = TRUE)) %>%
#'   spread(Time, mean) %>%
#'   column_to_rownames("Plant")
#'
#' gg_cor(tmp, label_size = 5)
#'
#' tmp %>%
#'   cor(use = "pairwise.complete.obs") %>%
#'   as.data.frame() %>%
#'   rownames_to_column(var = "Time") %>%
#'   gather("DAP2", "corr", -1) %>%
#'   type.convert(as.is = FALSE) %>%
#'   mutate(corr = ifelse(Time < DAP2, NA, corr)) %>%
#'   mutate(DAP2 = as.factor(DAP2)) %>%
#'   ggplot(
#'     aes(x = Time, y = corr, group = DAP2, color = DAP2)
#'   ) +
#'   geom_point() +
#'   geom_line() +
#'   theme_minimal(base_size = 15) +
#'   color_palette(palette = "jco") +
#'   labs(color = "Time", y = "Pearson Correlation")
#'
#' # Modeling ----------------------------------------------------------------
#'
#' # Identity variance model.
#' model_0 <- asreml(
#'   fixed = y ~ Time + Tmt + Tmt:Time,
#'   residual = ~ id(Plant):idv(Time),
#'   data = grassUV
#' )
#'
#' # Simple correlation model; homogeneous variance form.
#' model_1 <- asreml(
#'   fixed = y ~ Time + Tmt + Tmt:Time,
#'   residual = ~ id(Plant):corv(Time),
#'   data = grassUV
#' )
#'
#' # Exponential (or power) model; homogeneous variance form.
#' model_2 <- asreml(
#'   fixed = y ~ Time + Tmt + Tmt:Time,
#'   residual = ~ id(Plant):expv(Time),
#'   data = grassUV
#' )
#'
#' # Exponential (or power) model; heterogeneous variance form.
#' model_3 <- asreml(
#'   fixed = y ~ Time + Tmt + Tmt:Time,
#'   residual = ~ id(Plant):exph(Time),
#'   data = grassUV
#' )
#'
#' # Antedependence variance model of order 1
#' model_4 <- asreml(
#'   fixed = y ~ Time + Tmt + Tmt:Time,
#'   residual = ~ id(Plant):ante(Time),
#'   data = grassUV
#' )
#'
#' # Autoregressive model of order 1; homogeneous variance form.
#' model_5 <- asreml(
#'   fixed = y ~ Time + Tmt + Tmt:Time,
#'   residual = ~ id(Plant):ar1v(Time),
#'   data = grassUV
#' )
#'
#' # Autoregressive model of order 1; heterogeneous variance form.
#' model_6 <- asreml(
#'   fixed = y ~ Time + Tmt + Tmt:Time,
#'   residual = ~ id(Plant):ar1h(Time),
#'   data = grassUV
#' )
#'
#' # Unstructured variance model.
#' model_7 <- asreml(
#'   fixed = y ~ Time + Tmt + Tmt:Time,
#'   residual = ~ id(Plant):us(Time),
#'   data = grassUV
#' )
#'
#' # Model Comparison --------------------------------------------------------
#'
#' models <- list(
#'   "id" = model_0,
#'   "cor" = model_1,
#'   "exp" = model_2,
#'   "exph" = model_3,
#'   "ante" = model_4,
#'   "ar1" = model_5,
#'   "ar1h" = model_6,
#'   "us" = model_7
#' )
#'
#' summary_models <- data.frame(
#'   model = names(models),
#'   aic = unlist(lapply(models, \(x) summary(x)$aic)),
#'   bic = unlist(lapply(models, \(x) summary(x)$bic)),
#'   loglik = unlist(lapply(models, \(x) summary(x)$loglik)),
#'   nedf = unlist(lapply(models, \(x) summary(x)$nedf)),
#'   param = unlist(lapply(models, \(x) attr(summary(x)$aic, "param"))),
#'   row.names = NULL
#' )
#'
#' summary_models %>%
#'   ggplot(
#'     aes(x = reorder(model, -bic), y = bic, group = 1)
#'   ) +
#'   geom_point(size = 2) +
#'   geom_text(aes(x = model, y = bic + 5, label = param)) +
#'   geom_line() +
#'   theme_minimal(base_size = 15) +
#'   labs(x = NULL)
#'
#' # Extracting Variance Covariance Matrix -----------------------------------
#'
#' extract_rcov(model_4)
#'
#' covcor_heat(
#'   matrix = extract_rcov(model_1)$corr,
#'   legend = "none",
#'   size = 5
#' ) + ggtitle(label = "Uniform Correlation (corv)")
#' covcor_heat(
#'   matrix = extract_rcov(model_2)$corr,
#'   legend = "none",
#'   size = 5
#' ) + ggtitle(label = "Exponetial (expv)")
#' }
extract_rcov <- function(model = NULL,
                         time = NULL,
                         plot = NULL,
                         vc_error = NULL) {
  stopifnot(inherits(x = model, what = "asreml"))
  options <- c(
    "corv", "corh", "corgh",
    "us", "expv", "exph",
    "ar1v", "ar1h", "ante"
  )
  str_res <- as.character(model$call$residual)[2]
  vc_check <- as.character(model$call$residual[[2]][[3]][[1]])
  if (is.null(vc_error)) {
    vc_error <- vc_check
    if (!vc_error %in% options) {
      stop(
        "Variance covariance Not Available: \n\n\t",
        "vc_error = '", vc_error, "'\n\t",
        "residual = ~", str_res
      )
    }
  } else {
    stopifnot(vc_error %in% options)
    if (vc_error != vc_check) {
      stop(
        "Check that the argument vc_error matches your structure: \n\n\t",
        "vc_error = '", vc_error, "'\n\t",
        "residual = ~", str_res
      )
    }
  }
  if (is.null(plot) || is.null(time)) {
    plot <- as.character(model$call$residual[[2]][[2]][[2]])
    time <- as.character(model$call$residual[[2]][[3]][[2]])
  }
  lvls <- levels(data.frame(model$mf)[, time])
  s <- length(lvls)
  corr <- matrix(1, ncol = s, nrow = s, dimnames = list(lvls, lvls))
  vc <- summary(model)$varcomp
  pt <- paste0(plot, ":", time, "!", time)
  vc <- vc[grep(pt, rownames(vc), fixed = TRUE), ]
  if (vc_error == "corv") {
    corr <- vc[1, 1] * corr
    diag(corr) <- rep(1, s)
    D <- diag(rep(vc[2, 1], s))
    colnames(D) <- rownames(D) <- lvls
    vcov <- sqrt(D) %*% corr %*% sqrt(D)
    objt <- list(corr_mat = corr, vcov_mat = vcov, vc = vc_error)
  }
  if (vc_error == "corh") {
    corr <- vc[1, 1] * corr
    diag(corr) <- rep(1, s)
    D <- diag(vc[2:(s + 1), 1])
    colnames(D) <- rownames(D) <- lvls
    vcov <- sqrt(D) %*% corr %*% sqrt(D)
    objt <- list(corr_mat = corr, vcov_mat = vcov, vc = vc_error)
  }
  if (vc_error == "corgh") {
    vc_corr <- vc[grep(".cor", rownames(vc)), ]
    vc_var <- vc[-grep(".cor", rownames(vc)), ]
    k <- 1
    for (i in 1:s) {
      for (j in 1:i) {
        if (i != j) {
          corr[i, j] <- vc_corr[k, 1]
          corr[j, i] <- vc_corr[k, 1]
          k <- k + 1
        }
      }
    }
    D <- diag(vc_var[1:s, 1])
    colnames(D) <- rownames(D) <- lvls
    vcov <- sqrt(D) %*% corr %*% sqrt(D)
    objt <- list(corr_mat = corr, vcov_mat = vcov, vc = vc_error)
  }
  if (vc_error == "us") {
    vcov <- matrix(0, ncol = s, nrow = s, dimnames = list(lvls, lvls))
    k <- 1
    for (i in 1:s) {
      for (j in 1:i) {
        vcov[i, j] <- vc[k, 1]
        k <- k + 1
      }
    }
    vcov[upper.tri(vcov)] <- t(vcov)[upper.tri(vcov)]
    corr <- cov2cor(vcov)
    objt <- list(corr_mat = corr, vcov_mat = vcov, vc = vc_error)
  }
  if (vc_error == "expv") {
    x <- as.numeric(lvls)
    D <- diag(vc[2, 1], nrow = s, ncol = s)
    colnames(D) <- rownames(D) <- lvls
    for (i in 1:s) {
      for (j in 1:i) {
        corr[i, j] <- (vc[1, 1]^(abs(x[i] - x[j])))
      }
    }
    corr[upper.tri(corr)] <- t(corr)[upper.tri(corr)]
    vcov <- sqrt(D) %*% corr %*% sqrt(D)
    objt <- list(corr_mat = corr, vcov_mat = vcov, vc = vc_error)
  }
  if (vc_error == "exph") {
    x <- as.numeric(lvls)
    D <- diag(vc[2:(1 + s), 1], nrow = s, ncol = s)
    colnames(D) <- row.names(D) <- lvls
    for (i in 1:s) {
      for (j in 1:i) {
        corr[i, j] <- (vc[1, 1]^(abs(x[i] - x[j])))
      }
    }
    corr[upper.tri(corr)] <- t(corr)[upper.tri(corr)]
    vcov <- sqrt(D) %*% corr %*% sqrt(D)
    objt <- list(corr_mat = corr, vcov_mat = vcov, vc = vc_error)
  }
  if (vc_error == "ar1v") {
    D <- diag(rep(vc[2, 1], s), nrow = s, ncol = s)
    colnames(D) <- row.names(D) <- lvls
    for (i in 1:s) {
      for (j in 1:i) {
        corr[i, j] <- (vc[1, 1]^(i - j))
      }
    }
    corr[upper.tri(corr)] <- t(corr)[upper.tri(corr)]
    vcov <- sqrt(D) %*% corr %*% sqrt(D)
    objt <- list(corr_mat = corr, vcov_mat = vcov, vc = vc_error)
  }
  if (vc_error == "ar1h") {
    D <- diag(vc[2:(1 + s), 1], nrow = s, ncol = s)
    colnames(D) <- row.names(D) <- lvls
    for (i in 1:s) {
      for (j in 1:i) {
        corr[i, j] <- (vc[1, 1]^(abs(i - j)))
      }
    }
    corr[upper.tri(corr)] <- t(corr)[upper.tri(corr)]
    vcov <- sqrt(D) %*% corr %*% sqrt(D)
    objt <- list(corr_mat = corr, vcov_mat = vcov, vc = vc_error)
  }
  if (vc_error == "ante") {
    U <- diag(1, nrow = s, ncol = s)
    D <- diag(vc[seq(1, nrow(vc), 2), 1])
    colnames(D) <- row.names(D) <- lvls
    colnames(U) <- row.names(U) <- lvls
    U[row(U) + 1 == col(U)] <- vc[seq(2, nrow(vc), 2), 1]
    vcov <- solve(U %*% D %*% t(U))
    corr <- cov2cor(vcov)
    objt <- list(corr_mat = corr, vcov_mat = vcov, vc = vc_error, U = U, D = D)
  }
  return(objt)
}
