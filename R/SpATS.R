#' pseudo r square SpATS
#'
#' @param model SpATS object
#'
#' @return value
#' @noRd
#'
#' @examples
#' # in progress
r_square <- function(model) {
  response <- model$data[, model$model$response]
  mean.response <- mean(response, na.rm = TRUE)
  fitted <- model$fitted
  SS.fitted <- sum((response - fitted)^2, na.rm = TRUE)
  SS.response <- sum((response - mean.response)^2, na.rm = TRUE)
  R <- 1 - SS.fitted / SS.response
  names(R) <- "r.square"
  return(round(R, 3))
}

#' coefficient of variation SpATS
#'
#' @param model SpATS object
#'
#' @return value
#' @noRd
#'
#' @examples
#' # in progress
CV_SpATS <- function(model) {
  response <- model$data[, model$model$response]
  mean.response <- mean(response, na.rm = TRUE)
  fitted <- model$fitted
  MSE <- mean((response - fitted)^2, na.rm = TRUE)
  RMSE <- sqrt(MSE)
  NRMSE <- RMSE / mean.response
  cv_prcnt <- NRMSE * 100
  names(cv_prcnt) <- "CV"
  return(round(cv_prcnt, 2))
}


#' residuals SpATS
#'
#' @param model SpATS object
#' @param k number of standard desviations to define values as extremes (default = 3)
#'
#' @return data.frame
#' @noRd
#'
#' @examples
#' # in progress
res_spats <- function(model, k = 3) {
  dt <- model$data
  VarE <- model$psi[1]
  Data <- data.frame(
    Index = seq_along(stats::residuals(model)),
    Residuals = stats::residuals(model)
  )
  u <- +k * sqrt(VarE)
  l <- -k * sqrt(VarE)
  Data$Classify <- NA
  Data$Classify[which(abs(Data$Residuals) >= u)] <- "Outlier"
  Data$Classify[which(abs(Data$Residuals) < u)] <- "Normal"
  Data$l <- l
  Data$u <- u
  Data$gen <- dt[, model$model$geno$genotype]
  Data$col <- dt[, model$terms$spatial$terms.formula$x.coord]
  Data$row <- dt[, model$terms$spatial$terms.formula$y.coord]
  Data$fit <- stats::fitted.values(model)
  Data$response <- dt[, model$model$response]
  return(Data)
}

#' Interactive residuals plot SpATS
#'
#' @param data_out res_spats data.frame
#'
#' @return ggplot
#' @noRd
#'
#' @examples
#' # in progress
plot_res_index <- function(data_out) {
  k <- dplyr::filter(data_out, !is.na(Classify)) %>%
    ggplot2::ggplot(ggplot2::aes(x = Index, y = Residuals, color = Classify)) +
    ggplot2::geom_point(size = 2, alpha = 0.5, na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c("grey80", "red")) +
    ggplot2::geom_hline(yintercept = data_out$u, color = "red") +
    ggplot2::geom_hline(yintercept = data_out$l, color = "red") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
  k
}

#' Interactive residuals plot in the field SpATS
#'
#' @param data_out res_spats data.frame
#'
#' @return ggplot
#' @noRd
#'
#' @examples
#' # in progress
plot_res_map <- function(data_out) {
  k <- dplyr::filter(data_out, !is.na(Classify)) %>%
    ggplot2::ggplot(ggplot2::aes(x = col, y = row, color = Classify)) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c("grey80", "red"))
  k
}

#' Interactive residuals vs fitted values SpATS
#'
#' @param data_out res_spats data.frame
#'
#' @return ggplot
#' @noRd
#'
#' @examples
#' # in progress
plot_res_fitted <- function(data_out) {
  k <- dplyr::filter(data_out, !is.na(Classify)) %>%
    ggplot2::ggplot(ggplot2::aes(x = fit, y = Residuals, color = Classify)) +
    ggplot2::geom_point(size = 2, alpha = 0.5, na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c("grey80", "red")) +
    ggplot2::xlab("Fitted Values") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
  k
}

#' @noRd
res_hist <- function(data_out) {
  hi <- graphics::hist(data_out[, "Residuals"], plot = FALSE)
  br <- hi$breaks
  p <- ggplot(data_out, aes(x = Residuals)) +
    geom_histogram(aes(y = ..density..), alpha = 0.8, breaks = c(br), na.rm = TRUE) +
    theme_bw() +
    geom_density(alpha = 0.5, na.rm = TRUE) +
    geom_vline(xintercept = c(data_out$u, data_out$l), linetype = 2, color = "red")
  p
}

#' @noRd
res_compare <- function(Model, variable, factor) {
  data <- Model$data
  data$Residuals <- stats::residuals(Model)
  data <- type.convert(data)
  if (factor) {
    data[, variable] <- as.factor(data[, variable])
    p <- ggplot(data, aes_string(x = variable, y = "Residuals", fill = variable)) +
      geom_boxplot(na.rm = TRUE) +
      theme_bw()
  } else {
    data[, variable] <- as.numeric(data[, variable])
    p <- ggplot(data, aes_string(x = variable, y = "Residuals")) +
      geom_point(size = 2, alpha = 0.5, color = "grey80", na.rm = TRUE) +
      theme_bw()
  }
  p
}

#' @noRd
check_gen_SpATS <- function(gen, data, check_gen = c("ci", "st", "wa")) {
  data <- as.data.frame(data)
  indx <- sum(check_gen %in% data[, gen]) >= 1
  if (indx) {
    genotypes_id <- as.character(data[, gen])
    data[, gen] <- as.factor(ifelse(genotypes_id %in% check_gen, NA, genotypes_id))
    data$checks <- as.factor(ifelse(genotypes_id %in% check_gen, genotypes_id, "_NoCheck"))
  } else {
    message("No checks in this trial")
  }
  return(data)
}

# MSA
#' @noRd
VarG_SpATS <- function(model) {
  gen <- model$model$geno$genotype
  gen_ran <- model$model$geno$as.random
  if (gen_ran) {
    vargen <- round(model$var.comp[gen], 3)
    names(vargen) <- "Var_Gen"
    return(vargen)
  } else {
    CV <- NA
    return(CV)
  }
}

#' @noRd
VarE_msa <- function(model) {
  v <- round(model$psi[1], 3)
  names(v) <- "Var_Res"
  return(v)
}

#' @noRd
weight_SpATS <- function(model) {
  rand <- model$model$geno$as.random
  if (rand) {
    return()
  }

  C_inv <- as.matrix(rbind(
    cbind(model$vcov$C11_inv, model$vcov$C12_inv), # Combine components into one matrix C
    cbind(model$vcov$C21_inv, model$vcov$C22_inv)
  ))
  gen_mat <- colnames(model$vcov$C11_inv)

  genotype <- model$model$geno$genotype
  dt <- SpATS::predict.SpATS(model, which = genotype, predFixed = "marginal") %>%
    droplevels() %>%
    dplyr::mutate_if(is.numeric, round, 3)
  gen_lvls <- as.factor(unique(as.character(dt[, genotype])))

  intc <- intersect(gen_mat, gen_lvls)
  diff <- setdiff(gen_lvls, gen_mat)

  vcov <- C_inv[c("Intercept", intc), c("Intercept", intc)]
  colnames(vcov)[1] <- rownames(vcov)[1] <- diff
  diag_vcov <- diag(vcov)

  L <- diag(ncol(vcov))
  dimnames(L) <- list(colnames(vcov), rownames(vcov))
  L[, 1] <- 1
  Se2 <- diag(L %*% vcov %*% t(L))

  data_weights <- data.frame(
    gen = names(diag_vcov),
    vcov = diag_vcov,
    inv_vcov = 1 / diag_vcov,
    weights = 1 / Se2
  )
  data_weights <- merge(
    x = dt,
    y = data_weights,
    by.x = genotype,
    by.y = "gen",
    sort = FALSE
  )
  data_weights <- data_weights[, c(genotype, "predicted.values", "standard.errors", "weights")] # "vcov","inv_vcov",

  return(list(
    vcov = vcov,
    diag = diag_vcov,
    diag_inv = 1 / diag_vcov,
    se2 = Se2,
    se = sqrt(Se2),
    data_weights = data_weights
  ))
}

#' @noRd
lrt_SpATS <- function(Model_nested, Model_full) {
  lo.lik1 <- Model_nested / -2
  lo.lik2 <- Model_full / -2

  d <- 2 * (lo.lik2 - lo.lik1)

  p.value1 <- round(1 - stats::pchisq(d, 1), 3)

  siglevel <- 0
  if (abs(p.value1) < 0.05) {
    siglevel <- "*"
  } else {
    siglevel <- "Not signif"
  }
  if (abs(p.value1) < 0.01) {
    siglevel <- "**"
  }
  if (abs(p.value1) < 0.001) {
    siglevel <- "***"
  }

  names(p.value1) <- "p-value"
  cat("\n========================================================")
  cat("\n", "Likelihood Ratio Test")
  cat("\n", "p-value =", p.value1, siglevel)
  cat("\n", "Sig.level: 0'***' 0.001 '**' 0.01 '*' 0.05 'Not signif' 1\n")
  cat("========================================================")
}

var_comp_SpATS <- function(object, which = "variances") {
  which <- match.arg(which)
  var.comp <- object$var.comp
  psi <- object$psi[1]
  nterms <- length(var.comp)
  model <- names(var.comp)
  col.names <- c("Variance", "SD", "log10(lambda)")
  row.names <- c(model, NA, "Residual")
  vc <- matrix(ncol = 3, nrow = nterms + 2, dimnames = list(
    row.names,
    col.names
  ))
  vc[, 1] <- c(sprintf("%.3e", var.comp), NA, sprintf(
    "%.3e",
    psi
  ))
  vc[, 2] <- c(sprintf("%.3e", sqrt(var.comp)), NA, sprintf(
    "%.3e",
    sqrt(psi)
  ))
  vc[, 3] <- c(sprintf("%.5f", log10(psi / var.comp)), NA, NA)
  eff.dim <- object$eff.dim
  dim <- object$dim
  dim.nom <- object$dim.nom
  tot_ed <- sum(eff.dim, na.rm = TRUE)
  tot_dim <- sum(dim, na.rm = TRUE)
  dim.new <- dim[match(names(eff.dim), names(dim))]
  dim.nom <- dim.nom[match(names(eff.dim), names(dim.nom))]
  type <- rep(NA, length = length(dim.new))
  type[(attr(dim, "random") & !attr(dim, "spatial"))[match(
    names(eff.dim),
    names(dim)
  )]] <- "R"
  type[(!attr(dim, "random") & !attr(dim, "spatial"))[match(
    names(eff.dim),
    names(dim)
  )]] <- "F"
  type[is.na(type)] <- "S"
  eff.dim.new <- eff.dim
  smooth.comp <- attr(object$terms$spatial, "term")
  if (paste(smooth.comp, "Global") %in% names(dim)) {
    dim.new <- c(dim.new, dim[paste(smooth.comp, "Global")])
    dim.nom <- c(dim.nom, dim[paste(smooth.comp, "Global")])
    eff.dim.new <- c(eff.dim.new, sum(eff.dim.new[grep(smooth.comp,
      names(eff.dim.new),
      fixed = TRUE
    )]))
    names(eff.dim.new)[length(eff.dim.new)] <- paste(
      smooth.comp,
      "Global"
    )
    type <- c(type, "S")
  }
  ord <- c(which(type == "F"), which(type == "R"), which(type ==
    "S"))
  eff.dim.new <- eff.dim.new[ord]
  dim.new <- dim.new[ord]
  dim.nom <- dim.nom[ord]
  model <- model[ord]
  type <- type[ord]
  nterms <- length(eff.dim.new)
  model <- names(eff.dim.new)
  Nobs <- object$nobs
  col.names <- c(
    "Effective", "Model", "Nominal", "Ratio",
    "Type"
  )
  row.names <- c(model, NA, "Total", "Residual", "Nobs")
  m <- matrix(ncol = 5, nrow = nterms + 4, dimnames = list(
    row.names,
    col.names
  ))
  m[, 1] <- c(sprintf("%.1f", eff.dim.new), NA, sprintf(
    "%.1f",
    tot_ed
  ), sprintf("%.1f", Nobs - tot_ed), sprintf(
    "%.0f",
    Nobs
  ))
  m[, 2] <- c(sprintf("%.0f", dim.new), NA, sprintf(
    "%.0f",
    tot_dim
  ), NA, NA)
  m[, 3] <- c(sprintf("%.0f", dim.nom), NA, sprintf(
    "%.0f",
    sum(dim.nom, na.rm = TRUE)
  ), NA, NA)
  m[, 4] <- c(sprintf("%.2f", eff.dim.new / dim.nom), NA, sprintf(
    "%.2f",
    tot_ed / sum(dim.nom, na.rm = TRUE)
  ), NA, NA)
  m[, 5] <- c(type, NA, NA, NA, NA)
  object$p.table.vc <- vc
  object$p.table.dim <- m
  class(object) <- "summary.SpATS"

  v_name <- as.character(na.omit(row.names(vc)))
  vc <- vc %>%
    as.data.frame() %>%
    type.convert() %>%
    dplyr::mutate_if(is.numeric, round, 3)
  vc <- vc[-c(nrow(vc) - 1), ]
  vc$Component <- v_name
  vc <- vc[, c(4, 1:3)]
  return(vc)
}

# SpATS coefficients
#' @noRd
coef_SpATS <- function(model) {
  # coefficients
  coef_spats <- model$coeff
  coef_random <- attr(coef_spats, "random")
  coef_spats <- data.frame(level = names(coef_spats), solution = coef_spats, coef_random, row.names = NULL)

  # diagonal c_Inverse
  C_inv <- as.matrix(rbind(cbind(model$vcov$C11_inv, model$vcov$C12_inv), cbind(model$vcov$C21_inv, model$vcov$C22_inv)))
  se <- sqrt(diag(C_inv))
  C_inv <- data.frame(level = names(se), std.error = se, row.names = NULL)
  coef_spats <- merge(coef_spats, C_inv, by = "level", all = TRUE)
  coef_spats <- coef_spats %>%
    dplyr::mutate(z.ratio = solution / std.error) %>%
    dplyr::mutate_if(is.numeric, round, 4) %>%
    dplyr::arrange(coef_random)
  coef_spats <- coef_spats[, c("level", "solution", "std.error", "z.ratio", "coef_random")]
  return(coef_spats)
}
