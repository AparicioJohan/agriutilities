#' @noRd
#' @keywords internal
"circleFun" <- function(center = c(0, 0), diameter = 1, npoints = 100) {
  r <- diameter / 2
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


#' Correlation Covariance Heatmap
#'
#' @param matrix A numeric matrix.
#' @param corr A logical value indicating if the matrix is in a scaled form
#'  (\code{TRUE} by default, correlation matrix)
#' @param size A numeric value to define the letter size.
#' @param digits A numeric integer to define the number of digits to plot.
#'
#' @return A ggplot object showing the upper triangular elements of the matrix.
#' @export
#'
#' @examples
#' \donttest{
#' library(agriutilities)
#' data(iris)
#' M <- cor(iris[, -5])
#' covcor_heat(matrix = M, corr = TRUE)
#' }
covcor_heat <- function(matrix, corr = TRUE, size = 4, digits = 3) {
  matrix <- round(x = matrix, digits = 3)

  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat) {
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
  }

  reorder_cormat <- function(cormat) {
    # Use correlation between variables as distance
    dd <- stats::as.dist((1 - cormat) / 2)
    hc <- stats::hclust(dd)
    cormat <- cormat[hc$order, hc$order]
  }
  if (corr) {
    matrix <- reorder_cormat(matrix)
  }
  upper_tri <- get_upper_tri(matrix)
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

  u <- -1
  m <- 0
  l <- 1
  main <- "Correlation"
  col_pallete <- c("#db4437", "white", "#4285f4")
  col_letter <- "black"

  if (isFALSE(corr)) {
    u <- min(matrix, na.rm = TRUE)
    l <- max(matrix, na.rm = TRUE)
    m <- u + (l - u) / 2
    main <- "Covariance"
    col_pallete <- c("#440154", "#21908C", "#FDE725")
    col_letter <- "white"
  }

  melted_cormat$Var1 <- as.factor(melted_cormat$Var1)
  melted_cormat$Var2 <- as.factor(melted_cormat$Var2)

  ggheatmap <-
    ggplot2::ggplot(
      data = melted_cormat,
      mapping = ggplot2::aes(Var2, Var1, fill = value)
    ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = col_pallete[1],
      high = col_pallete[3],
      mid = col_pallete[2],
      midpoint = m,
      limit = c(u, l),
      space = "Lab",
      name = main
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        vjust = 1,
        size = 12,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = 12)
    )

  plot <- ggheatmap +
    ggplot2::geom_text(
      mapping = ggplot2::aes(Var2, Var1, label = value),
      color = col_letter,
      size = size
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        barwidth = 7,
        barheight = 1,
        title.position = "top",
        title.hjust = 0.5
      )
    )

  return(plot)
}


#' Factor Analytic Summary
#'
#' @param model Factor Analytic Model (ASReml object)
#' @param trial A character string indicating the column in data that contains
#' trials.
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param BLUEs_trial A data.frame containing BLUEs for each trial.
#' @param mult_fa1 A constant to multiply the first loading. Must be 1 or -1.
#'  (-1 by default)
#' @param mult_fa2 A constant to multiply the second loading. Must be 1 or -1.
#' (1 by default)
#' @param filter_score A numeric value to filter genotypes by the distance from
#' the origin.
#' @param k_biplot A numeric value to multiply the scores in the biplot.
#' @param size_label_var A numeric value to define the label size for the
#' variables.
#' @param alpha_label_var A numeric value between (0,1) to define the label for
#' the variables.
#' @param size_label_ind A numeric value to define the label size for the
#' individuals.
#' @param alpha_label_ind A numeric value between (0,1) to define the label for
#'  the individuals.
#' @param size_arrow A numeric value to define the arrow size.
#' @param alpha_arrow A numeric value between (0,1) to define the arrow.
#' @param base_size A numeric value to define the base size.
#'
#' @return list with  loadings = L, loading_star,
#' Gvar, Cmat, summary_loadings, paf_site, var_tot, scores,
#' plots = list(loadings, biplot,  biplot_scaled,  loadings_c)
#' @export
#'
#' @examples
#' \donttest{
#' library(dplyr)
#' library(asreml)
#' library(agridat)
#' data(besag.met)
#' dat <- besag.met
#'
#' dat <- dat %>% arrange(county)
#' model <- asreml(
#'   fixed = yield ~ 1 + county,
#'   random = ~ fa(county, 2):gen + county:rep + diag(county):rep:block,
#'   residual = ~ dsum(~ units | county),
#'   data = dat,
#'   na.action = list(x = "include", y = "include"),
#'   trace = 0
#' )
#'
#' pp <- predict(model, classify = "county")$pvals
#' fa_summary(
#'   model = model,
#'   trial = "county",
#'   genotype = "gen",
#'   BLUEs_trial = pp,
#'   mult_fa1 = -1,
#'   mult_fa2 = -1,
#'   filter_score = 1,
#'   k_biplot = 10,
#'   size_label_var = 3,
#'   alpha_label_var = 0.5,
#'   size_label_ind = 3,
#'   alpha_label_ind = 0.8,
#'   size_arrow = 0.2,
#'   alpha_arrow = 0.1
#' )
#' }
#' @import ggplot2 ggrepel
fa_summary <- function(model = NULL,
                       trial = "trial",
                       genotype = "genotype",
                       BLUEs_trial = NULL,
                       mult_fa1 = -1,
                       mult_fa2 = 1,
                       filter_score = 1.5,
                       k_biplot = 1,
                       size_label_var = 2,
                       alpha_label_var = 0.2,
                       size_label_ind = 2,
                       alpha_label_ind = 0.8,
                       size_arrow = 0.2,
                       alpha_arrow = 0.2,
                       base_size = 12) {
  vars <- summary(model)$varcomp
  vars <- data.frame(effect = rownames(vars), vars, check.names = FALSE)

  # Loading by Trial
  fa1_loadings <- vars[grepl("!fa1$", vars$effect), "component"]
  fa2_loadings <- vars[grepl("!fa2$", vars$effect), "component"]
  L <- as.matrix(cbind(fa1_loadings, fa2_loadings))
  svd_L <- svd(L)
  L_star <- L %*% svd_L$v
  psi <- vars[grepl("!var$", vars$effect), "component"]
  Gvar <- L_star %*% t(L_star) + diag(psi)
  Cmat <- cov2cor(Gvar)
  VarTot <- sum(diag(L_star %*% t(L_star))) / sum(diag(Gvar))
  ns <- nlevels(model$mf[, trial])
  k <- 2
  snam <- levels(model$mf[, trial])
  paf_site <- matrix(0, nrow = ns, ncol = k)
  dimnames(paf_site) <- list(snam, paste("fac", 1:k, sep = "_"))
  for (i in 1:k) {
    paf_site[, i] <- 100 * diag(L_star[, i] %*% t(L_star[, i])) /
      diag(L_star %*% t(L_star) + diag(psi))
  }
  if (k > 1) {
    all <- 100 * diag(L_star %*% t(L_star)) /
      diag(L_star %*% t(L_star) + diag(psi))
    paf_site <- cbind(paf_site, all)
  }
  VarTot <- VarTot
  VarGenEnv <- diag(L_star %*% t(L_star) + diag(psi))
  TotGenVar <- sum(VarGenEnv)
  VarFA1 <- sum(VarGenEnv * paf_site[, 1]) / 100
  VarFA2 <- sum(VarGenEnv * paf_site[, 2]) / 100
  PerVarFA1 <- VarFA1 / TotGenVar
  PerVarFA2 <- VarFA2 / TotGenVar
  percentg <- round(c(PerVarFA1, PerVarFA2) * 100, 2)
  L_star[, 1] <- L_star[, 1] * mult_fa1
  L_star[, 2] <- L_star[, 2] * mult_fa2
  faComp <- data.frame(
    site = snam,
    fa1 = L_star[, 1],
    fa2 = L_star[, 2],
    psi = psi,
    Vg = diag(Gvar),
    BLUE = BLUEs_trial$predicted.value
  )
  sqrt_diag <- sqrt(diag(Gvar))
  faComp[, c("fa1_scaled", "fa2_scaled")] <- diag(1 / sqrt_diag) %*% L_star
  row.names(L) <- row.names(L_star) <- snam
  row.names(Gvar) <- colnames(Gvar) <- row.names(Cmat) <- colnames(Cmat) <- snam
  # Scores by Genotype
  coef_fa <- coef(model)$random
  fa1_scores <- coef_fa[grep(paste0("Comp1:", genotype), rownames(coef_fa)), ]
  fa2_scores <- coef_fa[grep(paste0("Comp2:", genotype), rownames(coef_fa)), ]
  names(fa1_scores) <- sub(
    pattern = paste0("fa(", trial, ", 2)_Comp1:", genotype, "_"),
    replacement = "",
    x = names(fa1_scores),
    fixed = TRUE
  )
  names(fa2_scores) <- sub(
    pattern = paste0("fa(", trial, ", 2)_Comp2:", genotype, "_"),
    replacement = "",
    x = names(fa2_scores),
    fixed = TRUE
  )
  f_scores <- rbind(as.matrix(fa1_scores), as.matrix(fa2_scores))
  nGenotype <- nlevels(model$mf[, genotype])
  f_star <- kronecker(t(svd_L$v), diag(nGenotype)) %*% f_scores
  rownames(f_star) <- rownames(f_scores)
  fa12_scores <- merge(
    x = f_star[1:nGenotype, 1],
    y = f_star[(nGenotype + 1):(2 * nGenotype), 1],
    by = "row.names"
  )
  names(fa12_scores) <- c("Genotype", "fa1", "fa2")
  fa12_scores$fa1 <- fa12_scores$fa1 * mult_fa1
  fa12_scores$fa2 <- fa12_scores$fa2 * mult_fa2
  fa12_scores$distance_orig <- sqrt(fa12_scores$fa1^2 + fa12_scores$fa2^2)
  fa12_scores$Score <- ifelse(
    test = fa12_scores$distance_orig > filter_score,
    yes = 1,
    no = 0
  )
  # Plots
  d <- data.frame(
    x = rep(0, nrow(L_star)),
    y = rep(0, nrow(L_star)),
    vx = L_star[, 1],
    vy = L_star[, 2]
  )
  loadings <- faComp %>%
    ggplot(aes(x = fa1, y = fa2)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_text_repel(
      mapping = aes(label = site),
      nudge_y = 0.05,
      nudge_x = -0.03,
      force = 1,
      alpha = alpha_label_var,
      size = size_label_var
    ) +
    ggtitle(
      label = paste0(
        "Environment Factor Loadings ",
        "(", sum(percentg), "%)"
      ),
      subtitle = NULL
    ) +
    xlab(paste0("FA1 loading ", "(", percentg[1], "%)")) +
    ylab(paste0("FA2 loading ", "(", percentg[2], "%)")) +
    theme_bw(base_size = base_size) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(
      data = d,
      mapping = aes(
        x = x,
        y = y,
        xend = x + vx,
        yend = y + vy
      ),
      arrow = arrow(),
      color = "black",
      size = size_arrow,
      alpha = alpha_arrow
    )

  biplot <- faComp %>%
    ggplot(aes(x = fa1, y = fa2)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_text_repel(
      mapping = aes(label = site),
      nudge_y = 0.05,
      nudge_x = -0.03,
      force = 1,
      alpha = alpha_label_var,
      size = size_label_var
    ) +
    ggtitle("Environment Factor Loadings") +
    xlab("FA1 loading") +
    ylab("FA2 loading") +
    theme_bw(base_size = base_size) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(
      data = d,
      mapping = aes(x = x, y = y, xend = x + vx, yend = y + vy),
      arrow = arrow(),
      color = "black",
      size = size_arrow,
      alpha = alpha_arrow
    ) +
    geom_label_repel(
      data = subset(fa12_scores, Score == 1),
      mapping = aes(label = Genotype, x = k_biplot * fa1, y = k_biplot * fa2),
      colour = "red",
      segment.colour = "red",
      size = size_label_ind,
      alpha = alpha_label_ind
    )

  d_scaled <- data.frame(
    x = rep(0, nrow(faComp)),
    y = rep(0, nrow(faComp)),
    vx = faComp[, "fa1_scaled"],
    vy = faComp[, "fa2_scaled"]
  )
  biplot_scaled <- faComp %>%
    ggplot(aes(x = fa1_scaled, y = fa2_scaled)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_text_repel(
      mapping = aes(label = site),
      nudge_y = 0.05,
      nudge_x = -0.03,
      force = 1,
      alpha = alpha_label_var,
      size = size_label_var
    ) +
    ggtitle("Environment Factor Loadings") +
    xlab("FA1 loading") +
    ylab("FA2 loading") +
    theme_bw(base_size = base_size) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(
      data = d_scaled,
      mapping = aes(x = x, y = y, xend = x + vx, yend = y + vy),
      arrow = arrow(),
      color = "black",
      size = size_arrow,
      alpha = alpha_arrow
    ) +
    geom_label_repel(
      data = subset(fa12_scores, Score == 1),
      mapping = aes(label = Genotype, x = 1 * fa1, y = 1 * fa2),
      colour = "red",
      segment.colour = "red",
      size = size_label_ind,
      alpha = alpha_label_ind
    )

  circle <- circleFun(
    center = c(0, 0),
    diameter = 2,
    npoints = 100
  )
  var_centered <- faComp %>%
    ggplot(aes(x = fa1_scaled, y = fa2_scaled)) +
    geom_point(aes(colour = Vg, size = BLUE)) +
    scale_colour_gradient(low = "pink", high = "blue") +
    geom_text_repel(
      mapping = aes(label = site),
      nudge_y = 0.05,
      nudge_x = -0.03,
      force = 1,
      size = size_label_var,
      alpha = alpha_label_var
    ) +
    ggtitle(label = "Environment Factor Loadings") +
    xlab(label = "FA1 loading") +
    ylab(label = "FA2 loading") +
    theme_bw(base_size = base_size) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(
      data = d_scaled,
      mapping = aes(x = x, y = y, xend = x + vx, yend = y + vy),
      arrow = arrow(),
      color = "black",
      size = size_arrow,
      alpha = alpha_arrow
    ) +
    geom_path(data = circle, mapping = aes(x, y)) +
    coord_fixed()

  results <- list(
    loadings = L,
    loading_star = L_star,
    Gvar = Gvar,
    Cmat = Cmat,
    summary_loadings = faComp,
    paf_site = paf_site,
    var_tot = percentg,
    scores = fa12_scores,
    plots = list(
      loadings = loadings,
      biplot = biplot,
      biplot_scaled = biplot_scaled,
      loadings_c = var_centered
    )
  )
  return(results)
}
