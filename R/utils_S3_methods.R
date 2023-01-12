#' Print an object of class \code{checkAgri}
#'
#' @description Prints information about \code{check_design_MET()} function.
#'
#' @aliases print.checkAgri
#' @usage \method{print}{checkAgri}(x, ...)
#' @param x An object fitted with the function \code{check_design_MET()}.
#' @param ... Options used by the tibble package to format the output. See
#' `tibble::print()` for more details.
#' @author Johan Aparicio [aut]
#' @method print checkAgri
#' @importFrom utils str
#' @export
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_MET(
#'   data = dat,
#'   genotype = "gen",
#'   trial = "county",
#'   traits = c("yield"),
#'   rep = "rep",
#'   block = "block",
#'   col = "col",
#'   row = "row"
#' )
#' print(results)
#' }
print.checkAgri <- function(x, ...) {
  cat("---------------------------------------------------------------------\n")
  cat("Summary Traits by Trial:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$summ_traits, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("Experimental Design Detected:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$exp_design_list)
  cat("\n---------------------------------------------------------------------\n")
  cat("Summary Experimental Design:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$exp_design_resum[, 1:9], ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("Connectivity Matrix:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$connectivity_matrix)
  cat("\n---------------------------------------------------------------------\n")
  cat("Filters Applied:\n")
  cat("---------------------------------------------------------------------\n")
  str(x$filter)
  cat("\n")
}

#' Print an object of class \code{smaAgri}
#'
#' @description Prints information about \code{single_trial_analysis} function.
#'
#' @aliases print.smaAgri
#' @usage \method{print}{smaAgri}(x, ...)
#' @param x An object fitted with the function \code{single_trial_analysis()}.
#' @param ... Options used by the tibble package to format the output. See
#' `tibble::print()` for more details.
#' @author Johan Aparicio [aut]
#' @method print smaAgri
#' @importFrom utils head
#' @export
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_MET(
#'   data = dat,
#'   genotype = "gen",
#'   trial = "county",
#'   traits = c("yield"),
#'   rep = "rep",
#'   block = "block",
#'   col = "col",
#'   row = "row"
#' )
#' out <- single_trial_analysis(results, progress = FALSE)
#' print(out)
#' }
print.smaAgri <- function(x, ...) {
  cat("---------------------------------------------------------------------\n")
  cat("Summary Fitted Models:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$resum_fitted_model, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("Outliers Removed:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$outliers, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("First Predicted Values and Standard Errors (BLUEs/BLUPs):\n")
  cat("---------------------------------------------------------------------\n")
  print(head(x$blues_blups))
  cat("\n")
}

#' Print an object of class \code{metAgri}
#'
#' @description Prints information about \code{met_analysis()} function.
#'
#' @aliases print.metAgri
#' @usage \method{print}{metAgri}(x, ...)
#' @param x An object fitted with the function \code{met_analysis()}.
#' @param ... Options used by the tibble package to format the output. See
#' `tibble::print()` for more details.
#' @author Johan Aparicio [aut]
#' @method print metAgri
#' @importFrom utils head
#' @export
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_MET(
#'   data = dat,
#'   genotype = "gen",
#'   trial = "county",
#'   traits = c("yield"),
#'   rep = "rep",
#'   block = "block",
#'   col = "col",
#'   row = "row"
#' )
#' out <- single_trial_analysis(results, progress = FALSE)
#' met_results <- met_analysis(out)
#' print(met_results)
#' }
print.metAgri <- function(x, ...) {
  cat("---------------------------------------------------------------------\n")
  cat("Trial Effects (BLUEs):\n")
  cat("---------------------------------------------------------------------\n")
  print(x$trial_effects, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("Heritability:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$heritability, ...)
  cat("\n---------------------------------------------------------------------\n")
  cat("First Overall Predicted Values and Standard Errors (BLUPs):\n")
  cat("---------------------------------------------------------------------\n")
  print(head(x$overall_BLUPs))
  cat("\n---------------------------------------------------------------------\n")
  cat("Variance-Covariance Matrix:\n")
  cat("---------------------------------------------------------------------\n")
  lt <- names(x$VCOV)
  for (i in lt) {
    vc_str <- x$VCOV[[i]]$vc.model
    cat("\nCorrelation Matrix ", "('", vc_str, "'): ", i, "\n", sep = "")
    print(round(x$VCOV[[i]]$CORR, 2))
    cat("\nCovariance Matrix ", "('", vc_str, "'): ", i, "\n", sep = "")
    print(round(x$VCOV[[i]]$VCOV, 2))
    cat("\n---------------------------------------------------------------------")
  }
  cat("\nFirst Stability Coefficients:\n")
  cat("---------------------------------------------------------------------\n")
  print(head(x$stability))
  cat("\n")
}


#' Plot an object of class \code{checkAgri}
#'
#' @description Create several plots for an object of class \code{checkAgri}
#' @aliases plot.checkAgri
#' @param x An object inheriting from class \code{checkAgri} resulting of
#' executing the function \code{check_design_MET()}
#' @param type A character string specifiying the type of plot. "connectivity"
#' or "missing".
#' @param ... Further graphical parameters. For future improvements.
#' @param axis_size Numeric input to define the axis size.
#' @param text_size Numeric input to define the text size.
#' @author Johan Aparicio [aut]
#' @method plot checkAgri
#' @export
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_MET(
#'   data = dat,
#'   genotype = "gen",
#'   trial = "county",
#'   traits = c("yield"),
#'   rep = "rep",
#'   block = "block",
#'   col = "col",
#'   row = "row"
#' )
#' plot(results, type = "missing")
#' }
#'
plot.checkAgri <- function(x,
                           type = c("connectivity", "missing"),
                           axis_size = 15,
                           text_size = 5, ...) {
  type <- match.arg(type)
  if (type == "connectivity") {
    reorder_cormat <- function(cormat) {
      dd <- stats::dist(cormat)
      hc <- stats::hclust(dd)
      cormat <- cormat[hc$order, hc$order]
    }
    conn <- reorder_cormat(x$connectivity_matrix)
    order_trials_x <- rownames(conn)
    order_trials_y <- colnames(conn)
    MM <- conn %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "trial_x") %>%
      tidyr::gather(data = ., key = "trial_y", value = "n", -trial_x) %>%
      dplyr::mutate(n = ifelse(n == 0, NA, n)) %>%
      dplyr::filter(!is.na(n)) %>%
      dplyr::mutate(
        trial_x = factor(trial_x, levels = order_trials_x),
        trial_y = factor(trial_y, levels = order_trials_y)
      )
    colours <- c("#db4437", "white", "#4285f4")
    g_plot <- MM %>%
      ggplot(
        aes(
          x = trial_x,
          y = trial_y,
          fill = n
        )
      ) +
      geom_tile(color = "gray") +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = axis_size) +
      geom_text(
        aes(
          x = trial_x,
          y = trial_y,
          label = n
        ),
        color = "black",
        size = text_size
      ) +
      scale_fill_gradient2(
        low = colours[1],
        mid = colours[2],
        high = colours[3],
        midpoint = median(MM$n),
      ) +
      theme(
        axis.text.x = element_text(angle = 40, hjust = 1, size = axis_size),
        axis.text.y = element_text(size = axis_size),
        legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.major = element_blank()
      ) +
      labs(title = "Connectivity Matrix")
    g_plot

  }

  if (type == "missing") {
    g_plot <- x$summ_traits %>%
      ggplot(
        aes(x = .data[[x$inputs$trial]], y = miss_perc)
      ) +
      geom_point(size = 2) +
      geom_segment(
        aes(
          x = .data[[x$inputs$trial]],
          xend = .data[[x$inputs$trial]],
          y = 0,
          yend = miss_perc
        )
      ) +
      labs(title = "", y = "Proportion of Missing Data", x = "Trials") +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.6,
          color = "black",
          size = axis_size
        ),
        axis.text.y = element_text(size = axis_size)
      ) +
      facet_wrap(~traits, scales = "free_x") +
      theme(
        plot.title = element_text(
          color = "black",
          size = axis_size,
          hjust = 0.5
        ),
        axis.title.x = element_text(color = "black", size = axis_size),
        axis.title.y = element_text(color = "black", size = axis_size),
        legend.position = "bottom"
      ) +
      geom_hline(yintercept = 0.5, color = "black", linetype = 2) +
      ylim(c(0, 1))
    g_plot
  }
  return(g_plot)
}
