#' Print an object of class \code{checkAgri}
#'
#' @description Prints information about \code{check_design_met()} function.
#'
#' @aliases print.checkAgri
#' @usage \method{print}{checkAgri}(x, ...)
#' @param x An object fitted with the function \code{check_design_met()}.
#' @param ... Options used by the tibble package to format the output. See
#' `tibble::print()` for more details.
#' @author Johan Aparicio [aut]
#' @method print checkAgri
#' @return an object inheriting from class \code{checkAgri}.
#' @importFrom utils str
#' @export
#' @examples
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_met(
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
print.checkAgri <- function(x, ...) {
  cat("---------------------------------------------------------------------\n")
  cat("Summary Traits by Trial:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$summ_traits, ...)
  cat(
    "\n---------------------------------------------------------------------\n"
  )
  cat("Experimental Design Detected:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$exp_design_list)
  cat(
    "\n---------------------------------------------------------------------\n"
  )
  cat("Summary Experimental Design:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$exp_design_resum[, 1:9], ...)
  cat(
    "\n---------------------------------------------------------------------\n"
  )
  cat("Connectivity Matrix:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$connectivity_matrix)
  cat(
    "\n---------------------------------------------------------------------\n"
  )
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
#' @return an object inheriting from class \code{smaAgri}.
#' @importFrom utils head
#' @export
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_met(
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
  cat(
    "\n---------------------------------------------------------------------\n"
  )
  cat("Outliers Removed:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$outliers, ...)
  cat(
    "\n---------------------------------------------------------------------\n"
  )
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
#' @return an object inheriting from class \code{metAgri}.
#' @importFrom utils head
#' @export
#' @examples
#' \dontrun{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_met(
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
#' met_results <- met_analysis(out, progress = FALSE)
#' print(met_results)
#' }
print.metAgri <- function(x, ...) {
  cat("---------------------------------------------------------------------\n")
  cat("Trial Effects (BLUEs):\n")
  cat("---------------------------------------------------------------------\n")
  print(x$trial_effects, ...)
  cat(
    "\n---------------------------------------------------------------------\n"
  )
  cat("Heritability:\n")
  cat("---------------------------------------------------------------------\n")
  print(x$heritability, ...)
  cat(
    "\n---------------------------------------------------------------------\n"
  )
  cat("First Overall Predicted Values and Standard Errors (BLUPs):\n")
  cat("---------------------------------------------------------------------\n")
  print(head(x$overall_BLUPs))
  cat(
    "\n---------------------------------------------------------------------\n"
  )
  cat("Variance-Covariance Matrix:\n")
  cat("---------------------------------------------------------------------\n")
  lt <- names(x$VCOV)
  for (i in lt) {
    vc_str <- x$VCOV[[i]]$vc_model
    cat("\nCorrelation Matrix ", "('", vc_str, "'): ", i, "\n", sep = "")
    print(round(x$VCOV[[i]]$CORR, 2))
    cat("\nCovariance Matrix ", "('", vc_str, "'): ", i, "\n", sep = "")
    print(round(x$VCOV[[i]]$VCOV, 2))
    cat(
      "\n---------------------------------------------------------------------"
    )
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
#' executing the function \code{check_design_met()}
#' @param type A character string specifiying the type of plot. "connectivity",
#' "missing" or "boxplot".
#' @param ... Further graphical parameters. For future improvements.
#' @param axis_size Numeric input to define the axis size.
#' @param text_size Numeric input to define the text size.
#' @author Johan Aparicio [aut]
#' @method plot checkAgri
#' @return A ggplot object.
#' @export
#' @examples
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_met(
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
#' plot(results, type = "boxplot")
#'
plot.checkAgri <- function(x,
                           type = c("connectivity", "missing", "boxplot"),
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

  if (type == "boxplot") {
    g_plot <- x$data_design %>%
      dplyr::select(.data[[x$inputs$trial]], all_of(x$inputs$traits)) %>%
      tidyr::gather(
        data = .,
        key = "traits",
        value = "value",
        -.data[[x$inputs$trial]]
      ) %>%
      ggplot(
        aes(
          x = .data[[x$inputs$trial]],
          y = value,
          fill = .data[[x$inputs$trial]]
        )
      ) +
      geom_boxplot() +
      labs(title = "", y = "Boxplot", x = "Trials") +
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
      facet_wrap(~traits, scales = "free_y") +
      theme(
        plot.title = element_text(
          color = "black",
          size = axis_size,
          hjust = 0.5
        ),
        axis.title.x = element_text(color = "black", size = axis_size),
        axis.title.y = element_text(color = "black", size = axis_size),
        legend.position = "bottom"
      )
    g_plot
  }
  return(g_plot)
}


#' Plot an object of class \code{smaAgri}
#'
#' @description Create several plots for an object of class \code{smaAgri}
#' @aliases plot.smaAgri
#' @param x An object inheriting from class \code{smaAgri} resulting of
#' executing the function \code{single_trial_analysis()}
#' @param type A character string specifiying the type of plot. "summary" or
#' "correlation".
#' @param ... Further graphical parameters. For future improvements.
#' @param filter_traits An optional character vector to filter traits.
#' @param nudge_y_cv Vertical adjustment to nudge labels by when plotting CV
#' bars. Only works if the argument type is "summary". 3 by default.
#' @param nudge_y_h2 Vertical adjustment to nudge labels by when plotting h2
#' bars. Only works if the argument type is "summary". 0.07 by default.
#' @param horizontal If \code{FALSE}, the default, the labels are plotted
#' vertically. If \code{TRUE}, the labels are plotted horizontally.
#' @param theme_size Base font size, given in pts. 15 by default.
#' @param axis_size Numeric input to define the axis size.
#' @param text_size Numeric input to define the text size.
#' @author Johan Aparicio [aut]
#' @method plot smaAgri
#' @return A ggplot object.
#' @importFrom ggpubr ggarrange
#' @export
#' @examples
#' \donttest{
#' library(agridat)
#' library(agriutilities)
#' data(besag.met)
#' dat <- besag.met
#' results <- check_design_met(
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
#' plot(out, type = "summary", horizontal = TRUE)
#' plot(out, type = "correlation")
#' }
plot.smaAgri <- function(x,
                         type = c("summary", "correlation"),
                         filter_traits = NULL,
                         nudge_y_cv = 3,
                         nudge_y_h2 = 0.07,
                         horizontal = FALSE,
                         theme_size = 15,
                         axis_size = 8,
                         text_size = 4, ...) {
  type <- match.arg(type)
  if (type == "summary") {
    horizontal <- ifelse(test = isTRUE(horizontal), yes = 0, no = 90)

    tmp_out <- x$resum_fitted_model %>%
      {
        if (!is.null(filter_traits)) {
          filter(.data = ., trait %in% filter_traits)
        } else {
          .
        }
      } %>%
      select(-design) %>%
      tidyr::gather(key = "component", value = "value", -trait, -trial) %>%
      group_by(trait, component) %>%
      mutate(
        max = ifelse(
          test = value == max(value, na.rm = TRUE),
          yes = "max",
          no = ifelse(
            test = value == min(value, na.rm = TRUE),
            yes = "min",
            no = "med"
          )
        )
      ) %>%
      arrange(desc(value))

    A <- tmp_out %>%
      filter(component %in% c("CV")) %>%
      group_by(trial, trait, component) %>%
      ggplot(
        aes(
          x = trial,
          y = value,
          fill = max,
          label = paste0(round(value, 1), "%")
        )
      ) +
      geom_bar(stat = "identity", color = "black", alpha = 0.5) +
      facet_grid(
        facets = component ~ trait,
        scales = "free_y",
        switch = c("y")
      ) +
      theme_minimal(base_size = theme_size) +
      theme(
        axis.text.x = element_text(hjust = 1, angle = 90, size = axis_size),
        axis.text.y = element_text(size = axis_size),
        legend.position = "none"
      ) +
      scale_fill_manual(values = c("red", "grey", "steelblue")) +
      geom_text(
        aes(color = max),
        nudge_y = nudge_y_cv,
        angle = horizontal,
        size = text_size
      ) +
      scale_color_manual(values = c("darkred", "black", "blue")) +
      labs(x = "", y = "")

    B <- tmp_out %>%
      filter(component %in% c("heritability")) %>%
      group_by(trial, trait, component) %>%
      mutate(nudge = value + value * 0.15) %>%
      ggplot(
        aes(x = trial, y = value, fill = max, label = round(value, 2))
      ) +
      geom_bar(stat = "identity", color = "black", alpha = 0.5) +
      facet_grid(
        facets = component ~ trait,
        scales = "free_y",
        switch = c("y")
      ) +
      theme_minimal(base_size = theme_size) +
      theme(
        axis.text.x = element_text(hjust = 1, angle = 90, size = axis_size),
        axis.text.y = element_text(size = axis_size),
        legend.position = "none"
      ) +
      scale_fill_manual(values = c("red", "grey", "steelblue")) +
      geom_text(
        aes(color = max),
        nudge_y = nudge_y_h2,
        angle = horizontal,
        size = text_size
      ) +
      scale_color_manual(values = c("darkred", "black", "blue")) +
      labs(x = "", y = "") +
      ylim(c(NA, 1.2))

    C <- ggarrange(A, B, ncol = 1)
  }

  if (type == "correlation") {
    traits <- x$blues_blups %>%
      group_by(trait) %>%
      summarise(n = n_distinct(trial)) %>%
      {
        if (!is.null(filter_traits)) {
          filter(.data = ., trait %in% filter_traits)
        } else {
          .
        }
      } %>%
      filter(n > 1) %>%
      pull(trait)

    s <- list()
    for (i in traits) {
      tmp <- x$blues_blups %>%
        filter(trait %in% i) %>%
        droplevels() %>%
        select(genotype, trial, BLUPs) %>%
        spread(trial, BLUPs) %>%
        select(-genotype)
      h2 <- x$resum_fitted_model %>%
        filter(trait %in% i) %>%
        droplevels() %>%
        pull(heritability, name = trial)

      s[[i]] <- gg_cor(
        data = tmp,
        colours = c("#db4437", "white", "#4285f4"),
        Diag = h2,
        label_size = text_size
      ) +
        ggtitle(i)
    }
    C <- ggarrange(plotlist = s)
  }

  return(C)
}
