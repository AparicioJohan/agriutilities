#' Triangular Correlation Plot
#'
#' @description Return a ggplot object to plot a triangular correlation figure
#' between 2 or more variables.
#'
#' @param data A data.frame with numerical columns for each variable to be
#' compared.
#' @param colours A vector of size three with the colors to be used for values
#' -1, 0 and 1.
#' @param blackLabs A numeric vector of size two, with min and max correlation
#' coefficient.
#' @param show_signif Logical scalar. Display significance values ?
#' @param p_breaks Passed to function 'cut'. Either a numeric vector of two or
#' more unique cut points or a single number (greater than or equal to 2) giving
#' the number of intervals into which x is to be cut.
#' @param p_labels Passed to function 'cut'. labels for the levels of the
#' resulting category. By default, labels are constructed using "(a,b]" interval
#' notation. If \code{p_labels = FALSE}, simple integer codes are returned
#' instead of a factor.
#' @param show_diagonal Logical scalar. Display main diagonal values ?
#' @param diag A named vector of labels to display in the main diagonal. The
#' names are used to place each value in the corresponding coordinates of the
#' diagonal. Hence, these names must be the same as the colnames of data.
#' @param return_table Return the table to display instead of a ggplot object.
#' @param return_n Return plot with shared information.
#' @param adjusted Use the adjusted p values for multiple testing instead of
#' raw coeffs. \code{TRUE} by default.
#' @param label_size Numeric value indicating the label size. 3 by default.
#' @param method method="pearson" is the default value.
#'  The alternatives to be passed to cor are "spearman" and "kendall".
#'  These last two are much slower, particularly for big data sets.
#'
#' @return A ggplot object containing a triangular correlation figure with all
#' numeric variables in data. If return_table is \code{TRUE}, the table used to
#' produce the figure is returned instead.
#' @export
#'
#' @examples
#' library(agriutilities)
#' data(iris)
#' gg_cor(
#'   data = iris,
#'   colours = c("#db4437", "white", "#4285f4"),
#'   label_size = 6
#' )
#' @author Daniel Ariza, Johan Aparicio.
#' @importFrom stats na.omit
gg_cor <- function(data,
                   colours = c("#db4437", "white", "#4285f4"),
                   blackLabs = c(-0.7, 0.7),
                   show_signif = TRUE,
                   p_breaks = c(0, .001, .01, .05, Inf),
                   p_labels = c("***", "**", "*", "ns"),
                   show_diagonal = FALSE,
                   diag = NULL,
                   return_table = FALSE,
                   return_n = FALSE,
                   adjusted = TRUE,
                   label_size = 3,
                   method = "pearson") {
  # Drop non numeric columns in the dataset
  if (sum(!sapply(data, is.numeric))) {
    message(
      "Dropping non-numeric columns in the dataset:\n",
      paste(names(which(!sapply(data, is.numeric))), collapse = "\t")
    )
    data <- data[, sapply(data, is.numeric)]
  }
  # Calculate corr-coeffs and p values
  cors <- psych::corr.test(data, use = "pairwise.complete.obs", method = method)
  # Use the adjusted p values for multiple testing instead of raw coeffs
  if (adjusted) cors$p <- t(cors$p)
  # Keep only the matrices with correlation coefficients, p values and N shared
  # samples
  cors <- cors[c(1, 2, 4)]
  # Make sure you have a full matrix of N shared samples
  if (is.vector(cors$n)) {
    cors$n <- matrix(
      data = cors$n,
      ncol = ncol(cors$p),
      nrow = nrow(cors$p),
      dimnames = dimnames(cors$p)
    )
  }
  # For each matrix, do ...
  cors <- lapply(cors, function(x) {
    # Keep the upper triangle of the matrix
    x[upper.tri(x)] <- NA
    # Transpose the matrix to plot the lower triangle
    x <- as.data.frame(t(x))
    # Reshape the matrix to tidy format
    x[, "col"] <- colnames(x)
    x <- tidyr::gather(data = x, key = "row", value = "value", -col)
    colnames(x) <- c("col", "row", "value")
    # Round coefficients
    x$name <- round(x$value, 2)
    # Sort the x axis according to data column order
    x$col <- factor(x$col, levels = colnames(data))
    # Reverse the y axis for a triangle plot from top-left to bottom-right
    x$row <- factor(x$row, levels = rev(colnames(data)))
    # Remove NAs
    x <- na.omit(x)
  })

  # Combine both dataframes with p values and corr coefficients
  cors <- merge(
    x = merge(x = cors$r, y = cors$p, by = c("col", "row")),
    y = cors$n,
    by = c("col", "row")
  )
  # Keep x, y, p val and corr-coefficients columns
  cors <- cors[, c(1, 2, 4, 5, 7)]

  if (return_n) {
    if (return_table) {
      return(cors)
    }
    cors$cols <- scale(cors$value, center = TRUE, scale = TRUE)
    cors$cols <- ifelse(abs(cors$cols) < 2, "black", "white")
    p <- ggplot2::ggplot(
      data = cors,
      ggplot2::aes(x = col, y = row, fill = value)
    ) +
      ggplot2::geom_tile(color = "gray") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::geom_text(
        ggplot2::aes(x = col, y = row, label = value),
        color = cors$cols,
        size = label_size
      ) +
      ggplot2::scale_fill_gradient(
        low = colours[2],
        high = colours[3],
        limits = c(0, max(cors$value))
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
        legend.position = "none",
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank()
      )
    return(p)
  }

  if (show_signif) {
    # Create a categorical variable for p values as defined by p_breaks
    cors$signi <- cut(
      x = cors$value.y, right = FALSE,
      breaks = p_breaks, labels = p_labels
    )
    # Join corr-coeff and p-value to display it as a label for each tile
    cors$label <- paste(cors$name.x, cors$sign, sep = "\n")
  } else {
    # The label for each tile is the corr-coeff only
    cors$label <- cors$name.x
  }

  # If there are user-specified values to display in the diagonal
  if (!is.null(diag)) {
    # Check the names in diag are the same than colnames of data
    if (sum(!names(diag) %in% colnames(data))) {
      warning(
        "These elements in 'diag' do not correspond to column names in
        'data':\n",
        paste(names(diag)[!names(diag) %in% colnames(data)],
          collapse = "\t"
        )
      )
    }
    # The tiles of the diagonal are gray
    cors[cors$col == cors$row, "name.x"] <- NA
    # Get the name of x and y levels
    d <- as.character(cors[cors$col == cors$row, "row"])
    # Modify the elements of the diagonal and make sure they are displayed
    cors[cors$col == cors$row, "label"] <- diag[d]
    show_diagonal <- TRUE
  }

  # Remove the elements of the main diagonal if you don't want to display
  if (!show_diagonal) cors <- cors[cors$col != cors$row, ]

  # Show darker tiles with white labels for clarity
  cors$txtCol <- ifelse(
    test = cors$name.x > blackLabs[1] & cors$name.x < blackLabs[2],
    yes = "black",
    no = "white"
  )
  # Do not show tile labels for empty tiles.
  # Make tile labels of the diagonal white
  cors$txtCol[is.na(cors$txtCol)] <- "white"

  if (return_table) {
    return(cors)
  }

  p <- ggplot2::ggplot(
    data = cors,
    ggplot2::aes(x = col, y = row, fill = name.x)
  ) +
    ggplot2::geom_tile(color = "gray") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::geom_text(
      ggplot2::aes(x = col, y = row, label = label),
      color = cors$txtCol,
      size = label_size
    ) +
    ggplot2::scale_fill_gradient2(
      low = colours[1],
      mid = colours[2],
      high = colours[3]
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      legend.position = "none",
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank()
    )
  return(p)
}
