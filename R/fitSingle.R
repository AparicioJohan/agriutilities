#' @noRd
#' @keywords internal
fit_STA <- function(results, trait, design, remove_outliers, engine, progress) {
  design_td <- switch(design,
    "res_row_col" = "res.rowcol",
    "row_col" = "rowcol",
    "alpha_lattice" = "res.ibd",
    "rcbd" = "rcbd"
  )
  if (design %in% "alpha_lattice" || design %in% "rcbd") {
    results$inputs$row <- NULL
    results$inputs$col <- NULL
    spatial <- FALSE
  } else {
    spatial <- TRUE
  }
  m_models <- list()
  results <- results
  exp_to_remove <- results$filter[[trait]]$trials_to_remove

  trials <- results$data_design %>%
    filter(
      exp_design == design &
        !.data[[results$inputs$trial]] %in% exp_to_remove
    ) %>%
    pull(results$inputs$trial) %>%
    as.character() %>%
    unique()

  data <- results$data_design %>%
    filter(.data[[results$inputs$trial]] %in% trials) %>%
    droplevels()

  if (nrow(data) <= 0) {
    message("There is no data to fit any model for some experiments in:", trait)
    return()
  }
  td <- createTD(
    data = data,
    genotype = results$inputs$genotype,
    trial = results$inputs$trial,
    repId = results$inputs$rep,
    subBlock = results$inputs$block,
    rowCoord = results$inputs$row,
    colCoord = results$inputs$col,
    trDesign = design_td
  )
  m_models <- fitTD(
    TD = td,
    traits = trait,
    what = c("fixed", "random"),
    spatial = spatial,
    progress = progress,
    engine = engine
  )
  # Residuals
  outliers_td <- outlierSTA(
    STA = m_models,
    traits = trait,
    rLimit = 3,
    verbose = FALSE
  )
  if (is.null(outliers_td$outliers) || !remove_outliers) {
    outliers_td <- NULL
  }
  # Cleaning
  if (!is.null(outliers_td$outliers) && remove_outliers) {
    outliers_td <- outliers_td %>%
      .[["outliers"]] %>%
      dplyr::select(trial, genotype, id, outlier)
    data_clean <- data %>%
      merge(
        x = .,
        y = outliers_td,
        by.x = c(results$inputs$trial, results$inputs$genotype, "id"),
        by.y = c("trial", "genotype", "id"),
        all = TRUE,
        sort = FALSE
      ) %>%
      mutate(
        outlier = ifelse(is.na(outlier), FALSE, outlier),
        !!trait := ifelse(
          test = outlier == TRUE,
          yes = NA,
          no = .data[[trait]]
        )
      ) %>%
      droplevels() %>%
      data.frame() %>%
      rename(outlier = outlier)

    td <- createTD(
      data = data_clean,
      genotype = results$inputs$genotype,
      trial = results$inputs$trial,
      repId = results$inputs$rep,
      subBlock = results$inputs$block,
      rowCoord = results$inputs$row,
      colCoord = results$inputs$col,
      trDesign = design_td
    )
    m_models <- fitTD(
      TD = td,
      traits = trait,
      what = c("fixed", "random"),
      spatial = spatial,
      progress = progress,
      engine = engine
    )
  }
  # Heritability
  h2_cullis <- extractSTA(
    STA = m_models,
    what = "heritability"
  )
  # VarComps
  var_gen <- extractSTA(
    STA = m_models,
    what = "varGen"
  )
  var_error <- extractSTA(
    STA = m_models,
    what = "varErr"
  )
  # CV
  coef_var <- extractSTA(
    STA = m_models,
    what = "CV"
  )
  names(h2_cullis)[2] <- "heritability"
  names(var_gen)[2] <- "VarGen"
  names(var_error)[2] <- "VarErr"
  names(coef_var)[2] <- "CV"
  # Summary
  resum_fitted_model <- merge(
    x = merge(x = h2_cullis, coef_var, by = "trial"),
    y = merge(x = var_gen, var_error, by = "trial"),
    by = "trial"
  ) %>%
    mutate(design = design)
  # Fitted Models
  fitted_models <- m_models
  # BLUES
  blues_td <- STAtoTD(m_models, keep = c("trial"), addWt = TRUE)
  blues_td <- data.table::rbindlist(blues_td)
  names(blues_td) <- c(
    "genotype", "trial", "BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "wt"
  )
  blues_blups <- blues_td
  # standardized residuals
  std_res <- extractSTA(
    STA = m_models,
    what = "stdResR"
  )
  std_residuals <- std_res
  out <- list(
    fitted_models = fitted_models,
    resum_fitted_model = resum_fitted_model,
    outliers = outliers_td,
    blues_blups = blues_blups,
    std_residuals = std_residuals
  )
  return(out)
}

#' Single Trial Analysis
#'
#' @param results Object of class \code{checkAgri} resulting of executing
#' \code{check_design_met()} function.
#' @param progress Should the progress of the modeling be printed.
#' If \code{TRUE}, for every trial a line is output indicating the traits fitted
#' for the particular trial.
#' @param engine A character string specifying the name of the mixed modeling
#' engine to use, either \code{lme4} or \code{asreml}. For spatial designs,
#' \code{SpATS} is always used, for other designs \code{asreml} as a default.
#' @param remove_outliers Should outliers be removed? \code{TRUE} by default.
#'
#' @return An object of class \code{smaAgri}, with a list of:
#' \item{fitted_models}{A list containing the fitted models. (Both models, the
#' one with Genotype as Random and the one with Genotype as Fixed)}
#' \item{resum_fitted_model}{A data.frame containing a summary of the fitted
#' models.}
#' \item{outliers}{A data.frame containing extreme observations. If
#' \code{remove_outliers} is \code{TRUE}, this data.frame will contain the
#' observations removed.}
#' \item{blues_blups}{A data.frame containing BLUPs/BLUEs for all the genotypes
#' in each trial.}
#' \item{std_residuals}{A data.frame containing the standardized residuals for
#' the model with genotype as random component.}
#' \item{data}{A data.frame containing the data used. If \code{remove_outliers}
#' is \code{TRUE}, data will have missing values for the outliers detected.}
#' @export
#'
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
#' @importFrom statgenSTA createTD fitTD outlierSTA extractSTA STAtoTD
single_trial_analysis <- function(results = NULL,
                                  progress = TRUE,
                                  engine = "asreml",
                                  remove_outliers = TRUE) {
  if (!inherits(results, "checkAgri")) {
    stop("The object should be of checkAgri class")
  }
  traits <- results$inputs$traits

  fitted_models <- list()
  resum_fitted_model <- list()
  outliers <- list()
  blues_blups <- list()
  std_residuals <- list()

  for (i in traits) {

    # Spatials - res_row_col ---------------------------------------------------

    if ("res_row_col" %in% results$exp_design_list$exp_design) {
      objt_res_row_col <- fit_STA(
        results = results,
        trait = i,
        design = "res_row_col",
        remove_outliers = remove_outliers,
        engine = "SpATS",
        progress = progress
      )
      fitted_models[[i]] <- c(
        fitted_models[[i]],
        objt_res_row_col$fitted_models
      )
      resum_fitted_model[[i]] <- data.table::rbindlist(
        l = list(
          resum_fitted_model[[i]],
          objt_res_row_col$resum_fitted_model
        )
      )
      outliers[[i]] <- data.table::rbindlist(
        l = list(
          outliers[[i]],
          objt_res_row_col$outliers
        )
      )
      blues_blups[[i]] <- data.table::rbindlist(
        l = list(
          blues_blups[[i]],
          objt_res_row_col$blues_blups
        )
      )
      std_residuals[[i]] <- dplyr::bind_rows(
        std_residuals[[i]],
        objt_res_row_col$std_residuals
      )
    }

    # Spatials - row_col -------------------------------------------------------

    if ("row_col" %in% results$exp_design_list$exp_design) {
      objt_row_col <- fit_STA(
        results = results,
        trait = i,
        design = "row_col",
        remove_outliers = remove_outliers,
        engine = "SpATS",
        progress = progress
      )
      fitted_models[[i]] <- c(
        fitted_models[[i]],
        objt_row_col$fitted_models
      )
      resum_fitted_model[[i]] <- data.table::rbindlist(
        l = list(
          resum_fitted_model[[i]],
          objt_row_col$resum_fitted_model
        )
      )
      outliers[[i]] <- data.table::rbindlist(
        l = list(
          outliers[[i]],
          objt_row_col$outliers
        )
      )
      blues_blups[[i]] <- data.table::rbindlist(
        l = list(
          blues_blups[[i]],
          objt_row_col$blues_blups
        )
      )
      std_residuals[[i]] <- dplyr::bind_rows(
        std_residuals[[i]],
        objt_row_col$std_residuals
      )
    }

    # Alpha lattices ---------------------------------------------------------

    if ("alpha_lattice" %in% results$exp_design_list$exp_design) {
      objt_alpha <- fit_STA(
        results = results,
        trait = i,
        design = "alpha_lattice",
        remove_outliers = remove_outliers,
        engine = engine,
        progress = progress
      )
      fitted_models[[i]] <- c(
        fitted_models[[i]],
        objt_alpha$fitted_models
      )
      resum_fitted_model[[i]] <- data.table::rbindlist(
        l = list(
          resum_fitted_model[[i]],
          objt_alpha$resum_fitted_model
        )
      )
      outliers[[i]] <- data.table::rbindlist(
        l = list(
          outliers[[i]],
          objt_alpha$outliers
        )
      )
      blues_blups[[i]] <- data.table::rbindlist(
        l = list(
          blues_blups[[i]],
          objt_alpha$blues_blups
        )
      )
      std_residuals[[i]] <- dplyr::bind_rows(
        std_residuals[[i]],
        objt_alpha$std_residuals
      )
    }

    # RCBD -------------------------------------------------------------------

    if ("rcbd" %in% results$exp_design_list$exp_design) {
      objt_rcbd <- fit_STA(
        results = results,
        trait = i,
        design = "rcbd",
        remove_outliers = remove_outliers,
        engine = engine,
        progress = progress
      )
      fitted_models[[i]] <- c(
        fitted_models[[i]],
        objt_rcbd$fitted_models
      )
      resum_fitted_model[[i]] <- data.table::rbindlist(
        l = list(
          resum_fitted_model[[i]],
          objt_rcbd$resum_fitted_model
        )
      )
      outliers[[i]] <- data.table::rbindlist(
        l = list(
          outliers[[i]],
          objt_rcbd$outliers
        )
      )
      blues_blups[[i]] <- data.table::rbindlist(
        l = list(
          blues_blups[[i]],
          objt_rcbd$blues_blups
        )
      )
      std_residuals[[i]] <- dplyr::bind_rows(
        std_residuals[[i]],
        objt_rcbd$std_residuals
      )
    }

    # CRD ---------------------------------------------------------------------

    if ("crd" %in% results$exp_design_list$exp_design) {
      m_models_crd <- list()

      exp_to_remove <- results$filter[[i]]$trials_to_remove

      crd_trials <- results$data_design %>%
        filter(
          exp_design == "crd" &
            !.data[[results$inputs$trial]] %in% exp_to_remove
        ) %>%
        pull(results$inputs$trial) %>%
        as.character() %>%
        unique()

      data_crd <- results$data_design %>%
        filter(.data[[results$inputs$trial]] %in% crd_trials) %>%
        droplevels()

      if (nrow(data_crd) > 0) {
        td_crd <- fit_crd(
          data = data_crd,
          trial = results$inputs$trial,
          genotype = results$inputs$genotype,
          response = i
        )
        m_models_crd <- list(
          mRand = td_crd$models_rand,
          mFix = td_crd$models_fixed
        )

        # Residuals Cleaning
        if (!is.null(td_crd$outliers) && remove_outliers) {
          outliers_crd <- td_crd %>%
            .[["outliers"]] %>%
            dplyr::select(trial, genotype, id, outlier)

          outliers[[i]] <- data.table::rbindlist(
            l = list(
              outliers[[i]],
              outliers_crd
            )
          )
          data_crd_clean <- data_crd %>%
            merge(
              x = .,
              y = outliers_crd,
              by.x = c(results$inputs$trial, results$inputs$genotype, "id"),
              by.y = c("trial", "genotype", "id"),
              all = TRUE,
              sort = FALSE
            ) %>%
            mutate(
              outlier = ifelse(is.na(outlier), FALSE, outlier),
              !!i := ifelse(
                test = outlier == TRUE,
                yes = NA,
                no = .data[[i]]
              )
            ) %>%
            droplevels() %>%
            data.frame() %>%
            rename(outlier_crd = outlier)
          td_crd <- fit_crd(
            data = data_crd_clean,
            trial = results$inputs$trial,
            genotype = results$inputs$genotype,
            response = i
          )
          m_models_crd <- list(
            mRand = td_crd$models_rand,
            mFix = td_crd$models_fixed
          )
        }
        fitted_models[[i]] <- c(
          fitted_models[[i]],
          m_models_crd
        )
        # summary fitted models
        resum_fitted_model_crd <- td_crd$resum_fitted_model
        resum_fitted_model[[i]] <- data.table::rbindlist(
          l = list(
            resum_fitted_model[[i]],
            resum_fitted_model_crd
          )
        )
        # BLUES
        blues_td_crd <- td_crd$blues_blups
        blues_blups[[i]] <- data.table::rbindlist(
          l = list(
            blues_blups[[i]],
            blues_td_crd
          )
        )
        # standardized residuals
        std_res_crd <- td_crd$residuals
        std_residuals[[i]] <- dplyr::bind_rows(
          std_residuals[[i]],
          std_res_crd
        )
      }
    }
  }
  # stacking tables
  blues_blups <- data.table::rbindlist(blues_blups, idcol = "trait")
  resum_fitted_model <- data.table::rbindlist(
    l = resum_fitted_model,
    idcol = "trait"
  )
  std_residuals <- data.table::rbindlist(
    l = std_residuals,
    idcol = "trait",
    fill = TRUE
  )
  outliers <- data.table::rbindlist(outliers, idcol = "trait")
  # data used
  if (nrow(outliers) != 0 && remove_outliers) {
    traits_out <- outliers %>%
      dplyr::pull(trait) %>%
      unique()
    data_objt <- results$data_design
    for (k in traits_out) {
      indx_out <- outliers %>%
        dplyr::filter(trait %in% k) %>%
        dplyr::pull(id)
      data_objt[data_objt$id %in% indx_out, k] <- NA
      data_objt[data_objt$id %in% indx_out, "outliers"] <- TRUE
    }
  } else {
    data_objt <- results$data_design
  }
  # Output
  results <- list(
    fitted_models = fitted_models,
    resum_fitted_model = resum_fitted_model,
    outliers = outliers,
    blues_blups = blues_blups,
    std_residuals = std_residuals,
    data = data_objt
  )
  class(results) <- "smaAgri"
  return(invisible(results))
}
