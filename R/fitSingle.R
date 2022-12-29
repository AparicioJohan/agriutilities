#' Single Trial Analysis
#'
#' @param results Object of class 'checkAgri' resulting of executing
#' check_design_MET function.
#' @param progress Should the progress of the modeling be printed.
#' If TRUE, for every trial a line is output indicating the traits fitted for
#' the particular trial.
#' @param engine A character string specifying the name of the mixed modeling
#' engine to use, either "lme4" or "asreml". For spatial designs, "SpaTS" is
#' always used, for other designs "asreml" as a default.
#' @param remove_outliers Should outliers be removed? TRUE by default.
#'
#' @return A list of data.frames.
#' @export
#'
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
#' out <- single_model_analysis(results, progress = FALSE)
#' out$resum_fitted_model
#' }
#' @importFrom statgenSTA createTD fitTD outlierSTA extractSTA STAtoTD
single_model_analysis <- function(results = NULL,
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
      m_models_res_row_col <- list()

      exp_to_remove <- results$filter[[i]]$trials_to_remove

      res_row_col_trials <- results$data_design %>%
        filter(
          exp_design == "res_row_col" &
            !.data[[results$inputs$trial]] %in% exp_to_remove
        ) %>%
        pull(results$inputs$trial) %>%
        as.character() %>%
        unique()

      data_res_row_col <- results$data_design %>%
        filter(.data[[results$inputs$trial]] %in% res_row_col_trials) %>%
        droplevels()

      if (nrow(data_res_row_col) > 0) {
        td_res_row_col <- createTD(
          data = data_res_row_col,
          genotype = results$inputs$genotype,
          trial = results$inputs$trial,
          repId = results$inputs$rep,
          rowCoord = results$inputs$row,
          colCoord = results$inputs$col,
          trDesign = "res.rowcol"
        )
        m_models_res_row_col <- fitTD(
          TD = td_res_row_col,
          traits = i,
          what = c("fixed", "random"),
          spatial = TRUE,
          progress = progress,
          engine = "SpATS"
        )

        # Residuals
        outliers_res_row_col <- outlierSTA(
          STA = m_models_res_row_col,
          traits = i,
          rLimit = 3,
          verbose = FALSE
        )

        # Cleaning
        if (!is.null(outliers_res_row_col$outliers) && remove_outliers) {
          outliers_res_row_col <- outliers_res_row_col %>%
            .[["outliers"]] %>%
            dplyr::select(trial, genotype, id, outlier)

          outliers[[i]] <- rbind(outliers[[i]], outliers_res_row_col)

          data_res_row_col_clean <- data_res_row_col %>%
            merge(
              x = .,
              y = outliers_res_row_col,
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
            rename(outlier_res_row_col = outlier)

          td_res_row_col <- createTD(
            data = data_res_row_col_clean,
            genotype = results$inputs$genotype,
            trial = results$inputs$trial,
            repId = results$inputs$rep,
            rowCoord = results$inputs$row,
            colCoord = results$inputs$col,
            trDesign = "res.rowcol"
          )
          m_models_res_row_col <- fitTD(
            TD = td_res_row_col,
            traits = i,
            what = c("fixed", "random"),
            spatial = TRUE,
            progress = progress,
            engine = "SpATS"
          )
        }
        # Heritability
        h2_cullis_res_row_col <- extractSTA(
          STA = m_models_res_row_col,
          what = "heritability"
        )
        # VarComps
        VarG_res_row_col <- extractSTA(
          STA = m_models_res_row_col,
          what = "varGen"
        )
        VarE_res_row_col <- extractSTA(
          STA = m_models_res_row_col,
          what = "varErr"
        )
        # CV
        CV_res_row_col <- extractSTA(
          STA = m_models_res_row_col,
          what = "CV"
        )

        names(h2_cullis_res_row_col)[2] <- "heritability"
        names(VarG_res_row_col)[2] <- "VarGen"
        names(VarE_res_row_col)[2] <- "VarErr"
        names(CV_res_row_col)[2] <- "CV"

        resum_fitted_model_res_row_col <- merge(
          x = merge(x = h2_cullis_res_row_col, CV_res_row_col, by = "trial"),
          y = merge(x = VarG_res_row_col, VarE_res_row_col, by = "trial"),
          by = "trial"
        ) %>%
          mutate(design = "res_row_col")

        fitted_models[[i]] <- m_models_res_row_col

        resum_fitted_model[[i]] <- rbind(
          resum_fitted_model[[i]],
          resum_fitted_model_res_row_col
        )

        # BLUES
        blues_TD_res_row_col <- STAtoTD(m_models_res_row_col,
          keep = c("trial"),
          addWt = TRUE
        )
        blues_TD_res_row_col <- do.call(rbind, lapply(blues_TD_res_row_col, as.data.frame))
        blues_blups[[i]] <- blues_TD_res_row_col

        # standardized residuals
        stdRes_res_row_col <- extractSTA(
          STA = m_models_res_row_col,
          what = "stdResR"
        )
        std_residuals[[i]] <- stdRes_res_row_col
      }
    }

    # Spatials - row_col -------------------------------------------------------

    if ("row_col" %in% results$exp_design_list$exp_design) {
      m_models_row_col <- list()

      exp_to_remove <- results$filter[[i]]$trials_to_remove

      row_col_trials <- results$data_design %>%
        filter(
          exp_design == "row_col" &
            !.data[[results$inputs$trial]] %in% exp_to_remove
        ) %>%
        pull(results$inputs$trial) %>%
        as.character() %>%
        unique()

      data_row_col <- results$data_design %>%
        filter(.data[[results$inputs$trial]] %in% row_col_trials) %>%
        droplevels()

      if (nrow(data_row_col) > 0) {
        td_row_col <- createTD(
          data = data_row_col,
          genotype = results$inputs$genotype,
          trial = results$inputs$trial,
          rowCoord = results$inputs$row,
          colCoord = results$inputs$col,
          trDesign = "rowcol"
        )
        m_models_row_col <- fitTD(
          TD = td_row_col,
          traits = i,
          what = c("fixed", "random"),
          spatial = TRUE,
          progress = progress,
          engine = "SpATS"
        )

        # Residuals
        outliers_row_col <- outlierSTA(
          STA = m_models_row_col,
          traits = i,
          rLimit = 3,
          verbose = FALSE
        )

        # Cleaning
        if (!is.null(outliers_row_col$outliers) && remove_outliers) {
          outliers_row_col <- outliers_row_col %>%
            .[["outliers"]] %>%
            dplyr::select(trial, genotype, id, outlier)

          outliers[[i]] <- rbind(outliers[[i]], outliers_row_col)

          data_row_col_clean <- data_row_col %>%
            merge(
              x = .,
              y = outliers_row_col,
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
            rename(outlier_row_col = outlier)

          td_row_col <- createTD(
            data = data_row_col_clean,
            genotype = results$inputs$genotype,
            trial = results$inputs$trial,
            repId = results$inputs$rep,
            rowCoord = results$inputs$row,
            colCoord = results$inputs$col,
            trDesign = "rowcol"
          )
          m_models_row_col <- fitTD(
            TD = td_row_col,
            traits = i,
            what = c("fixed", "random"),
            spatial = TRUE,
            progress = progress,
            engine = "SpATS"
          )
        }
        # Heritability
        h2_cullis_row_col <- extractSTA(
          STA = m_models_row_col,
          what = "heritability"
        )
        # VarComps
        VarG_row_col <- extractSTA(
          STA = m_models_row_col,
          what = "varGen"
        )
        VarE_row_col <- extractSTA(
          STA = m_models_row_col,
          what = "varErr"
        )
        # CV
        CV_row_col <- extractSTA(
          STA = m_models_row_col,
          what = "CV"
        )

        names(h2_cullis_row_col)[2] <- "heritability"
        names(VarG_row_col)[2] <- "VarGen"
        names(VarE_row_col)[2] <- "VarErr"
        names(CV_row_col)[2] <- "CV"

        resum_fitted_model_row_col <- merge(
          x = merge(x = h2_cullis_row_col, CV_row_col, by = "trial"),
          y = merge(x = VarG_row_col, VarE_row_col, by = "trial"),
          by = "trial"
        ) %>%
          mutate(design = "row_col")

        fitted_models[[i]] <- m_models_row_col

        resum_fitted_model[[i]] <- rbind(
          resum_fitted_model[[i]],
          resum_fitted_model_row_col
        )

        # BLUES
        blues_TD_row_col <- STAtoTD(m_models_row_col,
          keep = c("trial"),
          addWt = TRUE
        )
        blues_TD_row_col <- do.call(rbind, lapply(blues_TD_row_col, as.data.frame))
        blues_blups[[i]] <- blues_TD_row_col

        # standardized residuals
        stdRes_row_col <- extractSTA(
          STA = m_models_row_col,
          what = "stdResR"
        )
        std_residuals[[i]] <- stdRes_row_col
      }
    }

    # Alpha lattices ---------------------------------------------------------

    if ("alpha_lattice" %in% results$exp_design_list$exp_design) {
      m_models_alpha <- list()

      exp_to_remove <- results$filter[[i]]$trials_to_remove

      alpha_trials <- results$data_design %>%
        filter(
          exp_design == "alpha_lattice" &
            !.data[[results$inputs$trial]] %in% exp_to_remove
        ) %>%
        pull(results$inputs$trial) %>%
        as.character() %>%
        unique()

      data_alpha <- results$data_design %>%
        filter(.data[[results$inputs$trial]] %in% alpha_trials) %>%
        droplevels()

      if (nrow(data_alpha) > 0) {
        td_alpha <- createTD(
          data = data_alpha,
          genotype = results$inputs$genotype,
          trial = results$inputs$trial,
          repId = results$inputs$rep,
          subBlock = results$inputs$block,
          trDesign = "res.ibd"
        )
        m_models_alpha <- fitTD(
          TD = td_alpha,
          traits = i,
          what = c("fixed", "random"),
          spatial = FALSE,
          progress = progress,
          engine = engine
        )

        # Residuals
        outliers_alpha <- outlierSTA(
          STA = m_models_alpha,
          traits = i,
          rLimit = 3,
          verbose = FALSE
        )

        # Cleaning
        if (!is.null(outliers_alpha$outliers) && remove_outliers) {
          outliers_alpha <- outliers_alpha %>%
            .[["outliers"]] %>%
            dplyr::select(trial, genotype, id, outlier)

          outliers[[i]] <- rbind(outliers[[i]], outliers_alpha)

          data_alpha_clean <- data_alpha %>%
            merge(
              x = .,
              y = outliers_alpha,
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
            rename(outlier_alpha = outlier)

          td_alpha <- createTD(
            data = data_alpha_clean,
            genotype = results$inputs$genotype,
            trial = results$inputs$trial,
            repId = results$inputs$rep,
            subBlock = results$inputs$block,
            trDesign = "res.ibd"
          )
          m_models_alpha <- fitTD(
            TD = td_alpha,
            traits = i,
            what = c("fixed", "random"),
            spatial = FALSE,
            progress = progress,
            engine = engine
          )
        }
        # Heritability
        h2_cullis_alpha <- extractSTA(
          STA = m_models_alpha,
          what = "heritability"
        )
        # VarComps
        VarG_alpha <- extractSTA(
          STA = m_models_alpha,
          what = "varGen"
        )
        VarE_alpha <- extractSTA(
          STA = m_models_alpha,
          what = "varErr"
        )
        # CV
        CV_alpha <- extractSTA(
          STA = m_models_alpha,
          what = "CV"
        )

        names(h2_cullis_alpha)[2] <- "heritability"
        names(VarG_alpha)[2] <- "VarGen"
        names(VarE_alpha)[2] <- "VarErr"
        names(CV_alpha)[2] <- "CV"

        resum_fitted_model_alpha <- merge(
          x = merge(x = h2_cullis_alpha, CV_alpha, by = "trial"),
          y = merge(x = VarG_alpha, VarE_alpha, by = "trial"),
          by = "trial"
        ) %>%
          mutate(design = "alpha_lattice")

        fitted_models[[i]] <- m_models_alpha

        resum_fitted_model[[i]] <- rbind(
          resum_fitted_model[[i]],
          resum_fitted_model_alpha
        )

        # BLUES
        blues_TD_alpha <- STAtoTD(m_models_alpha,
          keep = c("trial"),
          addWt = TRUE
        )
        blues_TD_alpha <- do.call(rbind, lapply(blues_TD_alpha, as.data.frame))
        blues_blups[[i]] <- blues_TD_alpha

        # standardized residuals
        stdRes_alpha <- extractSTA(
          STA = m_models_alpha,
          what = "stdResR"
        )
        std_residuals[[i]] <- stdRes_alpha
      }
    }

    # RCBD -------------------------------------------------------------------

    if ("rcbd" %in% results$exp_design_list$exp_design) {
      m_models_rcbd <- list()

      exp_to_remove <- results$filter[[i]]$trials_to_remove

      rcbd_trials <- results$data_design %>%
        filter(
          exp_design == "rcbd" &
            !.data[[results$inputs$trial]] %in% exp_to_remove
        ) %>%
        pull(results$inputs$trial) %>%
        as.character() %>%
        unique()

      data_rcbd <- results$data_design %>%
        filter(.data[[results$inputs$trial]] %in% rcbd_trials) %>%
        droplevels()

      if (nrow(data_rcbd) > 0) {
        td_rcbd <- createTD(
          data = data_rcbd,
          genotype = results$inputs$genotype,
          trial = results$inputs$trial,
          repId = results$inputs$rep,
          subBlock = results$inputs$block,
          trDesign = "rcbd"
        )
        m_models_rcbd <- fitTD(
          TD = td_rcbd,
          traits = i,
          what = c("fixed", "random"),
          spatial = FALSE,
          progress = progress,
          engine = engine
        )

        # Residuals
        outliers_rcbd <- outlierSTA(
          STA = m_models_rcbd,
          traits = i,
          rLimit = 3,
          verbose = FALSE
        )

        # Cleaning
        if (!is.null(outliers_rcbd$outliers) && remove_outliers) {
          outliers_rcbd <- outliers_rcbd %>%
            .[["outliers"]] %>%
            dplyr::select(trial, genotype, id, outlier)

          outliers[[i]] <- rbind(outliers[[i]], outliers_rcbd)

          data_rcbd_clean <- data_rcbd %>%
            merge(
              x = .,
              y = outliers_rcbd,
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
            rename(outlier_rcbd = outlier)

          td_rcbd <- createTD(
            data = data_rcbd_clean,
            genotype = results$inputs$genotype,
            trial = results$inputs$trial,
            repId = results$inputs$rep,
            subBlock = results$inputs$block,
            trDesign = "rcbd"
          )
          m_models_rcbd <- fitTD(
            TD = td_rcbd,
            traits = i,
            what = c("fixed", "random"),
            spatial = FALSE,
            progress = progress,
            engine = engine
          )
        }
        # Heritability
        h2_cullis_rcbd <- extractSTA(
          STA = m_models_rcbd,
          what = "heritability"
        )
        # VarComps
        VarG_rcbd <- extractSTA(
          STA = m_models_rcbd,
          what = "varGen"
        )
        VarE_rcbd <- extractSTA(
          STA = m_models_rcbd,
          what = "varErr"
        )
        # CV
        CV_rcbd <- extractSTA(
          STA = m_models_rcbd,
          what = "CV"
        )

        names(h2_cullis_rcbd)[2] <- "heritability"
        names(VarG_rcbd)[2] <- "VarGen"
        names(VarE_rcbd)[2] <- "VarErr"
        names(CV_rcbd)[2] <- "CV"

        resum_fitted_model_rcbd <- merge(
          x = merge(x = h2_cullis_rcbd, CV_rcbd, by = "trial"),
          y = merge(x = VarG_rcbd, VarE_rcbd, by = "trial"),
          by = "trial"
        ) %>%
          mutate(design = "rcbd")

        fitted_models[[i]] <- m_models_rcbd

        resum_fitted_model[[i]] <- rbind(
          resum_fitted_model[[i]],
          resum_fitted_model_rcbd
        )

        # BLUES
        blues_TD_rcbd <- STAtoTD(m_models_rcbd,
          keep = c("trial"),
          addWt = TRUE
        )
        blues_TD_rcbd <- do.call(rbind, lapply(blues_TD_rcbd, as.data.frame))
        blues_blups[[i]] <- blues_TD_rcbd

        # standardized residuals
        stdRes_rcbd <- extractSTA(
          STA = m_models_rcbd,
          what = "stdResR"
        )
        std_residuals[[i]] <- stdRes_rcbd
      }
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

          outliers[[i]] <- rbind(outliers[[i]], outliers_crd)

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
            data = data_crd,
            trial = results$inputs$trial,
            genotype = results$inputs$genotype,
            response = i
          )
          m_models_crd <- list(
            mRand = td_crd$models_rand,
            mFix = td_crd$models_fixed
          )
        }

        resum_fitted_model_crd <- td_crd$resum_fitted_model

        fitted_models[[i]] <- m_models_crd

        resum_fitted_model[[i]] <- rbind(
          resum_fitted_model[[i]],
          resum_fitted_model_crd
        )

        # BLUES
        blues_TD_crd <- td_crd$blues_blups
        blues_blups[[i]] <- blues_TD_crd

        # standardized residuals
        stdRes_crd <- td_crd$residuals
        std_residuals[[i]] <- stdRes_crd
      }
    }
  }

  return(
    list(
      fitted_models = fitted_models,
      resum_fitted_model = resum_fitted_model,
      outliers = outliers,
      blues_blups = blues_blups,
      std_residuals = std_residuals
    )
  )
}
