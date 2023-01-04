#' Stability Coefficients
#'
#' @param predictions A data.frame with one value per GxE combination.
#' @param genotype A character string indicating the column in data that
#' contains genotypes.
#' @param trial A character string indicating the column in data that contains
#' trials.
#' @param response A character vector specifying the response variable.
#' @param best A character string specifying how to define the best genotype
#' by numeric value ("min", "max"). "max" by default.
#'
#' @return A data.frame with several stability measures.
#' "superiority" (cultivar-superiority measure), "static" (Shukla's stability
#' variance) and "wricke" (Wricke's ecovalence)
#' @export
#'
#' @examples
#' # in progress
stability <- function(predictions = NULL,
                      genotype = NULL,
                      trial = NULL,
                      response = NULL,
                      best = "max") {
  names_env <- predictions %>%
    type.convert(as.is = FALSE) %>%
    pull(.data[[trial]]) %>%
    levels()
  n_env <- length(names_env)
  E <- mean(predictions$predicted.value, na.rm = TRUE)
  best_var <- ifelse(test = best == "max", yes = "max_trial", no = "min_trial")
  predictions %>%
    group_by(.data[[trial]]) %>%
    mutate(
      max_trial = max(predicted.value, na.rm = TRUE),
      min_trial = min(predicted.value, na.rm = TRUE),
      mean_trial = mean(predicted.value, na.rm = TRUE)
    ) %>%
    group_by(.data[[genotype]]) %>%
    mutate(mean_gen = mean(predicted.value, na.rm = TRUE)) %>%
    summarise(
      superiority = sum(
        (predicted.value - .data[[best_var]])^2
      ) / (2 * n_env) %>% sqrt(),
      static = sum((predicted.value - mean_gen)^2) / (n_env - 1) %>% sqrt(),
      wricke = sum((predicted.value - mean_gen - mean_trial + E)^2) %>% sqrt(),
      predicted.value = mean(predicted.value, na.rm = TRUE)
    ) %>%
    arrange(predicted.value) %>%
    as.data.frame()
}


#' Multi-Environmental Trial Analysis
#'
#' @param sma_output Object of class 'smaAgri' resulting of executing
#' single_model_analysis function.
#' @param h2_filter Numeric value to filter trials with poor heritability.
#' 0.2 by default.
#' @param workspace Sets the workspace for the core REML routines in the form of
#' a number optionally followed directly by a valid measurement unit. "128mb" by
#' default.
#' @param trials_to_fit_fa Number of trials necessary to fit a Factor Analytic
#' structure for the GxE interaction term. 4 by default.
#'
#' @return A list with a summary of the fitted models.
#' @export
#'
#' @examples
#' # In progress
met_analysis <- function(sma_output = NULL,
                         h2_filter = 0.2,
                         workspace = "1gb",
                         trials_to_fit_fa = 4) {
  if (!inherits(sma_output, "smaAgri")) {
    stop("The object should be of smaAgri class")
  }
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop("The package asreml is not loaded.")
  }
  asreml::asreml.options(trace = FALSE, workspace = workspace)
  data_td <- sma_output$blues_blups
  traits <- data_td %>%
    dplyr::pull("trait") %>%
    unique() %>%
    as.character()
  trials <- data_td %>%
    dplyr::pull("trial") %>%
    unique() %>%
    as.character()
  n_trials <- length(trials)
  if (n_trials <= 1) {
    stop("There is only one trial to fit an MET model.")
  }
  conn <- connectivity_matrix(
    data = data_td,
    genotype = "genotype",
    trial = "trial"
  )
  ceros <- which(conn == 0, arr.ind = TRUE)
  if (nrow(ceros) > 1) {
    warning(
      "Some trials have zero connectivity: \t",
      paste(rownames(ceros), collapse = ", ")
    )
  }
  minimun_req <- which(conn < 20, arr.ind = TRUE)
  if (nrow(minimun_req) > 1) {
    message(
      "Some trials could have poor connectivity: \t",
      paste(rownames(ceros), collapse = ", ")
    )
  }
  met_models <- VCOV <- trial_effects <- overall_BLUPs <- BLUPs_GxE <- list()
  stab_list <- h2_list <- list()

  for (var in traits) {
    trials_to_keep <- sma_output$resum_fitted_model %>%
      dplyr::filter(trait %in% var & heritability > h2_filter) %>%
      droplevels() %>%
      dplyr::pull(trial) %>%
      as.character()
    dt <- data_td %>%
      dplyr::filter(trait %in% var & trial %in% trials_to_keep) %>%
      droplevels() %>%
      as.data.frame()
    equation_fix <- stats::reformulate("trial", response = "BLUEs")
    if (n_trials < trials_to_fit_fa) {
      vcov_selected <- "us"
      equation_ran <- stats::reformulate("us(trial):genotype")
    }
    if (n_trials >= trials_to_fit_fa) {
      vcov_selected <- "fa2"
      equation_ran <- stats::reformulate("fa(trial, 2):genotype")
    }
    met_mod <- suppressWarnings(
      asreml::asreml(
        fixed = equation_fix,
        random = equation_ran,
        data = dt,
        weights = wt,
        family = asreml::asr_gaussian(dispersion = 1),
        na.action = list(x = "include", y = "include"),
        trace = 0,
        maxiter = 200
      )
    )
    met_mod <- suppressWarnings(asreml::update.asreml(met_mod))
    met_models[[var]] <- met_mod
    VCOV[[var]] <- extractG(
      model = met_mod,
      gen = "genotype",
      env = "trial",
      vc.model = vcov_selected
    )
    # trial effects
    trial_effects[[var]] <- rbind(
      trial_effects[[var]],
      suppressWarnings(
        as.data.frame(
          asreml::predict.asreml(met_mod, classify = "trial")$pvals
        )
      )
    )
    # Overall BLUPs
    BLUPs <- suppressWarnings(
      asreml::predict.asreml(met_mod, classify = "genotype", sed = TRUE)
    )
    overall_BLUPs[[var]] <- rbind(
      overall_BLUPs[[var]],
      suppressWarnings(
        as.data.frame(
          BLUPs$pvals
        )
      )
    )
    # BLUPs GxE
    BLUPs_GxE[[var]] <- rbind(
      BLUPs_GxE[[var]],
      as.data.frame(
        suppressWarnings(
          asreml::predict.asreml(met_mod, classify = "genotype:trial")$pvals
        )
      )
    )
    # Stability
    stab_genotypes <- stability(
      predictions = BLUPs_GxE[[var]],
      trial = "trial",
      genotype = "genotype",
      best = "max"
    )
    stab_list[[var]] <- rbind(stab_list[[var]], stab_genotypes)
    # Heritability
    vcov_mat <- VCOV[[var]]$VCOV
    g_var <- mean(vcov_mat[upper.tri(vcov_mat, diag = FALSE)])
    vdBLUP.mat <- BLUPs$sed^2
    vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag = FALSE)])
    h2_tmp <- 1 - (vdBLUP.avg / 2 / g_var)
    h2_tmp <- data.frame(trait = var, h2 = h2_tmp)
    h2_list[[var]] <- rbind(h2_list[[var]], h2_tmp)
  }
  trial_effects <- dplyr::bind_rows(trial_effects, .id = "trait")
  overall_BLUPs <- dplyr::bind_rows(overall_BLUPs, .id = "trait")
  BLUPs_GxE <- dplyr::bind_rows(BLUPs_GxE, .id = "trait")
  stab_list <- dplyr::bind_rows(stab_list, .id = "trait")
  h2_list <- dplyr::bind_rows(h2_list, .id = "trait")
  objt_out <- list(
    trial_effects = trial_effects,
    overall_BLUPs = overall_BLUPs,
    BLUPs_GxE = BLUPs_GxE,
    VCOV = VCOV,
    stability = stab_list,
    heritability = h2_list
  )
  return(objt_out)
}
