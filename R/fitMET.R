#' Stability Coefficients
#'
#' @param predictions A data.frame with one value per GxE combination.
#' @param genotype A character string indicating the column in predictions that
#' contains genotypes.
#' @param trial A character string indicating the column in predictions that
#' contains trials.
#' @param response A character string specifying the response variable.
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
      superiority = sqrt(sum(
        (predicted.value - .data[[best_var]])^2
      ) / (2 * n_env)),
      static = sqrt(sum((predicted.value - mean_gen)^2) / (n_env - 1)),
      wricke = sqrt(sum((predicted.value - mean_gen - mean_trial + E)^2)),
      predicted.value = mean(predicted.value, na.rm = TRUE)
    ) %>%
    arrange(predicted.value) %>%
    as.data.frame()
}


#' Multi-Environmental Trial Analysis
#'
#' @param sma_output Object of class \code{smaAgri} resulting of executing
#' \code{single_trial_analysis()} function.
#' @param h2_filter Numeric value to filter trials with poor heritability.
#' 0.2 by default.
#' @param workspace Sets the workspace for the core \code{REML} routines in the
#' form of a number optionally followed directly by a valid measurement unit.
#' "128mb" by default.
#' @param vcov A character string specifying the Variance-Covariance structure
#' to be fitted. Can be "fa2", "fa1", "us", or "corh". If \code{NULL} the
#' function will try to fit an "us" Variance-Covariance and if it fails, it will
#' try with "fa2" and then with "fa1".
#'
#' @return  An object of class \code{metAgri}, with a list of:
#' \item{trial_effects}{A data.frame containing Trial BLUEs.}
#' \item{overall_BLUPs}{A data.frame containing Genotypic BLUPs across trials,
#' by trait.}
#' \item{BLUPs_GxE}{A data.frame containing Genotypic BLUPs by trial/trait.}
#' \item{VCOV}{A list by trait contanining the variance-covariance fitted.}
#' \item{stability}{A data.frame containing several Stability coefficients
#' resulting of executing the function \code{stability()}.}
#' \item{heritability}{A data.frame containing overall heritabilities by trait.}
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
#' out <- single_trial_analysis(results, progress = FALSE)
#' met_results <- met_analysis(out)
#' print(met_results)
#' covcor_heat(matrix = met_results$VCOV$yield$CORR)
#' }
met_analysis <- function(sma_output = NULL,
                         h2_filter = 0.2,
                         workspace = "1gb",
                         vcov = NULL) {
  if (!inherits(sma_output, "smaAgri")) {
    stop("The object should be of smaAgri class")
  }
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop("The package asreml is not loaded.")
  }
  asreml::asreml.options(trace = FALSE, workspace = workspace)
  met_models <- VCOV <- trial_effects <- overall_BLUPs <- BLUPs_GxE <- list()
  stab_list <- h2_list <- list()

  data_td <- sma_output$blues_blups %>%
    dplyr::mutate(
      trial = as.factor(trial),
      genotype = as.factor(genotype)
    )

  traits <- data_td %>%
    dplyr::pull("trait") %>%
    unique() %>%
    as.character()

  for (var in traits) {

    trials_to_keep <- sma_output$resum_fitted_model %>%
      dplyr::filter(heritability > h2_filter & trait %in% var) %>%
      droplevels() %>%
      dplyr::pull(trial) %>%
      as.character()

    n_trials <- length(trials_to_keep)
    if (n_trials <= 1) {
      stop("There is only one trial to fit an MET model in '", var, "'")
    }

    dt <- data_td %>%
      dplyr::filter(trait %in% var & trial %in% trials_to_keep) %>%
      droplevels() %>%
      as.data.frame()

    conn <- connectivity_matrix(
      data = data_td,
      genotype = "genotype",
      trial = "trial",
      response = "BLUEs"
    )
    ceros <- which(conn == 0, arr.ind = TRUE)
    if (nrow(ceros) > 1) {
      warning(
        "Some trials have zero connectivity: \n",
        paste(unique(rownames(ceros)), collapse = ", "),
        "\n"
      )
    }
    minimun_req <- which(conn < 20, arr.ind = TRUE)
    if (nrow(minimun_req) > 1) {
      warning(
        "Some trials could have poor connectivity: \n",
        paste(unique(rownames(ceros)), collapse = ", ")
      )
    }

    equation_fix <- stats::reformulate("trial", response = "BLUEs")
    if (is.null(vcov)) {
      vcov_selected <- "us"
      equation_ran <- stats::reformulate("us(trial):genotype")
      met_mod <- try(
        suppressWarnings(
          asreml::asreml(
            fixed = equation_fix,
            random = equation_ran,
            data = dt,
            weights = wt,
            family = asreml::asr_gaussian(dispersion = 1),
            na.action = list(x = "exclude", y = "include"),
            trace = 0,
            maxiter = 200
          )
        ),
      )
      if (inherits(met_mod, "try-error")) {
        vcov_selected <- "fa2"
        equation_ran <- stats::reformulate("fa(trial, 2):genotype")
        met_mod <- try(
          suppressWarnings(
            asreml::asreml(
              fixed = equation_fix,
              random = equation_ran,
              data = dt,
              weights = wt,
              family = asreml::asr_gaussian(dispersion = 1),
              na.action = list(x = "exclude", y = "include"),
              trace = 0,
              maxiter = 200
            )
          ),
        )
      }
      if (inherits(met_mod, "try-error")) {
        vcov_selected <- "fa1"
        equation_ran <- stats::reformulate("fa(trial, 1):genotype")
        met_mod <- try(
          suppressWarnings(
            asreml::asreml(
              fixed = equation_fix,
              random = equation_ran,
              data = dt,
              weights = wt,
              family = asreml::asr_gaussian(dispersion = 1),
              na.action = list(x = "exclude", y = "include"),
              trace = 0,
              maxiter = 200
            )
          ),
        )
      }
      if (inherits(met_mod, "try-error")) {
        stop(
          "Trait '", var, "'\n",
          "We couldn't fit any variance-covariance structure."
        )
      }
    } else if (vcov == "fa1") {
      vcov_selected <- "fa1"
      equation_ran <- stats::reformulate("fa(trial, 1):genotype")
    } else if (vcov == "fa2") {
      vcov_selected <- "fa2"
      equation_ran <- stats::reformulate("fa(trial, 2):genotype")
    } else if (vcov == "us") {
      vcov_selected <- "us"
      equation_ran <- stats::reformulate("us(trial):genotype")
    } else if (vcov == "corh") {
      vcov_selected <- "corh"
      equation_ran <- stats::reformulate("corh(trial):genotype")
    } else if (vcov == "corv") {
      vcov_selected <- "corv"
      equation_ran <- stats::reformulate("corv(trial):genotype")
    } else {
      stop(paste0("No '", vcov, "' variance-covariance structure found."))
    }
    met_mod <- try(
      suppressWarnings(
        asreml::asreml(
          fixed = equation_fix,
          random = equation_ran,
          data = dt,
          weights = wt,
          family = asreml::asr_gaussian(dispersion = 1),
          na.action = list(x = "exclude", y = "include"),
          trace = 0,
          maxiter = 200
        )
      )
    )
    if (inherits(met_mod, "try-error")) {
      stop(
        "Trait '", var, "'\n",
        "We couldn't fit any variance-covariance structure."
      )
    }
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
  class(objt_out) <- "metAgri"
  return(invisible(objt_out))
}
