#############################################################
# IRW PERSON-FIT ANALYSES 
#############################################################

require(ltm)
require(irtoys)
require(PerFit)
require(aberrance)

# ============================================================
# 1. GENERAL IRT FUNCTIONS FOR 1PL / 2PL / 3PL
# ============================================================
# irtoys::est() returns item parameters as a, b, c for all three
# models. Under 1PL and 2PL, c = 0; under 1PL, a is common.

Pi_IRT <- function(theta, item_parameters) {

  a <- item_parameters[, 1]
  b <- item_parameters[, 2]
  c <- item_parameters[, 3]

  logistic_part <- plogis(a * (theta - b))
  probability <- c + (1 - c) * logistic_part

  first_derivative <-
    (1 - c) *
    a *
    logistic_part *
    (1 - logistic_part)

  second_derivative <-
    (1 - c) *
    a^2 *
    logistic_part *
    (1 - logistic_part) *
    (1 - 2 * logistic_part)

  probability <- pmin(
    pmax(probability, 1e-10),
    1 - 1e-10
  )

  list(
    Pi = probability,
    dPi = first_derivative,
    d2Pi = second_derivative
  )
}


ri_IRT <- function(theta, item_parameters) {

  values <- Pi_IRT(theta, item_parameters)

  values$dPi /
    (values$Pi * (1 - values$Pi))
}


r0_WL <- function(theta, item_parameters) {

  values <- Pi_IRT(theta, item_parameters)
  r_values <- ri_IRT(theta, item_parameters)

  denominator <- 2 * sum(r_values * values$dPi)

  if (!is.finite(denominator) || abs(denominator) < 1e-12) {
    return(NA_real_)
  }

  sum(r_values * values$d2Pi) / denominator
}


item_weights <- function(theta, item_parameters, index) {

  index <- match.arg(index, c("U", "W"))

  probability <- Pi_IRT(theta, item_parameters)$Pi
  J <- nrow(item_parameters)

  switch(
    index,
    U = (1 - 2 * probability) /
      (J * probability * (1 - probability)),
    W = (1 - 2 * probability) /
      sum(probability * (1 - probability))
  )
}


Wn <- function(response_vector, theta, item_parameters, index) {

  probability <- Pi_IRT(theta, item_parameters)$Pi

  weights <- item_weights(
    theta = theta,
    item_parameters = item_parameters,
    index = index
  )

  sum((response_vector - probability) * weights)
}


cn <- function(theta, item_parameters, index) {

  values <- Pi_IRT(theta, item_parameters)

  weights <- item_weights(
    theta = theta,
    item_parameters = item_parameters,
    index = index
  )

  numerator <- sum(values$dPi * weights)
  denominator <- sum(values$dPi * ri_IRT(theta, item_parameters))

  if (!is.finite(denominator) || abs(denominator) < 1e-12) {
    return(NA_real_)
  }

  numerator / denominator
}


corrected_weights <- function(theta, item_parameters, index) {

  correction <- cn(
    theta = theta,
    item_parameters = item_parameters,
    index = index
  )

  if (!is.finite(correction)) {
    return(rep(NA_real_, nrow(item_parameters)))
  }

  item_weights(
    theta = theta,
    item_parameters = item_parameters,
    index = index
  ) -
    correction * ri_IRT(theta, item_parameters)
}


tau2n <- function(theta, item_parameters, index) {

  probability <- Pi_IRT(theta, item_parameters)$Pi

  weights <- corrected_weights(
    theta = theta,
    item_parameters = item_parameters,
    index = index
  )

  if (anyNA(weights)) {
    return(NA_real_)
  }

  sum(
    weights^2 *
      probability *
      (1 - probability)
  ) / nrow(item_parameters)
}


compute_corrected_index <- function(
    data,
    item_parameters,
    theta,
    index = c("U", "W")) {

  index <- match.arg(index)
  corrected_index <- rep(NA_real_, nrow(data))

  for (person in seq_len(nrow(data))) {

    th <- theta[person]

    if (!is.finite(th)) {
      next
    }

    correction_c <- cn(
      theta = th,
      item_parameters = item_parameters,
      index = index
    )

    bias_WL <- r0_WL(
      theta = th,
      item_parameters = item_parameters
    )

    expected_corrected <- -correction_c * bias_WL

    variance_corrected <-
      nrow(item_parameters) *
      tau2n(
        theta = th,
        item_parameters = item_parameters,
        index = index
      )

    observed_value <- Wn(
      response_vector = data[person, ],
      theta = th,
      item_parameters = item_parameters,
      index = index
    )

    if (
      is.finite(expected_corrected) &&
      is.finite(variance_corrected) &&
      variance_corrected > 0
    ) {
      corrected_index[person] <-
        (observed_value - expected_corrected) /
        sqrt(variance_corrected)
    }
  }

  corrected_index
}


# ============================================================
# 2. SMALL HELPER FUNCTIONS
# ============================================================

empirical_skewness <- function(x) {

  x <- x[is.finite(x)]

  if (length(x) < 3L || sd(x) == 0) {
    return(NA_real_)
  }

  mean((x - mean(x))^3) / sd(x)^3
}


empirical_excess_kurtosis <- function(x) {

  x <- x[is.finite(x)]

  if (length(x) < 4L || sd(x) == 0) {
    return(NA_real_)
  }

  mean((x - mean(x))^4) / sd(x)^4 - 3
}


cohen_kappa_binary <- function(x, y) {

  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  tab <- table(
    factor(x, levels = c(FALSE, TRUE)),
    factor(y, levels = c(FALSE, TRUE))
  )

  n <- sum(tab)
  observed_agreement <- sum(diag(tab)) / n
  expected_agreement <-
    sum(rowSums(tab) * colSums(tab)) / n^2

  if (!is.finite(expected_agreement) ||
      abs(1 - expected_agreement) < 1e-12) {
    return(NA_real_)
  }

  (observed_agreement - expected_agreement) /
    (1 - expected_agreement)
}


make_ability_bins <- function(theta, n_bins = 5L) {

  out <- rep(NA_integer_, length(theta))
  keep <- is.finite(theta)
  n <- sum(keep)

  if (n == 0L) {
    return(out)
  }

  ranks <- rank(theta[keep], ties.method = "first")
  out[keep] <- pmin(
    n_bins,
    ceiling(ranks * n_bins / n)
  )

  out
}


simulate_irt_data <- function(item_parameters, theta) {

  N <- length(theta)
  J <- nrow(item_parameters)

  a <- item_parameters[, 1]
  b <- item_parameters[, 2]
  c <- item_parameters[, 3]

  theta_matrix <- matrix(theta, nrow = N, ncol = J)
  a_matrix <- matrix(a, nrow = N, ncol = J, byrow = TRUE)
  b_matrix <- matrix(b, nrow = N, ncol = J, byrow = TRUE)
  c_matrix <- matrix(c, nrow = N, ncol = J, byrow = TRUE)

  probability <-
    c_matrix +
    (1 - c_matrix) *
    plogis(a_matrix * (theta_matrix - b_matrix))

  simulated <- matrix(
    rbinom(N * J, size = 1, prob = as.vector(probability)),
    nrow = N,
    ncol = J
  )

  simulated
}


# ============================================================
# 3. FIT ONE IRT MODEL FOR MODEL COMPARISON
# ============================================================

fit_ltm_model <- function(X, model) {

  X_data <- as.data.frame(X)

  tryCatch(
    switch(
      model,
      "1PL" = ltm::rasch(
        X_data,
        IRT.param = TRUE
      ),
      "2PL" = ltm::ltm(
        X_data ~ z1,
        IRT.param = TRUE
      ),
      "3PL" = ltm::tpm(
        X_data,
        type = "latent.trait",
        IRT.param = TRUE
      ),
      stop("Unknown IRT model: ", model)
    ),
    error = function(e) e
  )
}


# ============================================================
# 4. COMPUTE PERSON-FIT STATISTICS UNDER ONE IRT MODEL
# ============================================================

compute_pfs_one_model <- function(X, model, respondent_ID, alpha = 0.05) {

  N <- nrow(X)

  fit_irtoys <- irtoys::est(
    resp = X,
    model = model,
    engine = "ltm"
  )

  item_parameters <- as.matrix(fit_irtoys$est)

  if (ncol(item_parameters) != 3L) {
    stop(model, ": item-parameter matrix does not have 3 columns.")
  }

  colnames(item_parameters) <- c("a", "b", "c")

  if (!is.null(colnames(X))) {
    rownames(item_parameters) <- colnames(X)
  }

  if (anyNA(item_parameters) || any(!is.finite(item_parameters))) {
    stop(model, ": non-finite item-parameter estimate(s).")
  }

  theta_object <- irtoys::wle(
    resp = X,
    ip = item_parameters
  )

  theta_WL <- as.numeric(theta_object[, 1])

  if (length(theta_WL) != N) {
    stop(model, ": number of theta estimates does not equal N.")
  }

  xi <- matrix(
    theta_WL,
    ncol = 1,
    dimnames = list(rownames(X), "theta")
  )

  lzstar_scores <- tryCatch(
    {
      out <- PerFit::lzstar(
        matrix = X,
        IP = item_parameters,
        IRT.PModel = model,
        Ability = theta_WL,
        Ability.PModel = "WL"
      )
      as.numeric(unlist(out$PFscores))
    },
    error = function(e) {
      warning(model, " lz*: ", conditionMessage(e))
      rep(NA_real_, N)
    }
  )

  ECI_scores <- tryCatch(
    {
      out <- aberrance::detect_pm(
        method = c("ECI2_S_TS", "ECI4_S_TS"),
        psi = item_parameters,
        xi = xi,
        x = X,
        alpha = alpha
      )

      cbind(
        ECI_Z2_star = as.numeric(out$stat[, "ECI2_S_TS"]),
        ECI_Z4_star = as.numeric(out$stat[, "ECI4_S_TS"])
      )
    },
    error = function(e) {
      warning(model, " ECI statistics: ", conditionMessage(e))
      cbind(
        ECI_Z2_star = rep(NA_real_, N),
        ECI_Z4_star = rep(NA_real_, N)
      )
    }
  )

  U_star_scores <- compute_corrected_index(
    data = X,
    item_parameters = item_parameters,
    theta = theta_WL,
    index = "U"
  )

  W_star_scores <- compute_corrected_index(
    data = X,
    item_parameters = item_parameters,
    theta = theta_WL,
    index = "W"
  )

  person_fit <- data.frame(
    Person = respondent_ID,
    Model = model,
    theta_WL = theta_WL,
    lzstar = lzstar_scores,
    U_star = U_star_scores,
    W_star = W_star_scores,
    ECI_Z2_star = ECI_scores[, "ECI_Z2_star"],
    ECI_Z4_star = ECI_scores[, "ECI_Z4_star"],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  list(
    item_parameters = item_parameters,
    theta_WL = theta_WL,
    person_fit = person_fit,
    fit_irtoys = fit_irtoys
  )
}


# ============================================================
# 5. MAIN ANALYSIS FUNCTION
# ============================================================

analyze_irw_personfit <- function(
    X,
    Person_ID = NULL,
    dataset_name = "IRW",
    alpha = 0.05,
    model_selection = c("BIC", "AIC"),
    n_ability_bins = 5L,
    make_plots = TRUE,
    save_tables = TRUE,
    output_dir = getwd(),
    run_bootstrap = FALSE,
    bootstrap_reps = 200L,
    seed = 123,
    verbose = TRUE) {

  model_selection <- match.arg(model_selection)
  models <- c("1PL", "2PL", "3PL")

  score_names <- c(
    "lzstar",
    "U_star",
    "W_star",
    "ECI_Z2_star",
    "ECI_Z4_star"
  )

  index_labels <- c(
    lzstar = "lz*",
    U_star = "U*",
    W_star = "W*",
    ECI_Z2_star = "ECI_Z2*",
    ECI_Z4_star = "ECI_Z4*"
  )

  set.seed(seed)

  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 0.5) {
    stop("alpha must be a single number strictly between 0 and 0.5.")
  }

  n_ability_bins <- as.integer(n_ability_bins)
  if (!is.finite(n_ability_bins) || n_ability_bins < 2L) {
    stop("n_ability_bins must be an integer of at least 2.")
  }

  X <- as.matrix(X)
  storage.mode(X) <- "numeric"

  if (nrow(X) < 2L || ncol(X) < 2L) {
    stop("X must contain at least 2 persons and 2 items.")
  }

  if (anyNA(X)) {
    stop(
      "This script currently requires complete dichotomous data because ",
      "irtoys::est() does not accept missing responses."
    )
  }

  if (!all(X %in% c(0, 1))) {
    stop("X must contain only dichotomous responses coded 0/1.")
  }

  constant_items <- which(
    colSums(X) == 0 | colSums(X) == nrow(X)
  )

  if (length(constant_items) > 0L) {
    stop(
      "The following item(s) have no response variation and should be removed: ",
      paste(constant_items, collapse = ", ")
    )
  }

  if (is.null(Person_ID)) {
    if (!is.null(rownames(X))) {
      respondent_ID <- rownames(X)
    } else {
      respondent_ID <- seq_len(nrow(X))
    }
  } else {
    if (length(Person_ID) != nrow(X)) {
      stop("Person_ID must have exactly nrow(X) values.")
    }
    respondent_ID <- Person_ID
  }

  safe_dataset_name <- gsub(
    "[^A-Za-z0-9_-]+",
    "_",
    dataset_name
  )

  analysis_dir <- file.path(
    output_dir,
    paste0(safe_dataset_name, "_personfit")
  )

  if (save_tables || make_plots) {
    dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  }


  # ----------------------------------------------------------
  # 5.1 Compare 1PL, 2PL, and 3PL
  # ----------------------------------------------------------

  ltm_fits <- setNames(
    lapply(models, function(model) fit_ltm_model(X, model)),
    models
  )

  model_comparison <- do.call(
    rbind,
    lapply(models, function(model) {

      fit <- ltm_fits[[model]]

      if (inherits(fit, "error")) {
        return(
          data.frame(
            Model = model,
            Fit_OK = FALSE,
            Log_likelihood = NA_real_,
            N_parameters = NA_integer_,
            AIC = NA_real_,
            BIC = NA_real_,
            Error = conditionMessage(fit),
            stringsAsFactors = FALSE
          )
        )
      }

      ll <- logLik(fit)

      data.frame(
        Model = model,
        Fit_OK = TRUE,
        Log_likelihood = as.numeric(ll),
        N_parameters = as.integer(attr(ll, "df")),
        AIC = AIC(fit),
        BIC = BIC(fit),
        Error = "",
        stringsAsFactors = FALSE
      )
    })
  )

  criterion_values <- model_comparison[[model_selection]]
  valid_model_rows <- which(is.finite(criterion_values))

  if (length(valid_model_rows) == 0L) {
    stop("None of the 1PL/2PL/3PL models could be fitted.")
  }

  min_AIC <- if (any(is.finite(model_comparison$AIC))) {
    min(model_comparison$AIC[is.finite(model_comparison$AIC)])
  } else {
    NA_real_
  }

  min_BIC <- if (any(is.finite(model_comparison$BIC))) {
    min(model_comparison$BIC[is.finite(model_comparison$BIC)])
  } else {
    NA_real_
  }

  model_comparison$Delta_AIC <- model_comparison$AIC - min_AIC
  model_comparison$Delta_BIC <- model_comparison$BIC - min_BIC

  model_comparison$Preferred_by_AIC <-
    is.finite(model_comparison$AIC) &
    model_comparison$AIC == min_AIC

  model_comparison$Preferred_by_BIC <-
    is.finite(model_comparison$BIC) &
    model_comparison$BIC == min_BIC

  working_row <- valid_model_rows[
    which.min(criterion_values[valid_model_rows])
  ]

  working_model <- model_comparison$Model[working_row]

  model_comparison$Working_model <-
    model_comparison$Model == working_model


  # ----------------------------------------------------------
  # 5.2 Compute PFSs under ALL successfully calibrated models
  # ----------------------------------------------------------

  pfs_by_model <- setNames(vector("list", length(models)), models)
  pfs_errors <- setNames(rep("", length(models)), models)

  for (model in models) {

    pfs_by_model[[model]] <- tryCatch(
      compute_pfs_one_model(
        X = X,
        model = model,
        respondent_ID = respondent_ID,
        alpha = alpha
      ),
      error = function(e) {
        pfs_errors[model] <<- conditionMessage(e)
        NULL
      }
    )
  }

  model_comparison$PFS_OK <- vapply(
    models,
    function(model) !is.null(pfs_by_model[[model]]),
    logical(1)
  )

  model_comparison$PFS_Error <- unname(pfs_errors[models])

  if (is.null(pfs_by_model[[working_model]])) {
    stop(
      "The working model (", working_model,
      ") was selected by ", model_selection,
      " but its PFS computation failed: ",
      pfs_errors[working_model]
    )
  }

  successful_models <- models[
    vapply(models, function(m) !is.null(pfs_by_model[[m]]), logical(1))
  ]

  person_fit_results <- do.call(
    rbind,
    lapply(successful_models, function(model) {
      pfs_by_model[[model]]$person_fit
    })
  )

  rownames(person_fit_results) <- NULL


  # ----------------------------------------------------------
  # 5.3 Flagging at alpha under each model
  # ----------------------------------------------------------

  normal_cutoff <- qnorm(alpha)
  upper_cutoff <- qnorm(1 - alpha)

  for (index_name in score_names) {
    person_fit_results[[paste0(index_name, "_flag")]] <-
      if (index_name == "U_star") {
        person_fit_results[[index_name]] > upper_cutoff
      } else {
        person_fit_results[[index_name]] < normal_cutoff
      }
  }

  # Put flagged values back into each model-specific object.
  for (model in successful_models) {
    idx <- person_fit_results$Model == model
    pfs_by_model[[model]]$person_fit <- person_fit_results[idx, ]
  }

  flag_summary <- do.call(
    rbind,
    lapply(successful_models, function(model) {

      model_data <- pfs_by_model[[model]]$person_fit

      do.call(
        rbind,
        lapply(score_names, function(index_name) {

          values <- model_data[[index_name]]
          flags <- model_data[[paste0(index_name, "_flag")]]

          n_valid <- sum(is.finite(values))
          n_flagged <- sum(flags, na.rm = TRUE)

          data.frame(
            Model = model,
            Index = unname(index_labels[index_name]),
            N_valid = n_valid,
            N_flagged = n_flagged,
            Percent_flagged = if (n_valid > 0L) {
              100 * n_flagged / n_valid
            } else {
              NA_real_
            },
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )


  # ----------------------------------------------------------
  # 5.4 Empirical distribution summaries
  #     Concise set: mean, SD, skewness, kurtosis
  # ----------------------------------------------------------

  descriptive_statistics <- do.call(
    rbind,
    lapply(successful_models, function(model) {

      model_data <- pfs_by_model[[model]]$person_fit

      do.call(
        rbind,
        lapply(score_names, function(index_name) {

          values <- model_data[[index_name]]
          values <- values[is.finite(values)]

          data.frame(
            Model = model,
            Index = unname(index_labels[index_name]),
            N_valid = length(values),
            Mean = if (length(values) > 0L) mean(values) else NA_real_,
            SD = if (length(values) > 1L) sd(values) else NA_real_,
            Skewness = empirical_skewness(values),
            Excess_kurtosis = empirical_excess_kurtosis(values),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )


  # ----------------------------------------------------------
  # 5.5 PFS behavior conditional on ability
  # ----------------------------------------------------------

  ability_conditional_summary <- do.call(
    rbind,
    lapply(successful_models, function(model) {

      model_data <- pfs_by_model[[model]]$person_fit
      model_data$Ability_bin <- make_ability_bins(
        model_data$theta_WL,
        n_bins = n_ability_bins
      )

      do.call(
        rbind,
        lapply(score_names, function(index_name) {

          do.call(
            rbind,
            lapply(seq_len(n_ability_bins), function(bin_number) {

              use <- model_data$Ability_bin == bin_number
              values <- model_data[[index_name]][use]
              theta_values <- model_data$theta_WL[use]
              flags <- model_data[[paste0(index_name, "_flag")]][use]

              valid_score <- is.finite(values)
              n_valid <- sum(valid_score)

              data.frame(
                Model = model,
                Index = unname(index_labels[index_name]),
                Ability_bin = bin_number,
                N = sum(use, na.rm = TRUE),
                Mean_theta = if (any(is.finite(theta_values))) {
                  mean(theta_values[is.finite(theta_values)])
                } else {
                  NA_real_
                },
                Mean_PFS = if (n_valid > 0L) {
                  mean(values[valid_score])
                } else {
                  NA_real_
                },
                SD_PFS = if (n_valid > 1L) {
                  sd(values[valid_score])
                } else {
                  NA_real_
                },
                N_flagged = sum(flags, na.rm = TRUE),
                Percent_flagged = if (n_valid > 0L) {
                  100 * sum(flags, na.rm = TRUE) / n_valid
                } else {
                  NA_real_
                },
                stringsAsFactors = FALSE
              )
            })
          )
        })
      )
    })
  )


  # ----------------------------------------------------------
  # 5.6 Stability across IRT analysis models
  #     Score stability: Pearson and Spearman correlations
  #     Flag stability: agreement rate and Cohen's kappa
  # ----------------------------------------------------------

  model_pairs <- if (length(successful_models) >= 2L) {
    combn(successful_models, 2, simplify = FALSE)
  } else {
    list()
  }

  if (length(model_pairs) > 0L) {

    model_stability <- do.call(
      rbind,
      lapply(model_pairs, function(pair) {

        model_1 <- pair[1]
        model_2 <- pair[2]

        data_1 <- pfs_by_model[[model_1]]$person_fit
        data_2 <- pfs_by_model[[model_2]]$person_fit

        do.call(
          rbind,
          lapply(score_names, function(index_name) {

            x <- data_1[[index_name]]
            y <- data_2[[index_name]]

            keep_scores <- is.finite(x) & is.finite(y)

            flag_1 <- data_1[[paste0(index_name, "_flag")]]
            flag_2 <- data_2[[paste0(index_name, "_flag")]]
            keep_flags <- !is.na(flag_1) & !is.na(flag_2)

            agreement <- if (any(keep_flags)) {
              mean(flag_1[keep_flags] == flag_2[keep_flags])
            } else {
              NA_real_
            }

            data.frame(
              Model_1 = model_1,
              Model_2 = model_2,
              Index = unname(index_labels[index_name]),
              N_common = sum(keep_scores),
              Pearson_r = if (sum(keep_scores) > 2L) {
                cor(x[keep_scores], y[keep_scores], method = "pearson")
              } else {
                NA_real_
              },
              Spearman_r = if (sum(keep_scores) > 2L) {
                cor(x[keep_scores], y[keep_scores], method = "spearman")
              } else {
                NA_real_
              },
              Agreement_rate = agreement,
              Percent_changed = if (is.finite(agreement)) {
                100 * (1 - agreement)
              } else {
                NA_real_
              },
              Cohen_kappa = cohen_kappa_binary(flag_1, flag_2),
              stringsAsFactors = FALSE
            )
          })
        )
      })
    )

  } else {

    model_stability <- data.frame(
      Model_1 = character(0),
      Model_2 = character(0),
      Index = character(0),
      N_common = integer(0),
      Pearson_r = numeric(0),
      Spearman_r = numeric(0),
      Agreement_rate = numeric(0),
      Percent_changed = numeric(0),
      Cohen_kappa = numeric(0),
      stringsAsFactors = FALSE
    )
  }


  # ----------------------------------------------------------
  # 5.7 Combine item parameters across models
  # ----------------------------------------------------------

  item_parameters_all_models <- do.call(
    rbind,
    lapply(successful_models, function(model) {

      ip <- pfs_by_model[[model]]$item_parameters

      data.frame(
        Model = model,
        Item = if (!is.null(rownames(ip))) {
          rownames(ip)
        } else {
          seq_len(nrow(ip))
        },
        a = ip[, 1],
        b = ip[, 2],
        c = ip[, 3],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    })
  )


  # ----------------------------------------------------------
  # 5.8 Distribution plots: histogram/density + Q-Q plot
  # ----------------------------------------------------------

  distribution_plot_file <- NULL

  if (make_plots) {

    distribution_plot_file <- file.path(
      analysis_dir,
      "person_fit_distributions.pdf"
    )

    grDevices::pdf(
      file = distribution_plot_file,
      width = 10,
      height = 5
    )

    old_par <- graphics::par(no.readonly = TRUE)

    for (model in successful_models) {

      model_data <- pfs_by_model[[model]]$person_fit

      for (index_name in score_names) {

        values <- model_data[[index_name]]
        values <- values[is.finite(values)]

        graphics::par(mfrow = c(1, 2))

        if (length(values) > 1L) {

          graphics::hist(
            values,
            breaks = "FD",
            probability = TRUE,
            main = paste(model, "-", unname(index_labels[index_name])),
            xlab = unname(index_labels[index_name]),
            ylab = "Density"
          )

          if (length(unique(values)) > 1L) {
            graphics::lines(
              stats::density(values),
              lwd = 2
            )
          }

          graphics::curve(
            stats::dnorm(x),
            add = TRUE,
            lty = 2,
            lwd = 2
          )

          graphics::legend(
            "topright",
            legend = c("Empirical density", "N(0,1)"),
            lty = c(1, 2),
            lwd = c(2, 2),
            bty = "n"
          )

          stats::qqnorm(
            values,
            main = paste("Q-Q:", model, "-", unname(index_labels[index_name]))
          )
          stats::qqline(values, lwd = 2)

        } else {

          graphics::plot.new()
          graphics::title(main = paste(model, "-", unname(index_labels[index_name])))
          graphics::text(0.5, 0.5, "Insufficient valid values")

          graphics::plot.new()
          graphics::title(main = "Q-Q plot")
          graphics::text(0.5, 0.5, "Insufficient valid values")
        }
      }
    }

    graphics::par(old_par)
    grDevices::dev.off()
  }


  # ----------------------------------------------------------
  # 5.9 Ability-conditional plots
  # ----------------------------------------------------------

  ability_plot_file <- NULL

  if (make_plots) {

    ability_plot_file <- file.path(
      analysis_dir,
      "person_fit_by_ability.pdf"
    )

    grDevices::pdf(
      file = ability_plot_file,
      width = 10,
      height = 5
    )

    for (model in successful_models) {

      model_data <- pfs_by_model[[model]]$person_fit

      for (index_name in score_names) {

        values <- model_data[[index_name]]
        theta_values <- model_data$theta_WL
        keep <- is.finite(values) & is.finite(theta_values)

        graphics::par(mfrow = c(1, 2))

        if (sum(keep) > 2L) {

          graphics::plot(
            theta_values[keep],
            values[keep],
            pch = 16,
            cex = 0.45,
            xlab = expression(theta[WL]),
            ylab = unname(index_labels[index_name]),
            main = paste(model, "- PFS by ability")
          )

          lowess_fit <- stats::lowess(
            theta_values[keep],
            values[keep],
            f = 2 / 3
          )

          graphics::lines(
            lowess_fit$x,
            lowess_fit$y,
            lwd = 2
          )

          graphics::abline(h = 0, lty = 2)

        } else {

          graphics::plot.new()
          graphics::title(main = paste(model, "- PFS by ability"))
          graphics::text(0.5, 0.5, "Insufficient valid values")
        }

        binned <- ability_conditional_summary[
          ability_conditional_summary$Model == model &
          ability_conditional_summary$Index == unname(index_labels[index_name]),
        ]

        keep_bins <-
          is.finite(binned$Mean_theta) &
          is.finite(binned$Percent_flagged)

        if (sum(keep_bins) > 0L) {

          y_values <- binned$Percent_flagged[keep_bins] / 100
          y_max <- max(c(y_values, alpha, 0.10), na.rm = TRUE)

          graphics::plot(
            binned$Mean_theta[keep_bins],
            y_values,
            type = "b",
            pch = 16,
            ylim = c(0, min(1, y_max * 1.15)),
            xlab = expression(theta[WL]),
            ylab = "Empirical flagging rate",
            main = paste(model, "- Flagging by ability")
          )

          graphics::abline(h = alpha, lty = 2, lwd = 2)

        } else {

          graphics::plot.new()
          graphics::title(main = paste(model, "- Flagging by ability"))
          graphics::text(0.5, 0.5, "Insufficient valid values")
        }
      }
    }

    grDevices::dev.off()
  }


  # ----------------------------------------------------------
  # 5.10 Optional IRW-informed parametric bootstrap
  # ----------------------------------------------------------
  # The observed WL theta distribution and item parameters from
  # the working model are treated as the generating values.
  # Each replicated dataset is then recalibrated under that same
  # working model before the PFS flagging rates are computed.

  bootstrap_details <- NULL
  bootstrap_typeI <- NULL

  if (run_bootstrap) {

    bootstrap_reps <- as.integer(bootstrap_reps)

    if (!is.finite(bootstrap_reps) || bootstrap_reps < 1L) {
      stop("bootstrap_reps must be an integer of at least 1.")
    }

    generating_ip <- pfs_by_model[[working_model]]$item_parameters
    generating_theta <- pfs_by_model[[working_model]]$theta_WL

    if (any(!is.finite(generating_theta))) {
      stop(
        "The parametric bootstrap requires finite WL theta estimates for all persons."
      )
    }

    bootstrap_rows <- vector("list", bootstrap_reps)

    for (b in seq_len(bootstrap_reps)) {

      X_boot <- simulate_irt_data(
        item_parameters = generating_ip,
        theta = generating_theta
      )

      dimnames(X_boot) <- dimnames(X)

      boot_fit <- tryCatch(
        compute_pfs_one_model(
          X = X_boot,
          model = working_model,
          respondent_ID = respondent_ID,
          alpha = alpha
        ),
        error = function(e) e
      )

      if (inherits(boot_fit, "error")) {

        bootstrap_rows[[b]] <- data.frame(
          Replicate = b,
          Index = unname(index_labels[score_names]),
          N_valid = NA_integer_,
          N_flagged = NA_integer_,
          Flag_rate = NA_real_,
          Error = conditionMessage(boot_fit),
          stringsAsFactors = FALSE
        )

      } else {

        boot_data <- boot_fit$person_fit

        bootstrap_rows[[b]] <- do.call(
          rbind,
          lapply(score_names, function(index_name) {

            values <- boot_data[[index_name]]
            flags <- if (index_name == "U_star") {
              values > upper_cutoff
            } else {
              values < normal_cutoff
            }
            n_valid <- sum(is.finite(values))
            n_flagged <- sum(flags, na.rm = TRUE)

            data.frame(
              Replicate = b,
              Index = unname(index_labels[index_name]),
              N_valid = n_valid,
              N_flagged = n_flagged,
              Flag_rate = if (n_valid > 0L) {
                n_flagged / n_valid
              } else {
                NA_real_
              },
              Error = "",
              stringsAsFactors = FALSE
            )
          })
        )
      }

      if (verbose && (b %% max(1L, floor(bootstrap_reps / 10L)) == 0L)) {
        message(
          "Parametric bootstrap: ",
          b,
          "/",
          bootstrap_reps,
          " replicates completed."
        )
      }
    }

    bootstrap_details <- do.call(rbind, bootstrap_rows)

    bootstrap_typeI <- do.call(
      rbind,
      lapply(unname(index_labels[score_names]), function(index_label) {

        rates <- bootstrap_details$Flag_rate[
          bootstrap_details$Index == index_label
        ]
        rates <- rates[is.finite(rates)]

        data.frame(
          Working_model = working_model,
          Index = index_label,
          Nominal_alpha = alpha,
          N_successful_replicates = length(rates),
          Mean_TypeI = if (length(rates) > 0L) mean(rates) else NA_real_,
          SD_TypeI = if (length(rates) > 1L) sd(rates) else NA_real_,
          Monte_Carlo_SE = if (length(rates) > 1L) {
            sd(rates) / sqrt(length(rates))
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      })
    )
  }


  # ----------------------------------------------------------
  # 5.11 Save analysis tables
  # ----------------------------------------------------------

  if (save_tables) {

    utils::write.csv(
      model_comparison,
      file.path(analysis_dir, "model_comparison.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      item_parameters_all_models,
      file.path(analysis_dir, "item_parameters_all_models.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      person_fit_results,
      file.path(analysis_dir, "person_fit_scores_all_models.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      flag_summary,
      file.path(analysis_dir, "flagging_summary.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      descriptive_statistics,
      file.path(analysis_dir, "distribution_summary.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      ability_conditional_summary,
      file.path(analysis_dir, "ability_conditional_summary.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      model_stability,
      file.path(analysis_dir, "model_stability.csv"),
      row.names = FALSE
    )

    if (run_bootstrap) {

      utils::write.csv(
        bootstrap_details,
        file.path(analysis_dir, "bootstrap_details.csv"),
        row.names = FALSE
      )

      utils::write.csv(
        bootstrap_typeI,
        file.path(analysis_dir, "bootstrap_typeI_summary.csv"),
        row.names = FALSE
      )
    }
  }


  # ----------------------------------------------------------
  # 5.12 Console summary
  # ----------------------------------------------------------

  if (verbose) {

    cat("\n============================================================\n")
    cat("IRW PERSON-FIT ANALYSIS:", dataset_name, "\n")
    cat("============================================================\n\n")

    print(model_comparison)

    cat(
      "\nWorking model according to ",
      model_selection,
      ": ",
      working_model,
      "\n",
      sep = ""
    )

    cat("\nEmpirical flagging rates:\n")
    print(flag_summary)

    cat("\nEmpirical distribution summaries:\n")
    print(descriptive_statistics)

    if (nrow(model_stability) > 0L) {
      cat("\nStability across IRT models:\n")
      print(model_stability)
    }

    if (run_bootstrap) {
      cat("\nIRW-informed parametric bootstrap Type I error:\n")
      print(bootstrap_typeI)
    }

    if (save_tables || make_plots) {
      cat("\nOutput directory:\n", analysis_dir, "\n", sep = "")
    }
  }


  invisible(
    list(
      dataset_name = dataset_name,
      model_selection_criterion = model_selection,
      working_model = working_model,
      model_comparison = model_comparison,
      item_parameters = item_parameters_all_models,
      person_fit = person_fit_results,
      flagging = flag_summary,
      distributions = descriptive_statistics,
      ability_conditional = ability_conditional_summary,
      model_stability = model_stability,
      bootstrap_details = bootstrap_details,
      bootstrap_typeI = bootstrap_typeI,
      model_objects = ltm_fits,
      pfs_objects = pfs_by_model,
      distribution_plot_file = distribution_plot_file,
      ability_plot_file = ability_plot_file,
      output_dir = analysis_dir
    )
  )
}


