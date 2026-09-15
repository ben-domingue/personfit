#############################################################
# IRW PERSON-FIT ANALYSES
#
# Revision of 2026-09-14 (from IRW_personfit_Sept09092026.R),
# responding to collaborator feedback:
#   (1) Absolute fit. AIC/BIC are complemented by M2, RMSEA and
#       SRMSR for every calibrated model, so the working model can
#       be checked for acceptable fit, not only relative fit
#       (section 5.2b).
#   (2) The parametric bootstrap now keeps the replicate PFS
#       scores (not just flag rates) and can rescore each
#       replicate under every analysis model, not only the
#       working model (section 5.10).
#   (3) Empirical critical values: the alpha (or 1 - alpha for U*)
#       quantile of the pooled bootstrap scores replaces the
#       normal-theory cutoff, and observed respondents are
#       reflagged with it (section 5.11).
#   (4) Ability-conditional calibration: the same critical values
#       derived within ability quintiles, with respondents
#       reflagged using the cutoff for their own quintile
#       (section 5.11). Null rates for both calibrated methods are
#       estimated out of sample (odd/even replicate split).
# Everything else is unchanged from the Sept 9 engine.
#############################################################

require(ltm)
require(irtoys)
require(PerFit)
require(aberrance)
require(mirt)
require(parallel)

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


# Flag a PFS vector against a cutoff (scalar, or one per person).
# U* flags in the upper tail; the other four indices in the lower.
flag_scores <- function(values, index_name, cutoff) {

  if (index_name == "U_star") {
    values > cutoff
  } else {
    values < cutoff
  }
}


# Quantile of the null distribution that serves as the cutoff.
tail_probability <- function(index_name, alpha) {

  if (index_name == "U_star") 1 - alpha else alpha
}


# Proportion of finite scores flagged against `cutoff`.
flag_rate <- function(values, index_name, cutoff) {

  cutoff <- rep_len(cutoff, length(values))
  keep <- is.finite(values) & is.finite(cutoff)

  if (!any(keep)) {
    return(NA_real_)
  }

  mean(flag_scores(values[keep], index_name, cutoff[keep]))
}


# Ability-bin cut points are taken from the OBSERVED theta
# distribution and then applied unchanged to bootstrap replicates,
# so that "quintile 1" means the same theta range in both.
# Cut points rather than ranks are used so that tied theta values
# (e.g. equal sum scores under the 1PL) always share a bin; bins
# can therefore be somewhat unequal in size.
ability_cut_points <- function(theta, n_bins) {

  stats::quantile(
    theta[is.finite(theta)],
    probs = seq_len(n_bins - 1L) / n_bins,
    names = FALSE
  )
}


assign_ability_bins <- function(theta, cut_points) {

  out <- findInterval(theta, cut_points) + 1L
  out[!is.finite(theta)] <- NA_integer_
  out
}


# ============================================================
# 3b. ABSOLUTE FIT (M2, RMSEA, SRMSR) FOR ONE CALIBRATED MODEL
# ============================================================
# AIC/BIC only rank the models against each other. The limited-
# information M2 statistic (Maydeu-Olivares & Joe, 2006) and the
# RMSEA and SRMSR derived from it describe absolute fit.
#
# The irtoys item parameters actually used for person-fit scoring
# are loaded into mirt as starting values and not re-estimated
# (TOL = NaN), so the indices describe exactly the calibration the
# PFSs are computed from. Degrees of freedom still count the
# model's free parameters (J for the 1PL's d's plus one common a;
# 2J for the 2PL; 3J for the 3PL).

absolute_fit_one_model <- function(X, item_parameters, model) {

  J <- ncol(X)
  X_data <- as.data.frame(X)
  item_names <- paste0("Item", seq_len(J))
  colnames(X_data) <- item_names

  itemtype <- if (model == "3PL") "3PL" else "2PL"

  model_syntax <- if (model == "1PL") {
    mirt::mirt.model(
      sprintf("F = 1-%d\nCONSTRAIN = (1-%d, a1)", J, J)
    )
  } else {
    1
  }

  values <- mirt::mirt(
    X_data,
    model_syntax,
    itemtype = itemtype,
    pars = "values",
    verbose = FALSE
  )

  a <- item_parameters[, 1]
  b <- item_parameters[, 2]
  c <- item_parameters[, 3]

  for (j in seq_len(J)) {

    rows <- values$item == item_names[j]

    values$value[rows & values$name == "a1"] <- a[j]
    values$value[rows & values$name == "d"] <- -a[j] * b[j]

    if (model == "3PL") {
      values$value[rows & values$name == "g"] <- c[j]
    }
  }

  fit <- mirt::mirt(
    X_data,
    model_syntax,
    itemtype = itemtype,
    pars = values,
    TOL = NaN,
    verbose = FALSE
  )

  m2 <- mirt::M2(fit)

  data.frame(
    M2 = m2$M2,
    M2_df = m2$df,
    M2_p = m2$p,
    RMSEA = m2$RMSEA,
    RMSEA_5 = m2$RMSEA_5,
    RMSEA_95 = m2$RMSEA_95,
    SRMSR = m2$SRMSR,
    CFI = m2$CFI,
    TLI = m2$TLI,
    stringsAsFactors = FALSE
  )
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
    verbose = TRUE,
    compute_absolute_fit = TRUE,
    rmsea_cutoff = 0.05,
    srmsr_cutoff = 0.05,
    bootstrap_analysis_models = c("all", "working"),
    bootstrap_cores = 1L,
    min_boot_scores_per_bin = 200L,
    max_modal_share_per_bin = 0.5,
    return_bootstrap_scores = FALSE) {

  model_selection <- match.arg(model_selection)
  bootstrap_analysis_models <- match.arg(bootstrap_analysis_models)
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
  # 5.2b Absolute fit of each calibrated model: M2, RMSEA, SRMSR
  # ----------------------------------------------------------
  # Acceptable_fit uses RMSEA <= rmsea_cutoff and SRMSR <=
  # srmsr_cutoff (both 0.05 by default, following Maydeu-Olivares
  # & Joe, 2014). The working model is still chosen by AIC/BIC;
  # this only reports whether it fits acceptably.

  absolute_fit_columns <- c(
    "M2", "M2_df", "M2_p", "RMSEA", "RMSEA_5", "RMSEA_95",
    "SRMSR", "CFI", "TLI"
  )

  for (column in absolute_fit_columns) {
    model_comparison[[column]] <- NA_real_
  }
  model_comparison$Acceptable_fit <- NA
  model_comparison$Absolute_fit_Error <- ""

  if (compute_absolute_fit) {

    for (model in successful_models) {

      row <- which(model_comparison$Model == model)

      absolute_fit <- tryCatch(
        absolute_fit_one_model(
          X = X,
          item_parameters = pfs_by_model[[model]]$item_parameters,
          model = model
        ),
        error = function(e) e
      )

      if (inherits(absolute_fit, "error")) {
        model_comparison$Absolute_fit_Error[row] <-
          conditionMessage(absolute_fit)
        next
      }

      for (column in absolute_fit_columns) {
        model_comparison[[column]][row] <- absolute_fit[[column]]
      }

      model_comparison$Acceptable_fit[row] <-
        is.finite(absolute_fit$RMSEA) &&
        is.finite(absolute_fit$SRMSR) &&
        absolute_fit$RMSEA <= rmsea_cutoff &&
        absolute_fit$SRMSR <= srmsr_cutoff
    }
  }

  working_model_acceptable <-
    model_comparison$Acceptable_fit[model_comparison$Model == working_model]


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
  # Each replicated dataset is then recalibrated and rescored under
  # every analysis model in `bootstrap_models`: only the working
  # model if bootstrap_analysis_models = "working", or all
  # successfully calibrated models if "all". For a non-working
  # analysis model this gives the null distribution of the PFSs
  # when the data follow the working model but are scored under
  # another model.
  #
  # Replicates run in parallel via parallel::mclapply when
  # bootstrap_cores > 1 (Unix only). Parallel runs use the
  # L'Ecuyer-CMRG generator, so results are reproducible for a
  # given seed AND bootstrap_cores, but differ from a serial run.

  bootstrap_details <- NULL
  bootstrap_typeI <- NULL
  bootstrap_scores <- NULL
  critical_values <- NULL
  calibrated_flagging <- NULL
  calibrated_ability_conditional <- NULL
  bootstrap_models <- character(0)

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

    bootstrap_models <- if (bootstrap_analysis_models == "all") {
      successful_models
    } else {
      working_model
    }

    score_columns <- c("theta_WL", score_names)

    run_one_replicate <- function(b) {

      X_boot <- simulate_irt_data(
        item_parameters = generating_ip,
        theta = generating_theta
      )

      dimnames(X_boot) <- dimnames(X)

      lapply(
        setNames(bootstrap_models, bootstrap_models),
        function(model) {

          boot_fit <- tryCatch(
            compute_pfs_one_model(
              X = X_boot,
              model = model,
              respondent_ID = respondent_ID,
              alpha = alpha
            ),
            error = function(e) e
          )

          if (inherits(boot_fit, "error")) {
            return(list(error = conditionMessage(boot_fit), scores = NULL))
          }

          list(
            error = "",
            scores = as.matrix(boot_fit$person_fit[, score_columns])
          )
        }
      )
    }

    bootstrap_cores <- max(1L, as.integer(bootstrap_cores))
    use_parallel <- bootstrap_cores > 1L && .Platform$OS.type == "unix"

    if (verbose) {
      message(
        "Parametric bootstrap: ", bootstrap_reps, " replicates x ",
        length(bootstrap_models), " analysis model(s), ",
        if (use_parallel) paste(bootstrap_cores, "cores") else "serial",
        "."
      )
    }

    if (use_parallel) {

      old_rng_kind <- RNGkind()[1]
      RNGkind("L'Ecuyer-CMRG")
      set.seed(seed)

      replicate_output <- parallel::mclapply(
        seq_len(bootstrap_reps),
        run_one_replicate,
        mc.cores = bootstrap_cores,
        mc.set.seed = TRUE
      )

      RNGkind(old_rng_kind)

    } else {

      replicate_output <- lapply(seq_len(bootstrap_reps), run_one_replicate)
    }

    # A crashed forked worker returns a try-error for its replicate;
    # record it as a failure under every analysis model.
    replicate_output <- lapply(replicate_output, function(rep_out) {
      if (inherits(rep_out, "try-error")) {
        setNames(
          lapply(bootstrap_models, function(m) {
            list(error = as.character(rep_out), scores = NULL)
          }),
          bootstrap_models
        )
      } else {
        rep_out
      }
    })


    # -- Normal-theory Type I error per replicate ---------------

    bootstrap_details <- do.call(
      rbind,
      lapply(seq_len(bootstrap_reps), function(b) {
        do.call(
          rbind,
          lapply(bootstrap_models, function(model) {

            out <- replicate_output[[b]][[model]]

            if (is.null(out$scores)) {
              return(data.frame(
                Replicate = b,
                Analysis_model = model,
                Index = unname(index_labels[score_names]),
                N_valid = NA_integer_,
                N_flagged = NA_integer_,
                Flag_rate = NA_real_,
                Error = out$error,
                stringsAsFactors = FALSE
              ))
            }

            do.call(
              rbind,
              lapply(score_names, function(index_name) {

                values <- out$scores[, index_name]
                cutoff <- if (index_name == "U_star") upper_cutoff else normal_cutoff
                n_valid <- sum(is.finite(values))
                n_flagged <- sum(
                  flag_scores(values, index_name, cutoff),
                  na.rm = TRUE
                )

                data.frame(
                  Replicate = b,
                  Analysis_model = model,
                  Index = unname(index_labels[index_name]),
                  N_valid = n_valid,
                  N_flagged = n_flagged,
                  Flag_rate = if (n_valid > 0L) n_flagged / n_valid else NA_real_,
                  Error = "",
                  stringsAsFactors = FALSE
                )
              })
            )
          })
        )
      })
    )

    bootstrap_typeI <- do.call(
      rbind,
      lapply(bootstrap_models, function(model) {
        do.call(
          rbind,
          lapply(unname(index_labels[score_names]), function(index_label) {

            rates <- bootstrap_details$Flag_rate[
              bootstrap_details$Analysis_model == model &
                bootstrap_details$Index == index_label
            ]
            rates <- rates[is.finite(rates)]

            data.frame(
              Working_model = working_model,
              Analysis_model = model,
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
      })
    )


    # -- Pool replicate scores per analysis model --------------
    # One matrix per model: Replicate, theta_WL, five PFSs, with
    # N rows per successful replicate.

    bootstrap_scores <- lapply(
      setNames(bootstrap_models, bootstrap_models),
      function(model) {
        pieces <- lapply(seq_len(bootstrap_reps), function(b) {
          scores <- replicate_output[[b]][[model]]$scores
          if (is.null(scores)) return(NULL)
          cbind(Replicate = b, scores)
        })
        pieces <- Filter(Negate(is.null), pieces)
        if (length(pieces) == 0L) return(NULL)
        do.call(rbind, pieces)
      }
    )

    rm(replicate_output)
  }


  # ----------------------------------------------------------
  # 5.11 Bootstrap-calibrated critical values and reflagging
  # ----------------------------------------------------------
  # Three flagging methods, per analysis model and index:
  #   "Normal theory"        cutoff qnorm(alpha) (qnorm(1 - alpha)
  #                          for U*), as in 5.3.
  #   "Bootstrap"            cutoff = alpha (1 - alpha) quantile of
  #                          all pooled bootstrap scores.
  #   "Bootstrap by ability" the same quantile, computed separately
  #                          within each ability bin; each observed
  #                          respondent is flagged against the
  #                          cutoff for their own bin.
  # Ability bins use cut points from the observed theta_WL under
  # that analysis model; bootstrap persons are assigned by their
  # re-estimated theta_WL. A bin with fewer than
  # min_boot_scores_per_bin finite bootstrap scores falls back to
  # the unconditional bootstrap cutoff. So does a degenerate bin: one
  # where, on average over replicates, more than
  # max_modal_share_per_bin of the scores share a single value (e.g.
  # a top quintile made up entirely of perfect scores on a short
  # test). A quantile of a point mass would flag either everyone or
  # no one in that bin.
  #
  # Null rates. For normal theory the null rate is the bootstrap
  # flagging rate. For the two bootstrap methods, evaluating the
  # cutoff on the same scores it was derived from would return
  # alpha by construction, so the null rate is estimated out of
  # sample: cutoffs from odd-numbered replicates are applied to
  # even-numbered ones and vice versa, and the two rates averaged.
  # With discrete scores (short tests) this can differ from alpha.

  if (run_bootstrap) {

    for (method_suffix in c("_flag_boot", "_flag_bootcond")) {
      for (index_name in score_names) {
        person_fit_results[[paste0(index_name, method_suffix)]] <- NA
      }
    }
    person_fit_results$Ability_bin_cal <- NA_integer_

    method_labels <- c(
      normal = "Normal theory",
      boot = "Bootstrap",
      bootcond = "Bootstrap by ability"
    )

    critical_value_rows <- list()
    flagging_rows <- list()
    ability_rows <- list()

    for (model in bootstrap_models) {

      S <- bootstrap_scores[[model]]

      if (is.null(S)) {
        next
      }

      observed_idx <- which(person_fit_results$Model == model)
      observed <- person_fit_results[observed_idx, ]

      cut_points <- ability_cut_points(observed$theta_WL, n_ability_bins)
      observed_bin <- assign_ability_bins(observed$theta_WL, cut_points)
      boot_bin <- assign_ability_bins(S[, "theta_WL"], cut_points)
      bin_edges <- c(-Inf, cut_points, Inf)

      person_fit_results$Ability_bin_cal[observed_idx] <- observed_bin

      replicate_id <- S[, "Replicate"]
      odd_half <- replicate_id %% 2 == 1

      for (index_name in score_names) {

        index_label <- unname(index_labels[index_name])
        p <- tail_probability(index_name, alpha)
        tail <- if (index_name == "U_star") "upper" else "lower"
        nominal_cutoff <- qnorm(p)

        v <- S[, index_name]
        finite_v <- is.finite(v)

        quantile_or_na <- function(x, min_n = 1L) {
          x <- x[is.finite(x)]
          if (length(x) < min_n) return(NA_real_)
          stats::quantile(x, probs = p, names = FALSE)
        }

        # Monte Carlo SE of a pooled quantile: SD of the per-replicate
        # quantiles divided by sqrt(number of replicates).
        mc_se <- function(use) {
          per_rep <- tapply(v[use], replicate_id[use], quantile_or_na)
          per_rep <- per_rep[is.finite(per_rep)]
          if (length(per_rep) < 2L) return(NA_real_)
          sd(per_rep) / sqrt(length(per_rep))
        }

        # Unconditional cutoff: full pool and each half.
        cv_all <- quantile_or_na(v)
        cv_odd <- quantile_or_na(v[odd_half])
        cv_even <- quantile_or_na(v[!odd_half])

        # Conditional cutoffs: full pool and each half, per bin,
        # falling back to the matching unconditional cutoff.
        # Mean within-replicate share of the most common score value.
        modal_share <- function(use) {
          if (!any(use)) return(NA_real_)
          mean(tapply(v[use], replicate_id[use], function(z) {
            max(table(round(z, 6))) / length(z)
          }))
        }

        degenerate_bin <- vapply(seq_len(n_ability_bins), function(bin) {
          share <- modal_share(boot_bin == bin & finite_v & !is.na(boot_bin))
          is.finite(share) && share > max_modal_share_per_bin
        }, logical(1))

        bin_cutoffs <- function(use_rows, fallback) {
          vapply(seq_len(n_ability_bins), function(bin) {
            if (degenerate_bin[bin]) return(fallback)
            use <- use_rows & boot_bin == bin & finite_v
            cv <- quantile_or_na(v[use], min_n = min_boot_scores_per_bin)
            if (is.finite(cv)) cv else fallback
          }, numeric(1))
        }

        cv_bin_all <- bin_cutoffs(rep(TRUE, length(v)), cv_all)
        cv_bin_odd <- bin_cutoffs(odd_half, cv_odd)
        cv_bin_even <- bin_cutoffs(!odd_half, cv_even)

        n_boot_bin <- vapply(seq_len(n_ability_bins), function(bin) {
          sum(boot_bin == bin & finite_v, na.rm = TRUE)
        }, integer(1))

        # Per-score cutoffs for the out-of-sample null rates: each
        # score is judged against the cutoff from the other half.
        holdout_uncond <- ifelse(odd_half, cv_even, cv_odd)
        holdout_cond <- ifelse(
          odd_half,
          cv_bin_even[boot_bin],
          cv_bin_odd[boot_bin]
        )

        null_cutoffs <- list(
          normal = rep(nominal_cutoff, length(v)),
          boot = holdout_uncond,
          bootcond = holdout_cond
        )

        # Observed reflagging.
        observed_values <- observed[[index_name]]
        observed_cutoffs <- list(
          normal = rep(nominal_cutoff, length(observed_values)),
          boot = rep(cv_all, length(observed_values)),
          bootcond = cv_bin_all[observed_bin]
        )
        observed_flags <- lapply(observed_cutoffs, function(cutoff) {
          flag_scores(observed_values, index_name, cutoff)
        })

        person_fit_results[[paste0(index_name, "_flag_boot")]][observed_idx] <-
          observed_flags$boot
        person_fit_results[[paste0(index_name, "_flag_bootcond")]][observed_idx] <-
          observed_flags$bootcond

        # Critical-value table.
        critical_value_rows[[length(critical_value_rows) + 1L]] <- data.frame(
          Working_model = working_model,
          Analysis_model = model,
          Index = index_label,
          Tail = tail,
          Type = c("Unconditional", rep("By ability", n_ability_bins)),
          Ability_bin = c(NA_integer_, seq_len(n_ability_bins)),
          Theta_lower = c(-Inf, bin_edges[seq_len(n_ability_bins)]),
          Theta_upper = c(Inf, bin_edges[-1L]),
          N_boot_scores = c(sum(finite_v), n_boot_bin),
          Normal_cutoff = nominal_cutoff,
          Critical_value = c(cv_all, cv_bin_all),
          Monte_Carlo_SE = c(
            mc_se(finite_v),
            vapply(seq_len(n_ability_bins), function(bin) {
              mc_se(finite_v & boot_bin == bin & !is.na(boot_bin))
            }, numeric(1))
          ),
          Degenerate_bin = c(FALSE, degenerate_bin),
          Fallback_to_unconditional = c(
            FALSE,
            n_boot_bin < min_boot_scores_per_bin | degenerate_bin
          ),
          stringsAsFactors = FALSE
        )

        # Overall flagging table.
        for (method in names(method_labels)) {

          n_valid <- sum(is.finite(observed_values))
          n_flagged <- sum(observed_flags[[method]], na.rm = TRUE)
          null_rate <- flag_rate(v, index_name, null_cutoffs[[method]])

          flagging_rows[[length(flagging_rows) + 1L]] <- data.frame(
            Working_model = working_model,
            Analysis_model = model,
            Index = index_label,
            Method = unname(method_labels[method]),
            N_valid = n_valid,
            N_flagged = n_flagged,
            Percent_flagged = if (n_valid > 0L) 100 * n_flagged / n_valid else NA_real_,
            Null_rate = null_rate,
            Excess = if (n_valid > 0L) n_flagged / n_valid - null_rate else NA_real_,
            stringsAsFactors = FALSE
          )
        }

        # Ability-conditional flagging table (same bins for all methods).
        for (method in names(method_labels)) {
          for (bin in seq_len(n_ability_bins)) {

            use_obs <- !is.na(observed_bin) & observed_bin == bin
            use_boot <- !is.na(boot_bin) & boot_bin == bin
            values_bin <- observed_values[use_obs]
            n_valid <- sum(is.finite(values_bin))
            n_flagged <- sum(observed_flags[[method]][use_obs], na.rm = TRUE)

            ability_rows[[length(ability_rows) + 1L]] <- data.frame(
              Working_model = working_model,
              Analysis_model = model,
              Index = index_label,
              Method = unname(method_labels[method]),
              Ability_bin = bin,
              Theta_lower = bin_edges[bin],
              Theta_upper = bin_edges[bin + 1L],
              N = sum(use_obs),
              Mean_theta = if (any(use_obs)) mean(observed$theta_WL[use_obs]) else NA_real_,
              N_valid = n_valid,
              N_flagged = n_flagged,
              Percent_flagged = if (n_valid > 0L) 100 * n_flagged / n_valid else NA_real_,
              Null_rate = flag_rate(
                v[use_boot],
                index_name,
                null_cutoffs[[method]][use_boot]
              ),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    critical_values <- do.call(rbind, critical_value_rows)
    calibrated_flagging <- do.call(rbind, flagging_rows)
    calibrated_ability_conditional <- do.call(rbind, ability_rows)

    # Keep the model-specific objects in sync with the new columns.
    for (model in successful_models) {
      pfs_by_model[[model]]$person_fit <-
        person_fit_results[person_fit_results$Model == model, ]
    }

    if (!return_bootstrap_scores) {
      bootstrap_scores <- NULL
    }
  }


  # ----------------------------------------------------------
  # 5.12 Save analysis tables
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

      utils::write.csv(
        critical_values,
        file.path(analysis_dir, "bootstrap_critical_values.csv"),
        row.names = FALSE
      )

      utils::write.csv(
        calibrated_flagging,
        file.path(analysis_dir, "calibrated_flagging_summary.csv"),
        row.names = FALSE
      )

      utils::write.csv(
        calibrated_ability_conditional,
        file.path(analysis_dir, "calibrated_ability_conditional_summary.csv"),
        row.names = FALSE
      )
    }
  }


  # ----------------------------------------------------------
  # 5.13 Console summary
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

    if (compute_absolute_fit) {
      cat(
        "Working model absolute fit acceptable (RMSEA <= ", rmsea_cutoff,
        ", SRMSR <= ", srmsr_cutoff, "): ", working_model_acceptable,
        "\n",
        sep = ""
      )
    }

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

      cat("\nFlagging rates by method (normal theory vs bootstrap-calibrated):\n")
      print(calibrated_flagging)
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
      working_model_acceptable_fit = working_model_acceptable,
      model_comparison = model_comparison,
      item_parameters = item_parameters_all_models,
      person_fit = person_fit_results,
      flagging = flag_summary,
      distributions = descriptive_statistics,
      ability_conditional = ability_conditional_summary,
      model_stability = model_stability,
      bootstrap_details = bootstrap_details,
      bootstrap_typeI = bootstrap_typeI,
      bootstrap_analysis_models = bootstrap_models,
      critical_values = critical_values,
      calibrated_flagging = calibrated_flagging,
      calibrated_ability_conditional = calibrated_ability_conditional,
      bootstrap_scores = bootstrap_scores,
      model_objects = ltm_fits,
      pfs_objects = pfs_by_model,
      distribution_plot_file = distribution_plot_file,
      ability_plot_file = ability_plot_file,
      output_dir = analysis_dir
    )
  )
}


