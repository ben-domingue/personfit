#############################################################
# PILOT COMPUTE SCRIPT: PERSON-FIT STATISTICS UNDER THE 3PL,
# ACROSS A PILOT SET OF REAL IRW DATASETS
#
# For each dataset in `pilot_datasets`:
#   1. Fetch (irw::irw_fetch) + reshape to wide 0/1 matrix
#      (irw::irw_long2resp), dropping zero-variance items.
#   2. Fit a 3PL: irtoys::est(..., engine = "ltm") first,
#      falling back to mirt::mirt(itemtype = "3PL") if the
#      ltm fit fails or returns non-finite/degenerate params.
#      Skip the dataset (with a logged reason) if both fail.
#   3. Estimate theta via irtoys::wle().
#   4. Compute six person-fit indices (lz*, U*, W*, ECI_Z2*,
#      ECI_Z4*, Ht) and flag at alpha = 0.05.
#   5. Additionally fit a 2PL (ltm::ltm) and 3PL (ltm::tpm) on
#      the same item matrix and compare them by AIC/BIC, as a
#      diagnostic on whether the 3PL is actually preferred over
#      the simpler 2PL for each dataset. This is independent of
#      the ltm/mirt engine used for the person-fit calibration
#      above and does not gate the rest of the pipeline.
#   6. Cache the per-dataset result as its own .rds (skipped
#      if it already exists, so reruns after a partial failure
#      don't redo completed datasets).
#
# All per-dataset results are combined into a single long data
# frame and cached as persondata/person_fit_results.rds. A log
# of skipped datasets (name + reason) is cached as
# persondata/skipped_datasets.rds. Per-dataset 2PL vs 3PL model
# comparisons are cached as persondata/model_comparison_results.rds,
# and the alpha level used for flagging as persondata/alpha.rds.
#############################################################

# ============================================================
# SETUP
# ============================================================

suppressMessages({
  library(irw)
  library(ltm)
  library(irtoys)
  library(PerFit)
  library(aberrance)
  library(mirt)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(future)
})

cache_dir <- file.path("persondata")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05
normal_cutoff <- qnorm(alpha)

# Pilot dataset selection: irw_filter(n_categories = 2, n_items = c(15, 60),
# n_participants = c(300, 5000), density = c(0.8, 1.0)) returned 65 candidate
# dichotomous tables. The 24 below were hand-picked from that list for
# diversity of construct/source while avoiding near-duplicate families
# (e.g. only one of the many "gilbert_meta_*" simulation-adjacent tables).
pilot_datasets <- c(
  "balance_mokken",
  "bicb-j_ishiguro_2025",
  "blum_2018_imak_bin",
  "cdm_ecpe",
  "cdm_timss03",
  "difnlr_msatb",
  "dscore_barrera_weber_2019",
  "erf_breuer_2017_frmfr",
  "experimental_iq",
  "FACIT_YOUNT_2021_clinic10",
  "florida_twins_dbi",
  "frac20",
  "g308_sirt",
  "gcbs_brotherton_2013_vcl",
  "geiser_tam",
  "gilbert_meta_2",
  "himmelstein-shipley_abstraction-2025",
  "janssen2_tam",
  "KoreanNursing_Park_2017",
  "mcmi_mokken",
  "mpsycho_wilmer",
  "pks_probability",
  "psychtools_ability",
  "verbagg"
)


# ============================================================
# PERSON-FIT HELPER FUNCTIONS (U*/W* corrected indices)
# Verbatim from the supplied single-dataset example code.
# ============================================================

Pi_3PL <- function(theta, item_parameters) {

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


ri_3PL <- function(theta, item_parameters) {

  values <- Pi_3PL(theta, item_parameters)

  values$dPi /
    (values$Pi * (1 - values$Pi))
}


r0_WL <- function(theta, item_parameters) {

  values <- Pi_3PL(theta, item_parameters)
  r_values <- ri_3PL(theta, item_parameters)

  denominator <- 2 * sum(r_values * values$dPi)

  if (!is.finite(denominator) || abs(denominator) < 1e-12) {
    return(NA_real_)
  }

  sum(r_values * values$d2Pi) / denominator
}


item_weights <- function(theta, item_parameters, index) {

  index <- match.arg(index, c("lz", "U", "W"))

  probability <- Pi_3PL(theta, item_parameters)$Pi
  J <- nrow(item_parameters)

  switch(
    index,

    lz = log(probability / (1 - probability)),

    U = (1 - 2 * probability) /
      (J * probability * (1 - probability)),

    W = (1 - 2 * probability) /
      sum(probability * (1 - probability))
  )
}


Wn <- function(response_vector, theta, item_parameters, index) {

  probability <- Pi_3PL(theta, item_parameters)$Pi

  weights <- item_weights(
    theta = theta,
    item_parameters = item_parameters,
    index = index
  )

  sum((response_vector - probability) * weights)
}


cn <- function(theta, item_parameters, index) {

  values <- Pi_3PL(theta, item_parameters)

  weights <- item_weights(
    theta = theta,
    item_parameters = item_parameters,
    index = index
  )

  numerator <- sum(values$dPi * weights)
  denominator <- sum(values$dPi * ri_3PL(theta, item_parameters))

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
    correction * ri_3PL(theta, item_parameters)
}


tau2n <- function(theta, item_parameters, index) {

  probability <- Pi_3PL(theta, item_parameters)$Pi

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
# 3PL CALIBRATION: ltm (via irtoys) WITH mirt FALLBACK
# ============================================================

# Sanity-check an (a, b, c) item-parameter matrix. Returns TRUE only if it
# is usable for the person-fit math downstream (finite, positive
# discrimination, guessing strictly inside [0, 1)).
is_valid_IP <- function(IP) {
  if (is.null(IP) || !is.matrix(IP) || ncol(IP) != 3L) {
    return(FALSE)
  }
  if (anyNA(IP) || any(!is.finite(IP))) {
    return(FALSE)
  }
  if (any(IP[, 1] <= 0)) {
    return(FALSE)
  }
  if (any(IP[, 3] < 0) || any(IP[, 3] >= 1)) {
    return(FALSE)
  }
  TRUE
}

# Cap guessing parameters away from the degenerate boundaries that break
# the (1 - c)-type terms in the person-fit math, and floor discrimination
# just above zero.
regularize_IP <- function(IP, c_max = 0.97, a_min = 1e-3) {
  IP[, 1] <- pmax(IP[, 1], a_min)
  IP[, 3] <- pmin(pmax(IP[, 3], 0), c_max)
  IP
}

fit_3PL_ltm <- function(X) {
  fit <- irtoys::est(resp = X, model = "3PL", engine = "ltm")
  IP <- as.matrix(fit$est)
  if (ncol(IP) != 3L) {
    stop("ltm/irtoys item-parameter matrix does not have three columns.")
  }
  colnames(IP) <- c("a", "b", "c")
  rownames(IP) <- colnames(X)
  IP
}

fit_3PL_mirt <- function(X) {
  fit <- mirt::mirt(as.data.frame(X), 1, itemtype = "3PL", verbose = FALSE)
  co <- mirt::coef(fit, IRTpars = TRUE, simplify = TRUE)$items
  IP <- co[, c("a", "b", "g"), drop = FALSE]
  colnames(IP) <- c("a", "b", "c")
  rownames(IP) <- colnames(X)
  as.matrix(IP)
}

# Returns list(IP = matrix or NULL, engine = "ltm"/"mirt"/NA, reason = NULL/character)
fit_3PL_with_fallback <- function(X) {

  ltm_result <- tryCatch(fit_3PL_ltm(X), error = function(e) e)

  if (!inherits(ltm_result, "error") && is_valid_IP(ltm_result)) {
    return(list(IP = regularize_IP(ltm_result), engine = "ltm", reason = NA_character_))
  }

  ltm_reason <- if (inherits(ltm_result, "error")) {
    paste("ltm error:", conditionMessage(ltm_result))
  } else {
    "ltm returned non-finite/degenerate item parameters"
  }

  mirt_result <- tryCatch(fit_3PL_mirt(X), error = function(e) e)

  if (!inherits(mirt_result, "error") && is_valid_IP(mirt_result)) {
    return(list(
      IP = regularize_IP(mirt_result),
      engine = "mirt",
      reason = paste0("ltm fallback triggered (", ltm_reason, ")")
    ))
  }

  mirt_reason <- if (inherits(mirt_result, "error")) {
    paste("mirt error:", conditionMessage(mirt_result))
  } else {
    "mirt returned non-finite/degenerate item parameters"
  }

  list(
    IP = NULL,
    engine = NA_character_,
    reason = paste0("both engines failed -- ltm: ", ltm_reason, " | mirt: ", mirt_reason)
  )
}


# ============================================================
# 2PL vs 3PL MODEL COMPARISON (AIC/BIC) -- DIAGNOSTIC ONLY
# ============================================================
# Relative-fit comparison, independent of the ltm/mirt engine used for
# the person-fit calibration above. Smaller AIC/BIC indicates the
# preferred model. No likelihood-ratio test is used since the 2PL sits
# on a boundary of the 3PL (guessing parameters = 0).

compare_2PL_3PL <- function(X) {

  X_data <- as.data.frame(X)

  fit_2PL_ltm <- ltm::ltm(X_data ~ z1, IRT.param = TRUE)
  fit_3PL_ltm <- ltm::tpm(X_data, type = "latent.trait", IRT.param = TRUE)

  model_comparison <- data.frame(
    Model = c("2PL", "3PL"),
    Log_likelihood = c(as.numeric(logLik(fit_2PL_ltm)), as.numeric(logLik(fit_3PL_ltm))),
    N_parameters = c(attr(logLik(fit_2PL_ltm), "df"), attr(logLik(fit_3PL_ltm), "df")),
    AIC = c(AIC(fit_2PL_ltm), AIC(fit_3PL_ltm)),
    BIC = c(BIC(fit_2PL_ltm), BIC(fit_3PL_ltm)),
    stringsAsFactors = FALSE
  )

  model_comparison$Delta_AIC <- model_comparison$AIC - min(model_comparison$AIC)
  model_comparison$Delta_BIC <- model_comparison$BIC - min(model_comparison$BIC)
  model_comparison$Preferred_by_AIC <- model_comparison$AIC == min(model_comparison$AIC)
  model_comparison$Preferred_by_BIC <- model_comparison$BIC == min(model_comparison$BIC)

  model_comparison
}


# ============================================================
# PER-DATASET WORKFLOW
# ============================================================

process_dataset <- function(dataset_name) {

  cache_file <- file.path(cache_dir, paste0(dataset_name, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  result <- list(
    dataset = dataset_name,
    status = "failed",
    reason = NA_character_,
    engine = NA_character_,
    data = NULL
  )

  # ---- 1. fetch + reshape --------------------------------------------
  fetched <- tryCatch(irw::irw_fetch(dataset_name), error = function(e) e)
  if (inherits(fetched, "error")) {
    result$reason <- paste("irw_fetch failed:", conditionMessage(fetched))
    saveRDS(result, cache_file)
    return(result)
  }

  resp_wide <- tryCatch(irw::irw_long2resp(fetched), error = function(e) e)
  if (inherits(resp_wide, "error")) {
    result$reason <- paste("irw_long2resp failed:", conditionMessage(resp_wide))
    saveRDS(result, cache_file)
    return(result)
  }

  ids <- resp_wide$id
  resp_wide$id <- NULL
  X <- as.matrix(resp_wide)
  storage.mode(X) <- "numeric"
  rownames(X) <- ids

  # Recode to strict 0/1. irw_filter(n_categories = 2) guarantees two
  # response categories in the source data, but they may not already be
  # coded 0/1 (e.g. 1/2), and irw_long2resp's default "mean" aggregation
  # of duplicate id-item pairs can introduce a few fractional values.
  raw_vals <- sort(unique(as.vector(X[!is.na(X)])))
  if (length(raw_vals) == 2) {
    Xcmp <- (X == raw_vals[2])
    storage.mode(Xcmp) <- "numeric"
    X <- Xcmp
  } else if (length(raw_vals) > 2) {
    X <- round(X)
    X[X < 0] <- 0
    X[X > 1] <- 1
  }

  # drop zero-variance items (needed for stable calibration)
  item_var <- apply(X, 2, function(z) var(z, na.rm = TRUE))
  keep_items <- !is.na(item_var) & item_var > 0
  X <- X[, keep_items, drop = FALSE]

  if (ncol(X) < 5 || nrow(X) < 50) {
    result$reason <- sprintf(
      "too few usable items/persons after dropping zero-variance items (n_items=%d, n_persons=%d)",
      ncol(X), nrow(X)
    )
    saveRDS(result, cache_file)
    return(result)
  }

  n_items <- ncol(X)
  n_participants <- nrow(X)
  density <- mean(!is.na(X))

  # ---- 1b. 2PL vs 3PL model comparison (diagnostic only) ---------------
  model_cmp <- tryCatch(compare_2PL_3PL(X), error = function(e) e)
  if (inherits(model_cmp, "error")) {
    result$model_comparison_reason <- paste("2PL/3PL comparison failed:", conditionMessage(model_cmp))
    result$model_comparison <- NULL
  } else {
    model_cmp$dataset <- dataset_name
    result$model_comparison <- model_cmp
    result$model_comparison_reason <- NA_character_
  }

  # ---- 2. calibrate 3PL (ltm, falling back to mirt) -------------------
  cal <- fit_3PL_with_fallback(X)

  if (is.null(cal$IP)) {
    result$reason <- cal$reason
    result$n_items <- n_items
    result$n_participants <- n_participants
    result$density <- density
    saveRDS(result, cache_file)
    return(result)
  }

  IP_3PL <- cal$IP

  # ---- 3. theta via WLE -------------------------------------------------
  theta_object <- tryCatch(
    irtoys::wle(resp = X, ip = IP_3PL),
    error = function(e) e
  )
  if (inherits(theta_object, "error")) {
    result$reason <- paste("irtoys::wle failed:", conditionMessage(theta_object))
    result$n_items <- n_items
    result$n_participants <- n_participants
    result$density <- density
    saveRDS(result, cache_file)
    return(result)
  }

  theta_WL <- as.numeric(theta_object[, 1])
  xi_3PL <- matrix(theta_WL, ncol = 1, dimnames = list(rownames(X), "theta"))

  # ---- 4. person-fit indices --------------------------------------------

  pf <- tryCatch({

    result_lzstar <- PerFit::lzstar(
      matrix = X, IP = IP_3PL, IRT.PModel = "3PL",
      Ability = theta_WL, Ability.PModel = "WL"
    )
    lzstar_scores <- as.numeric(unlist(result_lzstar$PFscores))

    result_Ht <- PerFit::Ht(matrix = X)
    Ht_scores <- as.numeric(unlist(result_Ht$PFscores))

    result_ECI <- aberrance::detect_pm(
      method = c("ECI2_S_TS", "ECI4_S_TS"),
      psi = IP_3PL, xi = xi_3PL, x = X, alpha = alpha
    )
    ECI_Z2_star <- as.numeric(result_ECI$stat[, "ECI2_S_TS"])
    ECI_Z4_star <- as.numeric(result_ECI$stat[, "ECI4_S_TS"])

    U_star_scores <- compute_corrected_index(
      data = X, item_parameters = IP_3PL, theta = theta_WL, index = "U"
    )
    W_star_scores <- compute_corrected_index(
      data = X, item_parameters = IP_3PL, theta = theta_WL, index = "W"
    )

    list(
      lzstar = lzstar_scores,
      Ht = Ht_scores,
      ECI_Z2_star = ECI_Z2_star,
      ECI_Z4_star = ECI_Z4_star,
      U_star = U_star_scores,
      W_star = W_star_scores,
      Ht_object = result_Ht
    )

  }, error = function(e) e)

  if (inherits(pf, "error")) {
    result$reason <- paste("person-fit index computation failed:", conditionMessage(pf))
    result$n_items <- n_items
    result$n_participants <- n_participants
    result$density <- density
    saveRDS(result, cache_file)
    return(result)
  }

  person_fit_results <- data.frame(
    dataset = dataset_name,
    Person = rownames(X),
    theta_WL = theta_WL,
    lzstar = pf$lzstar,
    U_star = pf$U_star,
    W_star = pf$W_star,
    ECI_Z2_star = pf$ECI_Z2_star,
    ECI_Z4_star = pf$ECI_Z4_star,
    Ht = pf$Ht,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  person_fit_results$lzstar_flag <- person_fit_results$lzstar < normal_cutoff
  person_fit_results$U_star_flag <- person_fit_results$U_star < normal_cutoff
  person_fit_results$W_star_flag <- person_fit_results$W_star < normal_cutoff
  person_fit_results$ECI_Z2_star_flag <- person_fit_results$ECI_Z2_star < normal_cutoff
  person_fit_results$ECI_Z4_star_flag <- person_fit_results$ECI_Z4_star < normal_cutoff

  # Ht uses a simulated empirical cutoff, not the normal-theory cutoff.
  Ht_cutoff <- tryCatch({
    set.seed(123)
    Ht_cutoff_object <- PerFit::cutoff(
      x = pf$Ht_object, ModelFit = "NonParametric",
      Nreps = 1000, Blvl = alpha, Breps = 1000
    )
    as.numeric(Ht_cutoff_object$Cutoff)[1]
  }, error = function(e) NA_real_)

  person_fit_results$Ht_flag <- if (is.finite(Ht_cutoff)) {
    person_fit_results$Ht <= Ht_cutoff
  } else {
    NA
  }

  person_fit_results$n_items <- n_items
  person_fit_results$n_participants <- n_participants
  person_fit_results$density <- density
  person_fit_results$engine <- cal$engine

  result$status <- "ok"
  result$reason <- cal$reason
  result$engine <- cal$engine
  result$n_items <- n_items
  result$n_participants <- n_participants
  result$density <- density
  result$Ht_cutoff <- Ht_cutoff
  result$data <- person_fit_results

  saveRDS(result, cache_file)
  result
}


# ============================================================
# RUN PILOT: PARALLEL OUTER LOOP, SINGLE-THREADED 3PL FITS
# ============================================================

n_workers <- max(1, min(4, future::availableCores() - 1))
future::plan(future::multisession, workers = n_workers)

set.seed(123)

all_results <- furrr::future_map(
  pilot_datasets,
  process_dataset,
  .options = furrr::furrr_options(seed = TRUE),
  .progress = TRUE
)

future::plan(future::sequential)

names(all_results) <- pilot_datasets


# ============================================================
# COMBINE + LOG SKIPPED DATASETS
# ============================================================

ok_results <- purrr::keep(all_results, ~ .x$status == "ok")
failed_results <- purrr::keep(all_results, ~ .x$status != "ok")

person_fit_results <- dplyr::bind_rows(purrr::map(ok_results, "data"))

skipped_datasets <- purrr::map_dfr(failed_results, function(r) {
  data.frame(
    dataset = r$dataset,
    reason = r$reason,
    n_items = if (!is.null(r$n_items)) r$n_items else NA_integer_,
    n_participants = if (!is.null(r$n_participants)) r$n_participants else NA_integer_,
    density = if (!is.null(r$density)) r$density else NA_real_,
    stringsAsFactors = FALSE
  )
})

model_comparison_results <- dplyr::bind_rows(purrr::map(all_results, "model_comparison"))

saveRDS(person_fit_results, file.path(cache_dir, "person_fit_results.rds"))
saveRDS(skipped_datasets, file.path(cache_dir, "skipped_datasets.rds"))
saveRDS(model_comparison_results, file.path(cache_dir, "model_comparison_results.rds"))
saveRDS(alpha, file.path(cache_dir, "alpha.rds"))

cat(sprintf(
  "\nDone. %d datasets succeeded, %d skipped.\n",
  length(ok_results), length(failed_results)
))
if (nrow(skipped_datasets) > 0) {
  print(skipped_datasets[, c("dataset", "reason")])
}
