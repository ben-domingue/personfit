#############################################################
# PILOT COMPUTE SCRIPT (v2): PERSON-FIT STATISTICS ACROSS
# 1PL / 2PL / 3PL, ON A PILOT SET OF IRW TABLES
#
# This is a thin IRW harness around the single-dataset engine in
# IRW_personfit_Sept14092026.R (`analyze_irw_personfit()`), which
# is sourced verbatim -- no analysis logic is reimplemented here.
#
# For each table in `poc_datasets`:
#   1. Fetch (irw::irw_fetch) + reshape to a wide 0/1 matrix
#      (irw::irw_long2resp), recode to strict 0/1, drop
#      zero-variance items, and drop incomplete cases. The last
#      step is required because irtoys::est() -- and therefore
#      analyze_irw_personfit() -- does not accept missing data.
#   2. Hand the matrix to analyze_irw_personfit(), which fits
#      1PL, 2PL and 3PL, selects a working model by BIC, checks
#      absolute fit (M2, RMSEA, SRMSR), computes five person-fit
#      indices (lz*, U*, W*, ECI_Z2*, ECI_Z4*) under every model
#      that calibrates, runs a parametric bootstrap from the
#      working model rescored under all three analysis models,
#      and reflags respondents against bootstrap critical values,
#      both unconditional and within ability quintiles.
#   3. Cache the per-dataset result as its own .rds (skipped if
#      it already exists, so reruns don't redo finished tables).
#
# Per-table results are then pooled into one .rds per analysis
# table, all under persondata_v2/. Because process_dataset_v2()
# short-circuits on an existing cache file regardless of code
# changes, delete persondata_v2/*.rds before rerunning after an
# edit to this script or to the engine.
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
})

source("IRW_personfit_Sept14092026.R")

cache_dir <- file.path("persondata_v2")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05
model_selection <- "BIC"
n_ability_bins <- 5L
run_bootstrap <- TRUE
bootstrap_reps <- 200L
bootstrap_analysis_models <- "all"
bootstrap_cores <- max(1L, min(16L, parallel::detectCores() - 2L))
rmsea_cutoff <- 0.05
srmsr_cutoff <- 0.05
seed <- 123

# The 18 tables that calibrated successfully in the July 3PL-only
# pilot. The Sept 10 proof of concept used only the first four listed;
# the full set is used now so that flagging rates can be examined as a
# function of test length and sample size. 200 bootstrap replicates
# are needed for stable quantiles within ability quintiles.
poc_datasets <- c(
  "verbagg", "janssen2_tam", "frac20", "cdm_timss03",
  "blum_2018_imak_bin", "cdm_ecpe", "difnlr_msatb", "experimental_iq",
  "florida_twins_dbi", "g308_sirt", "gcbs_brotherton_2013_vcl",
  "geiser_tam", "gilbert_meta_2", "himmelstein-shipley_abstraction-2025",
  "mcmi_mokken", "mpsycho_wilmer", "pks_probability", "psychtools_ability"
)


# ============================================================
# IRW FETCH + RESHAPE
# ============================================================
# Returns list(X = matrix or NULL, reason = NA/character, prep = list()).
# `prep` records what the preparation steps had to change, since the
# engine's preconditions (complete data, no constant items) are
# stricter than what IRW tables generally satisfy.

fetch_irw_matrix <- function(dataset_name) {

  fail <- function(reason, prep = list()) {
    list(X = NULL, reason = reason, prep = prep)
  }

  fetched <- tryCatch(irw::irw_fetch(dataset_name), error = function(e) e)
  if (inherits(fetched, "error")) {
    return(fail(paste("irw_fetch failed:", conditionMessage(fetched))))
  }

  resp_wide <- tryCatch(irw::irw_long2resp(fetched), error = function(e) e)
  if (inherits(resp_wide, "error")) {
    return(fail(paste("irw_long2resp failed:", conditionMessage(resp_wide))))
  }

  ids <- resp_wide$id
  resp_wide$id <- NULL
  X <- as.matrix(resp_wide)
  storage.mode(X) <- "numeric"
  rownames(X) <- ids

  n_items_raw <- ncol(X)
  n_persons_raw <- nrow(X)
  density_raw <- mean(!is.na(X))

  # Recode to strict 0/1. Two response categories are guaranteed by the
  # IRW filter used to pick these tables, but they may not already be
  # coded 0/1 (e.g. 1/2), and irw_long2resp's default "mean" aggregation
  # of duplicate id-item pairs can introduce a few fractional values.
  raw_vals <- sort(unique(as.vector(X[!is.na(X)])))
  if (length(raw_vals) == 2) {
    Xcmp <- (X == raw_vals[2])
    storage.mode(Xcmp) <- "numeric"
    dimnames(Xcmp) <- dimnames(X)
    X <- Xcmp
  } else if (length(raw_vals) > 2) {
    X <- round(X)
    X[X < 0] <- 0
    X[X > 1] <- 1
  }

  # Drop incomplete cases: irtoys::est() cannot take missing responses.
  complete_rows <- stats::complete.cases(X)
  n_dropped_persons <- sum(!complete_rows)
  X <- X[complete_rows, , drop = FALSE]

  if (nrow(X) < 50) {
    return(fail(
      sprintf("too few complete cases (n_persons=%d)", nrow(X))
    ))
  }

  # Drop constant items, which the engine rejects outright. This has to
  # happen after the complete-case step, since dropping persons can make
  # a previously varying item constant.
  item_sums <- colSums(X)
  keep_items <- item_sums > 0 & item_sums < nrow(X)
  n_dropped_items <- sum(!keep_items)
  X <- X[, keep_items, drop = FALSE]

  if (ncol(X) < 5) {
    return(fail(
      sprintf("too few non-constant items (n_items=%d)", ncol(X))
    ))
  }

  list(
    X = X,
    reason = NA_character_,
    prep = list(
      n_persons_raw = n_persons_raw,
      n_items_raw = n_items_raw,
      density_raw = density_raw,
      n_dropped_persons = n_dropped_persons,
      n_dropped_items = n_dropped_items,
      n_persons = nrow(X),
      n_items = ncol(X)
    )
  )
}


# ============================================================
# PER-DATASET WORKFLOW
# ============================================================

process_dataset_v2 <- function(dataset_name) {

  cache_file <- file.path(cache_dir, paste0(dataset_name, ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  message("=== ", dataset_name, " ===")
  started <- Sys.time()

  result <- list(
    dataset = dataset_name,
    status = "failed",
    reason = NA_character_,
    prep = NULL,
    analysis = NULL,
    elapsed_secs = NA_real_
  )

  prepared <- fetch_irw_matrix(dataset_name)
  result$prep <- prepared$prep

  if (is.null(prepared$X)) {
    result$reason <- prepared$reason
    saveRDS(result, cache_file)
    return(result)
  }

  analysis <- tryCatch(
    analyze_irw_personfit(
      X = prepared$X,
      dataset_name = dataset_name,
      alpha = alpha,
      model_selection = model_selection,
      n_ability_bins = n_ability_bins,
      make_plots = FALSE,
      save_tables = TRUE,
      output_dir = cache_dir,
      run_bootstrap = run_bootstrap,
      bootstrap_reps = bootstrap_reps,
      seed = seed,
      verbose = TRUE,
      compute_absolute_fit = TRUE,
      rmsea_cutoff = rmsea_cutoff,
      srmsr_cutoff = srmsr_cutoff,
      bootstrap_analysis_models = bootstrap_analysis_models,
      bootstrap_cores = bootstrap_cores
    ),
    error = function(e) e
  )

  result$elapsed_secs <- as.numeric(
    difftime(Sys.time(), started, units = "secs")
  )

  if (inherits(analysis, "error")) {
    result$reason <- paste(
      "analyze_irw_personfit failed:",
      conditionMessage(analysis)
    )
    saveRDS(result, cache_file)
    return(result)
  }

  # The fitted ltm objects carry the full response data and are by far
  # the largest part of the returned list; drop them before caching,
  # since every quantity the report needs is already tabulated.
  analysis$model_objects <- NULL
  analysis$pfs_objects <- NULL

  result$status <- "ok"
  result$analysis <- analysis
  saveRDS(result, cache_file)
  result
}


# ============================================================
# RUN AND POOL
# ============================================================

all_results <- lapply(poc_datasets, process_dataset_v2)
names(all_results) <- poc_datasets

ok_results <- Filter(function(r) identical(r$status, "ok"), all_results)

# Stack one analysis table across datasets, tagging each row with its
# source dataset. Returns a zero-row data frame if nothing succeeded.
pool_table <- function(component) {
  rows <- lapply(names(ok_results), function(name) {
    tab <- ok_results[[name]]$analysis[[component]]
    if (is.null(tab) || nrow(tab) == 0L) return(NULL)
    cbind(dataset = name, tab, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(data.frame())
  # Tables cached by an earlier engine revision can lack columns added
  # later (e.g. critical_values$Degenerate_bin); pad those with NA.
  all_cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(tab) {
    for (col in setdiff(all_cols, names(tab))) tab[[col]] <- NA
    tab[, all_cols]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

pooled <- list(
  model_comparison    = pool_table("model_comparison"),
  item_parameters     = pool_table("item_parameters"),
  person_fit          = pool_table("person_fit"),
  flagging            = pool_table("flagging"),
  distributions       = pool_table("distributions"),
  ability_conditional = pool_table("ability_conditional"),
  model_stability     = pool_table("model_stability"),
  bootstrap_typeI     = pool_table("bootstrap_typeI"),
  bootstrap_details   = pool_table("bootstrap_details"),
  critical_values     = pool_table("critical_values"),
  calibrated_flagging = pool_table("calibrated_flagging"),
  calibrated_ability_conditional = pool_table("calibrated_ability_conditional")
)

# One row per attempted table: dimensions, what preparation changed,
# which model won, and runtime.
manifest <- do.call(
  rbind,
  lapply(all_results, function(r) {
    prep <- r$prep
    data.frame(
      dataset = r$dataset,
      status = r$status,
      reason = ifelse(is.na(r$reason), "", r$reason),
      n_persons_raw = if (is.null(prep)) NA_integer_ else prep$n_persons_raw,
      n_items_raw = if (is.null(prep)) NA_integer_ else prep$n_items_raw,
      density_raw = if (is.null(prep)) NA_real_ else prep$density_raw,
      n_dropped_persons = if (is.null(prep)) NA_integer_ else prep$n_dropped_persons,
      n_dropped_items = if (is.null(prep)) NA_integer_ else prep$n_dropped_items,
      n_persons = if (is.null(prep)) NA_integer_ else prep$n_persons,
      n_items = if (is.null(prep)) NA_integer_ else prep$n_items,
      working_model = if (is.null(r$analysis)) NA_character_ else r$analysis$working_model,
      working_RMSEA = if (is.null(r$analysis)) NA_real_ else
        with(r$analysis$model_comparison, RMSEA[Working_model]),
      working_SRMSR = if (is.null(r$analysis)) NA_real_ else
        with(r$analysis$model_comparison, SRMSR[Working_model]),
      working_acceptable_fit = if (is.null(r$analysis)) NA else
        isTRUE(r$analysis$working_model_acceptable_fit),
      elapsed_secs = r$elapsed_secs,
      stringsAsFactors = FALSE
    )
  })
)
rownames(manifest) <- NULL

settings <- list(
  alpha = alpha,
  model_selection = model_selection,
  n_ability_bins = n_ability_bins,
  run_bootstrap = run_bootstrap,
  bootstrap_reps = bootstrap_reps,
  bootstrap_analysis_models = bootstrap_analysis_models,
  bootstrap_cores = bootstrap_cores,
  rmsea_cutoff = rmsea_cutoff,
  srmsr_cutoff = srmsr_cutoff,
  seed = seed,
  poc_datasets = poc_datasets,
  engine_script = "IRW_personfit_Sept14092026.R"
)

saveRDS(pooled, file.path(cache_dir, "pooled_results.rds"))
saveRDS(manifest, file.path(cache_dir, "manifest.rds"))
saveRDS(settings, file.path(cache_dir, "settings.rds"))

cat("\n\n================ PILOT SUMMARY ================\n")
print(manifest)
cat("\nPooled table row counts:\n")
print(vapply(pooled, nrow, integer(1)))
