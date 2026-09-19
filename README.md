# Person fit on IRW data

This repo studies how person-fit statistics behave on real item response data from the [Item Response Warehouse](https://itemresponsewarehouse.org/) (IRW). For each IRW table it fits 1PL, 2PL and 3PL models, picks a *working model* by BIC, and computes five person-fit indices (lz\*, U\*, W\*, ECI_Z2\*, ECI_Z4\*) under every model. It then runs a parametric bootstrap to see how often each index flags respondents when the model is true, and reflags respondents against bootstrap critical values, both overall and within ability quintiles.

The current results are in `person_fit_v2.html`, a self-contained page you can open in a browser. Discussion of the analysis happens in [GitHub issues](https://github.com/ben-domingue/personfit/issues).

## What's in the repo

| File | What it is |
|---|---|
| `IRW_personfit_Sept14092026.R` | **The engine.** `analyze_irw_personfit(X, ...)` takes one complete 0/1 response matrix and does all of the analysis. It knows nothing about IRW. |
| `person_fit_v2_compute.R` | **The harness.** Fetches each IRW table, turns it into a matrix the engine accepts, calls the engine, and caches and pools the results in `persondata_v2/`. |
| `person_fit_v2.qmd` | **The report.** Reads only the cached results and renders `person_fit_v2.html`. It also writes `persondata_v2/flag_results.csv`. |
| `persondata_v2/flag_results.csv` | Tidy flagging results with one row per table × analysis model × statistic × method. This is the only computed output kept in git. |
| `IRW_personfit_Sept09092026.R`, `archive_v2_sept10/` | The previous engine and the 10 September report, kept for comparison. |
| `person_fit_compute.R`, `person_fit.qmd`, `person_fit.html` | The July 3PL-only pilot. Superseded by the current pipeline. |
| `Rcode.txt`, `Rcode with some updates_july272026.txt` | Original single-dataset scratch scripts the July pilot was built from. |

## Getting started

### 1. Install software

You need R (≥ 4.3; developed on 4.6.1) and [Quarto](https://quarto.org/docs/get-started/) (developed on 1.8) to render the report.

```r
install.packages(c("ltm", "irtoys", "PerFit", "aberrance", "mirt",
                   "dplyr", "tidyr", "ggplot2", "knitr", "scales",
                   "devtools"))
devtools::install_github("itemresponsewarehouse/Rpkg")   # the irw package
```

### 2. Get access to IRW data

IRW data is hosted on [Redivis](https://redivis.com/), so you need a free Redivis account. The first time you fetch a table, `irw` opens a browser window so you can log in. After that the credential is cached and later fetches run without asking. To check that access works:

```r
library(irw)
df <- irw_fetch("verbagg")   # opens a Redivis login on first use
head(df)
```

The data is in **long format**, with one row per person × item response: `id` (person), `item`, `resp`, plus other columns in some tables. `irw_long2resp(df)` reshapes it into a wide person × item matrix, which is what the harness does. The [IRW website](https://itemresponsewarehouse.org/) has a getting-started guide and a searchable catalogue of tables.

On a server or other machine without a browser, create an API token in your Redivis account settings and set the `REDIVIS_API_TOKEN` environment variable instead (for example in `~/.Renviron`, which is gitignored here).

### 3. Look at the existing results first

Open `person_fit_v2.html`. It has the method, all figures and tables, and the code that made them. `persondata_v2/flag_results.csv` has the same flagging results as a flat file.

### 4. Rerun the analysis

Run everything from the repo root, because paths are relative:

```sh
Rscript person_fit_v2_compute.R    # fetch + analyse all tables (~1.5 h, see below)
quarto render person_fit_v2.qmd    # rebuild the report from the cache
```

Keep these points in mind:

- **The cache is sticky.** `process_dataset_v2()` returns `persondata_v2/<table>.rds` if the file already exists, even if you changed the code. After you edit the engine or the harness settings, delete the caches first. Keep `flag_results.csv`, because it is tracked in git and the report rewrites it anyway:
  ```sh
  rm persondata_v2/*.rds
  ```
  To rerun just one table, delete only that table's `.rds`.
- **Runtime.** The last full run (18 tables, 200 bootstrap replicates, 16 cores) took about 86 minutes. Most of that is the bootstrap, and the time grows with sample size: `cdm_ecpe` (n = 2,922) alone took about 18 minutes. For a quick test, set `bootstrap_reps <- 20L` or trim `poc_datasets` to a few small tables such as `verbagg` or `pks_probability`.
- **Parallelism** uses `parallel::mclapply`, so it only works on macOS and Linux. On Windows the bootstrap runs on one core and will be much slower.
- The report only reads caches. After changing the `.qmd` you can re-render without recomputing anything.

### 5. Analyse different or more IRW tables

Edit `poc_datasets` near the top of `person_fit_v2_compute.R`. To find candidate tables:

```r
library(irw)
irw_filter(n_categories = 2)                  # dichotomous tables (~475)
irw_filter(n_categories = 2, n_items = c(10, 60), n_participants = c(200, 5000))
irw_info("verbagg")                           # metadata for one table
```

The harness prepares each table the same way. It recodes responses to strict 0/1, drops **incomplete cases** (because `irtoys::est()` cannot handle missing data), drops constant items, and skips tables with fewer than 50 complete persons or fewer than 5 items. `manifest.rds` records how much each step removed, as `n_dropped_persons` and `n_dropped_items`. Check it, because tables with sparse designs can lose most of their respondents.

To run the engine on a single matrix outside the harness:

```r
source("IRW_personfit_Sept14092026.R")
res <- analyze_irw_personfit(X, dataset_name = "mytable",
                             bootstrap_reps = 50, make_plots = FALSE)
names(res)   # model_comparison, person_fit, flagging, bootstrap_typeI, ...
```

## Outputs in `persondata_v2/` (not in git, regenerated)

| File | Contents |
|---|---|
| `<table>.rds` | Full result for one table: `status`, `reason` if it failed, `prep` (what data prep removed) and `analysis` (the engine's output). |
| `pooled_results.rds` | Each engine table stacked across datasets, with a `dataset` column added. |
| `manifest.rds` | One row per attempted table: status, size before and after prep, working model, RMSEA/SRMSR and runtime. |
| `settings.rds` | The alpha, number of bootstrap replicates and other settings used for the run. |
| `<table>_personfit/` | Per-table CSVs written by the engine. |

## Known issues and open questions

- **Two of the 18 tables fail** (`florida_twins_dbi`, `gcbs_brotherton_2013_vcl`). BIC picks the 3PL for both, and the person-fit calibration then fails with a singular fit. There is no fallback to `mirt` yet.
- **Complete-case deletion is not neutral.** It rules out tables with planned missingness and changes the sample in others.
- `mcmi_mokken` and `pks_probability` have results for only 2 of the 3 analysis models in `flag_results.csv`. The cause has not been investigated yet.
- The working model is chosen by BIC, and AIC often disagrees. The choice feeds into the bootstrap reference distribution.
- Ht (in the July pilot) was dropped from the current engine. Whether to bring it back is still open.
