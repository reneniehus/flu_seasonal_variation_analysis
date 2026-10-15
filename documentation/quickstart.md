# Quickstart — kicking off a project

A clone-and-go template for loading ERVISS / RespiCompass respiratory-surveillance data
(plus contact, vaccine and demography data) and turning it into tidy, model-ready tables
with a quality/dynamics report. See `data_overview.md` for what data is available.

## 1. Prerequisites
- **R 4.3+** (the renv lockfile pins R 4.3.3 and 149 packages).
- System libraries for two optional features: **pandoc** (to render the HTML report) and
  **libmagick++** (pulled in by `summarytools`). On Debian/Ubuntu:
  `apt-get install pandoc libmagick++-dev`.

## 2. Install dependencies
Open the project (`flu_seasonal_variation_analysis.Rproj`) and restore the pinned package versions:
```r
renv::restore()
```

## 3. Run the pipeline
Open `code/00_main.R` and run it top to bottom. In order it:
1. sources `setup.R` (libraries + helpers) and `settings_version0.R` (`params`),
2. `load_data(params, regenerate, new_from_online)` -> the `data` list,
3. `gen_model_input(params, data)` -> `models_in` (long / wide / season-summary / contacts),
4. `eyeballing(models_in, params, data)` -> a figure manifest for the report.

Two flags drive data loading:
- `regenerate` — `F` reuse cached `output/*.Rdata`; `T` rebuild them.
- `new_from_online` — `T` fetch from the internet and refresh the local `data/` snapshots;
  `F` use the snapshots (self-bootstrapping: fetches anything missing). First offline run
  works because the snapshots are committed.

## 4. Change what it does — settings
Edit `code/02_settings/settings_version0.R` (the `params` list — live knobs only; every entry has a
consumer in `code/`). Common knobs:
- `susc_*` — the fixed SIR assumptions (R0, infectious period, seed) + method-demo countries.
- `season_start_monthday` / `season_end_monthday` — the Aug→Jul season window.
- `respicompass_round` — RespiCompass hub round folder (bump for a new season; also stamps the
  vaccination-scenario season).

Demography comes from the committed snapshot (`data/population_pyramid.fst`); keep
`use_ecdc_db = FALSE`.

To version a different configuration, copy the file to `settings_version1.R` and source that.

## 5. The eyeballing report
Render the large-format quality + dynamics HTML report:
```r
rmarkdown::render("code/03_report/eyeballing_report.Rmd")
```
Figures and their narrative bullets are defined in `code/01_main_supporting/eyeballing.R`
(edit the bullet placeholders to match what you see). The rendered `.html` is gitignored.

## 6. Add a new data stream
- **Another ERVISS file:** add one row to the `erviss_registry` tribble in `load_data.R`
  (name, file, snapshot, schema). The loop and standardisation handle the rest.
- **A non-ERVISS source:** write a `load_data_<x>()` builder following the existing pattern
  (uses `load_or_build()` for caching) and call it from the `load_data()` mother function.
- To surface it in `models_in`, add a small extractor in `gen_model_input.R` and bind it into
  `make_data_timeseries_long()`.

## 7. Tests / shielding from upstream changes
Lightweight contracts in `code/01_main_supporting/validate.R` fail loudly if ECDC renames or
drops a required column / stream. Run the suite (offline, from the committed snapshots):
```r
Rscript run_tests.R        # or: testthat::test_dir("tests/testthat")
```
It checks the data contracts, the canonical tables, and the key invariant
(`ILI+ == ILI × positivity`, and agreement with the RespiCompass ILI+).

## 8. Modelling
The active modelling is the swappable **method framework** plus a reference Stan model (both below).
The data layer's contact-matrix transform lives in `flu_functions.R`.

Two worked SIR implementations of the same flu ILI+ model are included:

- **`stan/SIR_multiseason_age_vax_2.stan`** — the full age- and vaccination-structured,
  multi-season Bayesian SIR (HMC via rstan), with scenario projections. Note its priors are
  currently commented out; re-enable them before production use.
- **`code/01_main_supporting/` method framework** — a single-population R version of the same
  generative story (SIR → ILI+, neg-binomial-like noise) as a swappable set of methods: a shared
  engine (`sir_core.R`), one file per method (`methods/method_sir_*.R`) and a registry
  (`methods_registry.R`). Each method fixes R0 / infectious period / seed and reads a per-season
  **susceptibility** `S0` and **reporting fraction** `c` off the ILI+ waves. Base R. Fit a country:
  ```r
  sl  <- load_flu_iliplus_slim("DK")               # committed slim panel (no pipeline needed)
  fit <- run_method("deterministic", sl, params)   # any registered method, same call
  summarise_method_fit(fit)                         # S0, R_eff, c, peak/onset week, cor, ...
  ```
  Demo (every method, figure per method): `Rscript code/04_modelling/fit_methods_demo.R`.
  Tested in `tests/testthat/test-sir-deterministic.R` and `test-methods-registry.R`.

## Layout
```
code/00_main.R                 orchestrator (run this)
code/01_main_supporting/       setup, validate, load_data, gen_model_input, flu_functions, eyeballing,
                               sir_core, methods/, methods_registry, send_report
code/02_settings/              settings_version0.R (params)
code/03_report/                eyeballing_report.Rmd, data_availability.R
code/04_modelling/             build_slim_panel.R, fit_methods_demo.R, descriptive_overview.R, ekf_overview.R
code/05_analysis/              driver analysis (prepare_descriptors, analyse_patterns,
                               hierarchical_models, dominant_subtype, bayes_*, plot_*)
code/08_flu_hospitalisation/   influenza-hospitalisation indicators (build_all.sh -> data/flu_hospitalisation_indicators/)
stan/                          SIR_multiseason_age_vax_2.stan (Bayesian SIR)
data/                          committed raw snapshots (offline bootstrap); one folder per data need,
                               each with a README (index: data/README.md)
output/                        cached data lists (gitignored, regenerated)
tests/testthat/                contract + invariant tests
documentation/                 quickstart, data_overview, decisions, analysis_strategy,
                               findings_descriptors, documentation.Rmd
```
