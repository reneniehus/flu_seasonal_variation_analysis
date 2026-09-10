# flu_seasonal_variation_analysis

Characterising the **drivers of inter-seasonal variability** of seasonal influenza across the
EU/EEA — 166 country-seasons of real surveillance data turned into per-season wave features, then
related to candidate drivers (dominant subtype, vaccine effectiveness, 65+ coverage, prior-season
burden). Research aim and scope: `PROJECT_SCOPE.md`.

> Renamed from `flu-susceptible-fit` (2026-08) to match the project purpose. Two related
> repositories share history with this one:
>
> - **[respi_starter](https://github.com/reneniehus/respi_starter)** — the reusable ECDC
>   surveillance ingest template this project grew out of: real ERVISS / RespiCompass data in,
>   canonical model-ready tables plus a data-quality report out, with modelling deliberately
>   left as a scaffold.
> - **[pandemic_toolbook](https://github.com/reneniehus/pandemic_toolbook)** — a pandemic
>   simulation playground with known ground truth, paired with an outbreak-analysis toolbox.
>   It began on the `claude/pandemic-simulation-playground-jlg1rs` branch here before moving to
>   its own repo.

The season-summarising engine reads a per-season sense of population **susceptibility** off **flu ILI+** waves with a
susceptible-reconstruction **SIR model**, through a small **swappable-method framework** plus a
reference Bayesian model:

- **Method framework (R)** — `code/01_main_supporting/sir_core.R` (shared SIR engine) +
  `code/01_main_supporting/methods/` (one file per method) + `methods_registry.R` (registry +
  common per-season summary schema). Every method fixes R0, the infectious period and the seed
  `I0` (from settings) and fits, per season, the **susceptibility** `S0` and a **reporting
  fraction** `c`, with a shared baseline `b` and overdispersion `phi`. Methods so far:
  - `deterministic` — the season's wave *is* a single deterministic SIR; observations are
    overdispersed noise around it. `S0` is identified by the wave's rise rate. Transparent; misfit
    shows up as honest residuals.
  - `ekf` — the same model with an Extended Kalman Filter, admitting **process noise** so the
    epidemic can wander off the deterministic SIR (one shared `q_I` per country, fitted but
    regularised small). The two methods agree on wave shape but the EKF draws more contrast between
    seasons' `S0` — comparing them is exactly what the framework is for.
  - `descriptive` — non-mechanistic: **smooths** each curve (centered moving average) and reads off
    features (AUC, peak height, onset week, steepness) directly. It fits no SIR and reports no `S0` —
    the observed-shape features are compared *within* a country, not mapped onto the SIR
    susceptibility axis (that mechanistic-vs-phenomenological mapping is deliberately out of scope).
- **Bayesian SIR (Stan)** — `stan/SIR_multiseason_age_vax_2.stan`: an age- and vaccination-
  structured, multi-season SIR fit by HMC, with scenario projections (reference model).

The repo also carries the full data layer: ECDC **ERVISS** + **RespiCompass** surveillance
(ILI/ARI, typing/positivity), contact matrices, vaccination and demography, turned into tidy
model-ready tables, with a data-quality / dynamics report.

## Quick start

```r
renv::restore()                       # install pinned dependencies (renv.lock)
source("code/00_main.R")              # build the data and model inputs
```

Fit per-season **susceptibility** across a country's seasons (base R; no pipeline needed — uses the
committed slim panel `data/slim_flu_iliplus.csv`):

```r
source("code/02_settings/settings_version0.R"); params <- settings()
source("code/01_main_supporting/sir_core.R")
source("code/01_main_supporting/methods/method_sir_deterministic.R")
source("code/01_main_supporting/methods_registry.R")

sl  <- load_flu_iliplus_slim("DK")                  # one country's per-season weekly ILI+
fit <- run_method("deterministic", sl, params)      # any registered method, same call
summarise_method_fit(fit)                           # tidy: S0, R_eff, c, peak/onset week, cor, ...
```

`fit$params$S0` is the per-season susceptibility — **interpret the within-country RANKING across
seasons**, not the absolute level (conditional on the fixed R0/seed) and not the fine spacing of the
gaps (method-dependent, esp. for the EKF — see `documentation/decisions.md`).

Run every registered method and save a figure per method:

```sh
Rscript code/04_modelling/fit_methods_demo.R        # -> output/fit_<method>.png
```

## The model in one paragraph

A latent SIR in population proportions (`dS/dt = -beta S I`, `dI/dt = beta S I - gamma I`,
`beta = R0·gamma`) drives weekly ILI+, taken proportional to weekly new infections
(`E[y] = c · new infections + b`). Observation noise is neg-binomial-like (`Var = mu + mu²/phi`).
With R0, the infectious period and the seed `I0` fixed, the wave's **rise rate**
`r = gamma·(R0·S0 − 1)` identifies the per-season susceptibility `S0`; `c` carries the reporting
scale and `b` the off-season baseline. The deterministic method trusts the SIR fully (residuals are
noise); the EKF method adds process noise so the latent state can deviate from it. A single-season,
single-strain SIR does not capture secondary-strain or strongly NPI-shaped seasons, by design.

## Tests

```r
Rscript run_tests.R                   # offline, from the committed data/ snapshots
```

Checks the data contracts and canonical tables; that each method converges with plausible
parameters and reproduces the observed waves (`tests/testthat/test-sir-deterministic.R`); that
every registered method honours the common summary schema (`tests/testthat/test-methods-registry.R`);
that the committed slim panel obeys its construction rules and is reproduced exactly by the build
(`tests/testthat/test-slim-panel.R`); and that both mechanistic methods RECOVER known S0 values (rank
order + calibrated tolerance) from data simulated by the generative model itself
(`tests/testthat/test-sir-recovery.R` — the S0-identifiability claim, exercised on known truth).

## Layout

```
code/00_main.R                 build data + model inputs (run-order map for everything in its header)
code/01_main_supporting/       setup, validate, load_data, gen_model_input, eyeballing,
                               stitch_iliplus (THE panel-assembly rules, one implementation),
                               sir_core (shared SIR engine + multi-start optimiser + loaders),
                               methods/ (one file per fitting method),
                               methods_registry (registry + per-season summary schema),
                               send_report (manual email utility)
code/02_settings/              settings_version0.R (live params only, incl. susc_* fixed values)
code/03_report/                eyeballing_report.Rmd (data-quality / dynamics report),
                               data_availability.R + driver_availability.R (coverage heatmaps)
code/04_modelling/             build_slim_panel.R (write the committed panel via the shared stitch),
                               fit_methods_demo.R (every method, all seasons), descriptive_overview.R,
                               ekf_overview.R
code/06_comp_model/            the age x vaccination compartmental model (branch flu_comp_model):
                               ASSUMPTIONS.md (every assumption, with provenance and evidence),
                               comp_model_settings.R, contact_matrix.R (exact-R0 scaling),
                               comp_model_core.R (R reference) + comp_model_core.cpp (identical
                               C++ engine), comp_model_data.R, comp_model_fit.R (two-stage EKF),
                               comp_model_report.R (eyeballing figures), run_comp_model.R,
                               comp_model_joint.R + run_comp_model_joint.R (shared season R0)
code/05_analysis/              the driver analysis: analysis_helpers.R (shared Gibbs samplers, rhat,
                               VE-vs-dominant rule), prepare_descriptors.R, dominant_subtype.R,
                               analyse_patterns.R, plot_patterns.R, plot_vax_scatter.R,
                               hierarchical_models.R, bayes_subtype.R, subtype_8season.R,
                               bayes_prior_burden.R, bayes_precovid_ve_subtype.R,
                               precovid_predict_postcovid.R
stan/                          SIR_multiseason_age_vax_2.stan (Bayesian SIR)
data/                          committed ERVISS / RespiCompass snapshots + slim_flu_iliplus.csv
                               (+ data/external/ driver tables with provenance)
output/                        cached data lists + figures (gitignored, regenerated)
tests/testthat/                contract + method + panel + parameter-recovery tests
documentation/                 quickstart, data_overview, decisions, analysis_strategy,
                               findings_descriptors, documentation.Rmd (see table below)
```

## Documentation

Where each kind of information lives:

| File | Holds |
|---|---|
| `README.md` | what the repo does, quick start, the method framework, layout |
| `documentation/quickstart.md` | how to set up and run |
| `documentation/data_overview.md` | what data is present (`data`, `models_in`, indicators) |
| `documentation/documentation.Rmd` | the model maths / science (SIR, inference, contact matrix) |
| `documentation/decisions.md` | **why** — rationale for key modelling / method / data decisions |
| `documentation/analysis_strategy.md` | the driver analysis — strategy, principles, what we've learned, ranked next steps |
| `documentation/findings_descriptors.md` | results of the descriptor / driver analyses |
| `documentation/external_drivers.md` | externally-sourced drivers (subtype, vaccination, climate, VE) — provenance, values, caveats |
| `documentation/reflections.md` | project-relevance thoughts / interpretations / pondering (also inline `[REFLECTION]` tags) |
| `PROJECT_SCOPE.md` | project scope (in / out of scope), research aim, references |
| inline header comments | why each function / file works the way it does |

New design decisions go in `documentation/decisions.md` (append-only), so the reasoning is kept
with the code rather than only in commit messages.

## Reproducibility

`renv.lock` pins all dependencies (R 4.3.3); `renv::restore()` reproduces the environment.
Data is public ECDC ERVISS / RespiCompass surveillance data. The build chain regenerates every
analysis input from the committed snapshots: `code/00_main.R` writes the model-ready inputs to
`output/models_in.rds` (a gitignored cache), from which `code/04_modelling/build_slim_panel.R`
reproduces the committed panel `data/slim_flu_iliplus.csv` (verified byte-identical) and
`code/05_analysis/prepare_descriptors.R` the descriptor tables the analyses consume. See
`documentation/quickstart.md` and `documentation/data_overview.md` for more.

## Note on the Stan model

The model's verified indexing, accumulator and population-conservation bugs are now fixed, and a
minimal proper prior set is re-enabled. Still pending before production HMC use: (a) the per-season
effect parameters (`i_season`, `r_season`, `prop_ili_season`) remain disabled, so the model fits a
season-invariant initial state and reporting ratio; (b) the remaining commented-out priors need
Jacobian-aware re-derivation; (c) it is not compile-tested in this environment (no Stan toolchain).
The R method framework remains the ready-to-run path.
