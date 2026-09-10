# run_age_experiment.R -- WHICH AGE MECHANISM? (ASSUMPTIONS.md C7 / F2 / open decision 1)
#
# Age-invariant reporting over-predicted the medium group and under-predicted the young and the
# elderly in the Danish pilot. Two mechanisms can absorb that and they are NOT distinguishable on one
# country: age-specific REPORTING (ILI+ per infection differs by age; the dynamics are untouched) or
# age-specific SUSCEPTIBILITY (the attack-rate profile itself differs; enters the force of infection).
# The cross-country logic (owner, 2026-09): a susceptibility effect is biology and should show the SAME
# direction of over/under-prediction by age in every country when reporting is age-invariant;
# reporting effects are properties of surveillance systems and should differ between countries.
# The external anchor is the PHIRST cohort (South Africa; reporting-independent infection incidence
# by age: young > adult > elderly; see phirst_attack_by_age() in comp_model_report.R).
#
# Four variants per country (the two switches of comp_model_settings()):
#   A_none            c_by_age = FALSE, susc_by_age = FALSE   (the contact structure alone sets the age profile)
#   B_reporting       c_by_age = TRUE,  susc_by_age = FALSE
#   C_susceptibility  c_by_age = FALSE, susc_by_age = TRUE
#   D_both            c_by_age = TRUE,  susc_by_age = TRUE
# B and C have the same number of parameters, so their likelihoods compare directly.
#
# Usage:  Rscript code/06_comp_model/run_age_experiment.R [part] [CC,CC,...] [variant,variant]
# (no arguments: every country and variant in one process, ~1 min per fit). Each invocation writes
# output/comp_model/age_experiment/age_experiment_<part>.csv plus one fit .rds per country x variant;
# collect_age_experiment() merges the parts, and plot_cm_age_experiment() draws the figure.

suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_cpp.R")
source("code/06_comp_model/comp_model_data.R"); source("code/06_comp_model/comp_model_fit.R")
source("code/06_comp_model/comp_model_report.R")
cm_load_cpp()

age_experiment_countries = c("DK", "EE", "ES", "FR", "NO", "BE", "CZ", "IE", "IT", "PL", "HR", "NL")   # >= 6 age-complete seasons (F5)
age_experiment_variants  = list(A_none = c(c_by_age = FALSE, susc_by_age = FALSE), B_reporting = c(c_by_age = TRUE, susc_by_age = FALSE),
                                C_susceptibility = c(c_by_age = FALSE, susc_by_age = TRUE), D_both = c(c_by_age = TRUE, susc_by_age = TRUE))
age_experiment_dir = "output/comp_model/age_experiment"

# one row per country x variant: likelihoods, the fitted age effects, modelled attack rates (mean over
# seasons) and the obs/model ratio of cumulative ILI+ by age (mean over seasons, deterministic curve)
age_experiment_row = function(fit, cc, v){
  ratio = colMeans(do.call(rbind, lapply(seq_along(fit$seasons), function(k){ y = fit$y[[k]]; m = fit$mu_det[[k]]; ok = is.finite(y)
    vapply(seq_len(ncol(y)), function(a) sum(y[ok[, a], a]) / sum(m[ok[, a], a]), numeric(1)) })), na.rm = TRUE)
  att = colMeans(fit$attack)
  # negll_* are the PENALISED objectives (each N(0,1) age prior adds 0.92 nats at its centre, so they
  # are not comparable across variants); loglik_ekf is the pure EKF log-likelihood used for the gains
  data.frame(country = cc, variant = v, n_seasons = length(fit$seasons), negll_stage1 = if (is.null(fit$stage1)) NA_real_ else fit$stage1$negll, negll_ekf = fit$negll,
             loglik_ekf = fit$loglik,
             phi = fit$params$phi, S0 = fit$params$S0,
             c_young_rel = fit$params$c_age[1] / fit$params$c_age[2], c_elderly_rel = fit$params$c_age[3] / fit$params$c_age[2],
             se_log2c_young = unname(fit$se["log2c_young"])[1], se_log2c_elderly = unname(fit$se["log2c_elderly"])[1],
             susc_young = fit$params$sigma[1], susc_elderly = fit$params$sigma[3],
             se_log2susc_young = unname(fit$se["log2susc_young"])[1], se_log2susc_elderly = unname(fit$se["log2susc_elderly"])[1],
             attack_young = att[1], attack_medium = att[2], attack_elderly = att[3],
             ratio_young = ratio[1], ratio_medium = ratio[2], ratio_elderly = ratio[3], stringsAsFactors = FALSE)
}

run_age_experiment = function(countries = age_experiment_countries, variants = names(age_experiment_variants), part = "all",
                              dir = age_experiment_dir, models_in = readRDS("output/models_in.rds"), demo = NULL){
  if (is.null(demo)){ load("output/demography_respicast.Rdata"); demo = obj }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE); rows = list()
  for (cc in countries){
    cd = tryCatch(build_comp_data(cc, models_in, demo, comp_model_settings()), error = function(e){ cat("SKIP", cc, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(cd)) next
    for (v in variants){
      s = comp_model_settings(); s$c_by_age = unname(age_experiment_variants[[v]]["c_by_age"]); s$susc_by_age = unname(age_experiment_variants[[v]]["susc_by_age"])
      # WARM STARTS: B and C nest A and start at A's deterministic optimum (new slots at 0) -- from the
      # data-driven start alone the susceptibility variants wandered into worse optima than A. D nests
      # both B and C and the objective is BIMODAL (a reporting-like and a susceptibility-like optimum),
      # so D starts from each of the two and keeps the better; a missing parent falls back to A.
      parents = if (v == "D_both") c("B_reporting", "C_susceptibility", "D_both") else "A_none"   # D also from its own earlier fit, if any
      warms = lapply(parents, function(pv){ fp = file.path(dir, sprintf("fit_%s_%s.rds", cc, pv)); if (file.exists(fp)) readRDS(fp)$stage1$theta else NULL })
      warms = Filter(Negate(is.null), warms); if (!length(warms)) warms = list(NULL)
      t0 = Sys.time(); fit = NULL
      for (w in warms){
        fw = tryCatch(fit_comp_model(cd, s, R0_free = TRUE, verbose = FALSE, start = w), error = function(e){ cat("FAIL", cc, v, ":", conditionMessage(e), "\n"); NULL })
        if (!is.null(fw) && (is.null(fit) || fw$negll < fit$negll)) fit = fw
      }
      # a nested variant can never be worse than its parent: the parent's EKF optimum padded with zeros
      # IS a valid point of the larger model. Polish from there (EKF stage only) and keep the best --
      # this closes the ~4-nat optimiser noise that otherwise shows D 'losing' to B.
      if (v != "A_none"){
        s_polish = s; s_polish$two_stage = FALSE
        for (pv in parents[parents != v]){ fp = file.path(dir, sprintf("fit_%s_%s.rds", cc, pv)); if (!file.exists(fp)) next
          fw = tryCatch(fit_comp_model(cd, s_polish, R0_free = TRUE, n_starts = 1, verbose = FALSE, start = readRDS(fp)$theta), error = function(e) NULL)
          if (!is.null(fw) && (is.null(fit) || fw$negll < fit$negll)){ fw$stage1 = fit$stage1; fit = fw }   # keep a stage-1 record for the figures
        }
      }
      if (is.null(fit)) next
      fit$Cn = cd$Cn; saveRDS(fit, file.path(dir, sprintf("fit_%s_%s.rds", cc, v)))
      rows[[length(rows) + 1]] = r = age_experiment_row(fit, cc, v)
      cat(sprintf("%s %-16s negll1 %7.1f  ekf %7.1f  c_rel young %.2f eld %.2f | susc young %.2f eld %.2f | attack %.2f/%.2f/%.2f | obs/model %.2f/%.2f/%.2f | %.1f min\n",
                  cc, v, r$negll_stage1, r$negll_ekf, r$c_young_rel, r$c_elderly_rel, r$susc_young, r$susc_elderly,
                  r$attack_young, r$attack_medium, r$attack_elderly, r$ratio_young, r$ratio_medium, r$ratio_elderly,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      write.csv(do.call(rbind, rows), file.path(dir, sprintf("age_experiment_%s.csv", part)), row.names = FALSE)
    }
  }
  invisible(do.call(rbind, rows))
}

# the merged table, rebuilt from the saved fits (one .rds per country x variant; the CSV parts are
# only progress logs). Fits whose objective never became finite (negll >= 1e9) are dropped.
collect_age_experiment = function(dir = age_experiment_dir){
  fs = list.files(dir, "^fit_[A-Z]{2}_[A-Z]_.*\\.rds$", full.names = TRUE)
  rows = lapply(fs, function(fn){ m = regmatches(basename(fn), regexec("^fit_([A-Z]{2})_(.*)\\.rds$", basename(fn)))[[1]]
    fit = readRDS(fn); if (!is.finite(fit$negll) || fit$negll >= 1e9) return(NULL); age_experiment_row(fit, m[2], m[3]) })
  tab = do.call(rbind, rows); tab[order(tab$country, tab$variant), ]
}

if (sys.nframe() == 0){
  args = commandArgs(TRUE)
  part = if (length(args) >= 1) args[1] else "all"
  countries = if (length(args) >= 2) strsplit(args[2], ",")[[1]] else age_experiment_countries
  variants = if (length(args) >= 3) strsplit(args[3], ",")[[1]] else names(age_experiment_variants)
  run_age_experiment(countries, variants, part)
  cat("DONE", part, "\n")
}
