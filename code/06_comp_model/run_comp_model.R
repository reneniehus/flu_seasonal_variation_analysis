# run_comp_model.R -- per-country stage of the compartmental model (ASSUMPTIONS.md I1 stage ii)
#
# For each country: build the model data, fit S0_c + per-season R0_{c,s} (+ c, b, phi, q) with the
# EKF likelihood, save the fit, and write the eyeballing figures (comp_model_report.R). Countries
# default to those with age-complete ILI+ in every panel season (ASSUMPTIONS.md F5).
# Prerequisites: output/models_in.rds and output/demography_respicast.Rdata (code/00_main.R).
# Run from the repo root:  Rscript code/06_comp_model/run_comp_model.R [country ...]
#
# Runtime note (base-R reference): one likelihood evaluation is ~0.15 s for an 8-season country, a
# numerical gradient ~3.5 s, so a 4-start fit is ~20 min on 4 cores. The Rcpp port (stage iii) is
# for the joint fit, not for this stage.

suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R")
source("code/01_main_supporting/sir_core.R")                 # .fit_multistart
source("code/06_comp_model/contact_matrix.R")
source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R")
source("code/06_comp_model/comp_model_data.R")
source("code/06_comp_model/comp_model_fit.R")
source("code/06_comp_model/comp_model_report.R")

args = commandArgs(trailingOnly = TRUE)
countries = if (length(args)) args else c("DK", "FR", "EE", "ES", "NO")
settings  = comp_model_settings()
models_in = readRDS("output/models_in.rds")
load("output/demography_respicast.Rdata"); demo = obj
dir.create("output/comp_model", showWarnings = FALSE, recursive = TRUE)

fits = list()
for (cc in countries){
  cd  = build_comp_data(cc, models_in, demo, settings)
  cat(sprintf("\n== %s: %d seasons, groups N = %s, contacts = %s ==\n", cc, length(cd$seasons),
              paste(format(cd$N, big.mark = ","), collapse = " / "), cd$contact_source))
  elapsed = system.time(fit <- fit_comp_model(cd, settings, R0_free = TRUE))[["elapsed"]]
  cat(sprintf("   fitted in %.0f s\n", elapsed))
  fit$Cn = cd$Cn
  saveRDS(fit, sprintf("output/comp_model/fit_%s.rds", cc))
  fits[[cc]] = fit
}
save_cm_report(fits)
cat("\nfigures + summary -> output/comp_model/\n")
print(do.call(rbind, lapply(fits, summarise_comp_fit))[, c("country","season","S0","R0","R_eff","c","q","cor","attack_young","attack_medium","attack_elderly")],
      row.names = FALSE, digits = 3)
