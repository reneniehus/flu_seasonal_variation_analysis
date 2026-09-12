# run_joint_model.R -- fit the joint model to every country with at least 6 age-complete seasons,
# then write the diagnostics and the figures. See MODEL.md for what the model assumes.
#   Rscript code/07_joint_model/run_joint_model.R
setwd(here::here())
suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_data.R")
source("code/07_joint_model/joint_model.R"); source("code/07_joint_model/joint_report.R")

t_compile = system.time(jm_load_cpp())[["elapsed"]]
models_in = readRDS("output/models_in.rds"); load("output/demography_respicast.Rdata"); demo = obj

# every country the panel can age-structure; jm_build_data drops those with fewer than 6 seasons
candidates = c("DK", "EE", "ES", "FR", "NO", "BE", "CZ", "IE", "IT", "PL", "HR", "NL")
d = jm_build_data(candidates, models_in, demo, min_seasons = 6L)

cat("\n--- fitting ---\n")
fit = jm_fit(d, cores = max(1L, parallel::detectCores() - 1L))

cat("\n--- identifiability ---\n")
id = jm_identifiability(fit)
cat(sprintf("Hessians in %.0f s. Penalised Hessian positive definite: %s\n", id$seconds, id$post_pd))
cat(sprintf("likelihood-only spectrum: largest %.3g, smallest |eigenvalue| %.3g, ratio %.2g\n",
            max(abs(id$eig_lik)), min(abs(id$eig_lik)), id$flat_ratio))
cat(sprintf("directions the data cannot see (|eig| < 1e-8 x largest): %d ; negative directions: %d\n",
            id$n_near_flat, id$n_negative))
print(id$family, row.names = FALSE)

cat("\n--- model adequacy ---\n")
ad = jm_adequacy(fit); print(ad, row.names = FALSE)
cat(sprintf("median noise the fit needed / noise the data have: %.2fx\n", median(ad$excess, na.rm = TRUE)))

cat("\n--- parameters ---\n")
print(jm_summary_season(fit), row.names = FALSE)
print(jm_summary_country(fit), row.names = FALSE)

dir.create("output/joint_model", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(fit = fit, id = id, adequacy = ad, compile_seconds = t_compile),
        "output/joint_model/joint_fit.rds")
cat("\n--- figures ---\n")
save_jm_report(fit, id)
f = jm_fitted_cpp(fit$theta, d)
cors = vapply(seq_len(d$n_cs), function(i){ y = as.numeric(d$y[[i]]); m = as.numeric(f$mu[[i]])
  ok = is.finite(y) & is.finite(m); if (sum(ok) > 3) cor(y[ok], m[ok]) else NA_real_ }, numeric(1))
cat(sprintf("obs vs model correlation per country-season: median %.3f, range %.3f - %.3f\n",
            median(cors, na.rm = TRUE), min(cors, na.rm = TRUE), max(cors, na.rm = TRUE)))
cat(sprintf("\nTOTAL: compile %.1f s, fit %.1f s, Hessians %.0f s\n", t_compile, fit$seconds, id$seconds))
