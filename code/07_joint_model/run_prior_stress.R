# run_prior_stress.R -- is "seasons differ far more in visibility than in how easily they spread" the
# data's or the priors'? Refits the working model under wider and hostile season priors and reports the
# spread of each season effect: susceptibility on the log scale of R0 x S0 for a typical country, and
# visibility on the log scale (MODEL.md, section 6). About 15 minutes on one core.
setwd(here::here())
suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_data.R")
source("code/07_joint_model/joint_model.R"); jm_load_cpp()
models_in = readRDS("output/models_in.rds"); load("output/demography_respicast.Rdata"); demo = obj
o = readRDS("output/joint_model/joint_fit.rds"); cands = o$fit$d$countries
spread = function(fit){ p = jm_unpack(fit$theta, fit$d); med = median(qlogis(vapply(p$country, `[[`, 0, "S0")))
  c(Reff_sd_log = sd(log(fit$d$R0_fixed * plogis(med + p$x))), vis_sd_log = sd(p$delta)) }
rows = list(as_fitted = c(x_sd = 0.5, delta_sd = 0.5, spread(o$fit), loglik = o$fit$loglik))
for (cfg in list(c(x = 2.0, del = 0.5), c(x = 0.5, del = 0.125), c(x = 2.0, del = 2.0))){
  s = modifyList(jm_settings(), list(pr_x_sd = cfg[["x"]], pr_delta_sd = cfg[["del"]]))
  d = jm_build_data(cands, models_in, demo, set = s, verbose = FALSE)
  f = jm_fit(d, theta0 = o$fit$theta, cores = max(1L, parallel::detectCores() - 1L), verbose = FALSE)
  rows[[sprintf("season prior sd: S0 %.3g, visibility %.3g", cfg[["x"]], cfg[["del"]])]] =
    c(x_sd = cfg[["x"]], delta_sd = cfg[["del"]], spread(f), loglik = f$loglik)
}
tab = do.call(rbind, rows); print(tab, digits = 4)
saveRDS(tab, "output/joint_model/prior_stress_R0_2.0.rds")
