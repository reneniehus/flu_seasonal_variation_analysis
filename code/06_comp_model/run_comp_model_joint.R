# run_comp_model_joint.R -- stage iv: all countries and seasons at once, R0_s shared across countries
#
# Starts from the per-country fits in output/comp_model/fit_<CC>.rds (run_comp_model.R first), fits
# the joint model with the C++ engine (comp_model_joint.R), and writes:
#   output/comp_model/joint_fit.rds, joint_seasons.csv, joint_countries.csv
#   output/comp_model/joint_R0_seasons.png   shared R0_s with 95% intervals over the per-country lines
#   output/comp_model/joint_S0_countries.png  country S0_c with intervals
#   output/comp_model/joint_fit_<CC>.png      age x season fit grids under the JOINT parameters
# Run from the repo root:  Rscript code/06_comp_model/run_comp_model_joint.R [country ...]

suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R")
source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R")
source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R")
source("code/06_comp_model/comp_model_cpp.R")
source("code/06_comp_model/comp_model_data.R")
source("code/06_comp_model/comp_model_fit.R")
source("code/06_comp_model/comp_model_joint.R")
source("code/06_comp_model/comp_model_report.R")
cm_load_cpp()

args = commandArgs(trailingOnly = TRUE)
settings  = comp_model_settings()
fit_files = list.files("output/comp_model", "^fit_[A-Z]{2}\\.rds$", full.names = TRUE)
run_countries = if (length(args)) args else sub("^fit_([A-Z]{2})\\.rds$", "\\1", basename(fit_files))
stopifnot(length(run_countries) >= 2)
models_in = readRDS("output/models_in.rds"); load("output/demography_respicast.Rdata"); demo = obj

fits = setNames(lapply(run_countries, function(cc) readRDS(sprintf("output/comp_model/fit_%s.rds", cc))), run_countries)
cds  = setNames(lapply(run_countries, function(cc) build_comp_data(cc, models_in, demo, settings)), run_countries)

jf = fit_comp_model_joint(cds, fits, settings, engine = "cpp")
saveRDS(jf, "output/comp_model/joint_fit.rds")
sj = summarise_joint_fit(jf)
write.csv(sj$seasons, "output/comp_model/joint_seasons.csv", row.names = FALSE)
write.csv(sj$countries, "output/comp_model/joint_countries.csv", row.names = FALSE)
cat("\n== shared season transmissibility R0_s (joint fit, 95% Laplace intervals) ==\n"); print(sj$seasons, digits = 3, row.names = FALSE)
cat("\n== country susceptibility S0_c ==\n"); print(sj$countries, digits = 3, row.names = FALSE)

# ---- figures ----
per_country = do.call(rbind, lapply(fits, summarise_comp_fit))
p1 = ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5 * exp(-1.96 * settings$prior_logR0_sd), ymax = 1.5 * exp(1.96 * settings$prior_logR0_sd), fill = "#8da0cb", alpha = 0.25) +
  geom_hline(yintercept = 1.5, colour = "grey60", linetype = "dashed") +
  geom_line(data = per_country, aes(season, R0, group = country, colour = country), alpha = 0.5) +
  geom_point(data = per_country, aes(season, R0, colour = country), alpha = 0.5, size = 1.6) +
  geom_pointrange(data = sj$seasons, aes(season, R0, ymin = lo, ymax = hi), colour = "black", size = 0.6, linewidth = 0.9) +
  labs(title = "Season transmissibility: joint fit (black, shared across countries, 95% interval) over the per-country fits (colours)",
       subtitle = "blue band = the prior 95% around 1.5. Where the black points leave the band, the data pull the season away from 'typical' despite the strong prior.",
       x = NULL, y = "R0_s") + theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.subtitle = element_text(size = 8))
ggsave("output/comp_model/joint_R0_seasons.png", p1, width = 9.5, height = 5, dpi = 110)
p2 = ggplot(sj$countries, aes(S0, reorder(country, S0))) + geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#d95f02") +
  labs(title = "Country susceptibility S0_c under the joint fit (shared across seasons)", subtitle = "Laplace 95% intervals; relative comparison only", x = "S0", y = NULL) +
  theme_minimal(base_size = 10)
ggsave("output/comp_model/joint_S0_countries.png", p2, width = 6, height = 4, dpi = 110)

# per-country fit grids under the joint parameters: build fit-like objects from the joint filtered output
for (i in seq_along(run_countries)){
  cc = run_countries[i]; cd = cds[[cc]]; f = cm_fixed(cd$Cn, cd$N, settings); K = length(cd$seasons)
  R0s = jf$params$R0[match(cd$seasons, jf$seasons)]
  det = lapply(seq_len(K), function(s){ sim = cm_simulate_season(f, jf$params$S0[i], R0s[s], cd$n_weeks[s], cd$vax_day, cd$vax_frac[[s]], I0 = jf$params$I0[[i]][s])
    list(mu = cm_mu(sim$inc, f, jf$params$c[i], jf$params$b[i]), attack = sim$attack) })   # joint stage: scalar c per country (age/season deviations: TODO)
  fake = list(country = paste0(cc, " (joint)"), seasons = cd$seasons, groups = cd$groups, N = cd$N, settings = settings,
              params = list(S0 = jf$params$S0[i], R0 = R0s, I0 = jf$params$I0[[i]], c = jf$params$c[i], b = jf$params$b[i], phi = jf$params$phi[i], q = jf$params$q),
              y = cd$y, mu_filt = lapply(jf$filt[[i]], `[[`, "mu_pred"), mu_det = lapply(det, `[[`, "mu"),
              attack = do.call(rbind, lapply(det, `[[`, "attack")), R_eff = R0s * jf$params$S0[i], vax = cd$vax)
  ggsave(sprintf("output/comp_model/joint_fit_%s.png", cc), plot_cm_fit(fake), width = 2.6 * K + 2, height = 7.5, dpi = 110)
}
cat("\nfigures -> output/comp_model/joint_*.png\n")
