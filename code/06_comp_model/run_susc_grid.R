# run_susc_grid.R -- IS THERE ONE AGE-SUSCEPTIBILITY PROFILE FOR ALL COUNTRIES? (ASSUMPTIONS.md C7)
#
# Within one country, age-specific susceptibility (biology) and age-specific reporting (surveillance)
# trade off almost perfectly: run_age_experiment.R's variant D moves the elderly effect freely between
# the two. The two-way logic separates them ACROSS countries: biology should be the same everywhere,
# surveillance is not. This script computes the profile likelihood of a susceptibility profile
# (young, elderly; medium = 1) SHARED by every country, with the reporting offsets FREE per country:
# for each grid point the profile is fixed (settings$susc_fixed), every country is refitted (warm-
# started from its reporting-only optimum), and the negative log-likelihoods are summed. A profile
# that lowers the total by many nats over (1, 1) is what the weekly data, freed of each country's
# reporting level, say about biology -- to be read against the PHIRST attack-rate ordering.
#
# Usage: Rscript code/06_comp_model/run_susc_grid.R [part] [CC,CC,...]  (writes susc_grid_<part>.csv)

suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_cpp.R")
source("code/06_comp_model/comp_model_data.R"); source("code/06_comp_model/comp_model_fit.R")
source("code/06_comp_model/run_age_experiment.R")            # countries, dir, age_experiment_row
cm_load_cpp()

susc_grid = expand.grid(young = c(1, 1.4, 1.8, 2.4), elderly = c(1, 1.7, 2.8))
susc_grid_dir = "output/comp_model/susc_grid"

run_susc_grid = function(countries = age_experiment_countries, grid = susc_grid, part = "all", dir = susc_grid_dir,
                         models_in = readRDS("output/models_in.rds"), demo = NULL){
  if (is.null(demo)){ load("output/demography_respicast.Rdata"); demo = obj }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE); rows = list()
  for (cc in countries){
    cd = tryCatch(build_comp_data(cc, models_in, demo, comp_model_settings()), error = function(e){ cat("SKIP", cc, ":", conditionMessage(e), "\n"); NULL })
    if (is.null(cd)) next
    fb = file.path(age_experiment_dir, sprintf("fit_%s_B_reporting.rds", cc))
    warm = if (file.exists(fb)) readRDS(fb)$stage1$theta else NULL
    for (g in seq_len(nrow(grid))){
      s = comp_model_settings(); s$c_by_age = TRUE; s$susc_by_age = FALSE; s$susc_fixed = c(grid$young[g], 1, grid$elderly[g])
      t0 = Sys.time()
      fit = tryCatch(fit_comp_model(cd, s, R0_free = TRUE, n_starts = 2, verbose = FALSE, start = warm), error = function(e){ cat("FAIL", cc, g, ":", conditionMessage(e), "\n"); NULL })
      if (is.null(fit)) next
      r = age_experiment_row(fit, cc, sprintf("y%.1f_e%.1f", grid$young[g], grid$elderly[g])); r$susc_young = grid$young[g]; r$susc_elderly = grid$elderly[g]
      rows[[length(rows) + 1]] = r
      cat(sprintf("%s young %.1f elderly %.1f  negll1 %7.1f ekf %7.1f  c_rel young %.2f eld %.2f  attack %.2f/%.2f/%.2f  %.1f min\n", cc, grid$young[g], grid$elderly[g],
                  r$negll_stage1, r$negll_ekf, r$c_young_rel, r$c_elderly_rel, r$attack_young, r$attack_medium, r$attack_elderly, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      write.csv(do.call(rbind, rows), file.path(dir, sprintf("susc_grid_%s.csv", part)), row.names = FALSE)
    }
  }
  invisible(do.call(rbind, rows))
}

# merged table plus the profile: total negll over countries per grid point, relative to (1, 1)
collect_susc_grid = function(dir = susc_grid_dir){
  fs = list.files(dir, "^susc_grid_.*\\.csv$", full.names = TRUE)
  tab = do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
  tab = tab[!duplicated(tab[, c("country", "variant")], fromLast = TRUE), ]
  prof = tab %>% group_by(susc_young, susc_elderly) %>% summarise(n = n(), negll_ekf = sum(negll_ekf), negll_stage1 = sum(negll_stage1), .groups = "drop") %>%
    mutate(gain_ekf = negll_ekf[susc_young == 1 & susc_elderly == 1] - negll_ekf, gain_stage1 = negll_stage1[susc_young == 1 & susc_elderly == 1] - negll_stage1)
  list(tab = tab, profile = prof)
}

if (sys.nframe() == 0){
  args = commandArgs(TRUE)
  part = if (length(args) >= 1) args[1] else "all"
  countries = if (length(args) >= 2) strsplit(args[2], ",")[[1]] else age_experiment_countries
  run_susc_grid(countries, part = part)
  cat("DONE", part, "\n")
}
