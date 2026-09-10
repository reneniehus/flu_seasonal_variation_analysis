# comp_model_data.R -- assemble one country's inputs for the compartmental model
#
# Everything here follows ASSUMPTIONS.md sections 2, 4 (D2) and 6 (F1, F3, F5, F6):
#   - weekly ILI+ per age band from the panel stitch (stitch_iliplus_by_age: same source rules and
#     total-based alignment factor as the committed panel; only panel country-seasons);
#   - the four bands merged to young 0-14 | medium 15-64 | elderly 65+ by population-weighted rates;
#   - rates (per 100 000 of the group) converted to COUNTS with EACH group's own population;
#   - contact matrix collapsed to the three groups and scaled to spectral radius 1;
#   - vaccination: the season's 65+ coverage as one pulse on 1 October (season day 62), young and
#     medium 0; coverage from the observed history, else data/external post-COVID, else the
#     RespiCompass scenario midpoint, else the country's last known value (each provenance recorded).
# Requires setup.R (tidyverse, EU_long), stitch_iliplus.R, contact_matrix.R, comp_model_settings.R.

# ---- |-populations per model group from the fine pyramid (5-year bands) ----
.cm_populations = function(demo_fine, country_long){
  p = demo_fine[demo_fine$country == country_long, ]
  lo = suppressWarnings(as.integer(sub("[-+].*$", "", p$age_group)))     # "0-4" -> 0, "80+" -> 80
  c(young = sum(p$population[lo <= 10]), medium = sum(p$population[lo >= 15 & lo <= 60]),
    elderly = sum(p$population[lo >= 65]))
}

# ---- |-65+ coverage per season with provenance ----
.cm_vax_coverage = function(country, seasons, models_in, ext_path = "data/external/vaccination_coverage_65plus_postcovid.csv"){
  L = models_in$data_timeseries_long; cc = country                # cc: the external CSV has a 'country' COLUMN that would mask the argument
  hist = L %>% filter(stream == "vaccination_history_65plus", observed, country_short == cc) %>%
    distinct(season, value) %>% deframe()
  ext = if (file.exists(ext_path)) read.csv(ext_path, stringsAsFactors = FALSE) %>%
    filter(country_short == cc, !grepl("of invited", age_band)) %>% distinct(season, coverage_pct) %>%
    mutate(cov = coverage_pct / 100) %>% select(season, cov) %>% deframe() else c()
  scen = L %>% filter(stream == "vaccination_scenario", observed, country_short == cc,
                      scenario %in% c("higher_vax_coverage", "lower_vax_coverage")) %>% pull(value)
  scen_mid = if (length(scen)) mean(scen) else NA_real_
  out = data.frame(season = seasons, coverage = NA_real_, provenance = NA_character_, stringsAsFactors = FALSE)
  last = NA_real_
  for (i in seq_along(seasons)){
    s = seasons[i]
    if (!is.na(hist[s]))      { out$coverage[i] = hist[[s]];  out$provenance[i] = "observed history" }
    else if (!is.na(ext[s]))  { out$coverage[i] = ext[[s]];   out$provenance[i] = "data/external post-COVID" }
    else if (is.finite(scen_mid)) { out$coverage[i] = scen_mid; out$provenance[i] = "RespiCompass scenario midpoint" }
    else if (is.finite(last)) { out$coverage[i] = last;       out$provenance[i] = "carried forward" }
    else                      { out$coverage[i] = 0;          out$provenance[i] = "none available -> 0" }
    last = out$coverage[i]
  }
  out
}

# ---- |-the country's model data ----
build_comp_data = function(country, models_in, demo, settings,
                           panel = read.csv("data/slim_flu_iliplus.csv", stringsAsFactors = FALSE)){
  country_long = EU_long(country)
  N = .cm_populations(demo$population_pyramid_fine, country_long)
  N4 = demo$population_pyramid %>% filter(country == country_long) %>%
    { .$population[match(c("0-4","5-14","15-64","65+"), .$age_group)] }

  # age-specific weekly rates on the panel grid -> 3 groups -> counts
  by_age = stitch_iliplus_by_age(models_in, panel) %>% filter(country_short == country)
  seasons = sort(unique(by_age$season))
  stopifnot(length(seasons) > 0)
  band_pop = c(age_00_04 = N4[1], age_05_14 = N4[2], age_15_64 = N4[3], age_65_99 = N4[4])
  grp = c(age_00_04 = "young", age_05_14 = "young", age_15_64 = "medium", age_65_99 = "elderly")
  merged = by_age %>% mutate(group = grp[agegroup], w = band_pop[agegroup]) %>%
    group_by(season, week, group) %>%
    summarise(rate = if (all(is.na(value))) NA_real_ else sum(value * w, na.rm = TRUE) / sum(w[!is.na(value)]),
              .groups = "drop")                        # population-weighted mean of the bands present
  mk = function(s, col){
    d = merged %>% filter(season == s) %>% select(week, group, all_of(col)) %>%
      tidyr::pivot_wider(names_from = group, values_from = all_of(col)) %>% arrange(week)
    m = as.matrix(d[, c("young","medium","elderly")]); rownames(m) = d$week; m
  }
  rates = lapply(seasons, mk, col = "rate"); names(rates) = seasons
  y = lapply(rates, function(m) sweep(m, 2, N / settings$rate_per, "*")); # counts per group-week

  # contact matrix (4 -> 3 groups, spectral radius 1); EU average where the country has none
  ct = models_in$contacts
  C4 = if (!is.null(ct[[country_long]])) ct[[country_long]] else ct[["EU"]]
  cm = model_contact_matrix(C4, N4, settings)

  # vaccination pulse per season
  vax = .cm_vax_coverage(country, seasons, models_in)
  vax_day = as.integer(as.Date(paste0(substr(seasons[1], 1, 4), settings$vax_pulse_monthday)) -
                       as.Date(paste0(substr(seasons[1], 1, 4), settings$season_start_monthday))) + 1L
  vax_frac = lapply(seq_along(seasons), function(i) c(young = 0, medium = 0, elderly = vax$coverage[i]))

  list(country = country, country_long = country_long, seasons = seasons, groups = names(N),
       N = N, N4 = N4, Cn = cm$C, contact_source = if (!is.null(ct[[country_long]])) country_long else "EU average",
       rates = rates, y = y, n_weeks = vapply(y, nrow, integer(1)),
       vax = vax, vax_day = vax_day, vax_frac = vax_frac,
       source = panel %>% filter(country_short == country) %>% distinct(season, source))
}
