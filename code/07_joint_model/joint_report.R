# joint_report.R -- figures for the joint model. Each plot carries its own interpretation in the
# subtitle: what the parameter means, and what the pattern in front of you is saying. The figures are
# meant to be readable without the code open.

suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

.jm_theme = function(base = 10){
  theme_minimal(base_size = base) +
    theme(plot.title = element_text(face = "bold", size = base + 3),
          plot.subtitle = element_text(size = base - 0.5, colour = "grey30", lineheight = 1.25),
          plot.caption = element_text(size = base - 1.5, colour = "grey45", hjust = 0),
          strip.text = element_text(face = "bold", size = base - 1),
          panel.grid.minor = element_blank(),
          legend.position = "top")
}
.jm_grp = c("young", "medium", "elderly")
.jm_gcol = c(young = "#1B7F5C", medium = "#2B5D8A", elderly = "#C1541E")

# ---- |-1. fit vs data: the compact overview, all ages pooled ----
plot_jm_fit_overview = function(fit){
  d = fit$d; f = jm_fitted_cpp(fit$theta, d)
  rows = lapply(seq_len(d$n_cs), function(i){
    ic = d$cs_country[i] + 1L; N = d$N[[ic]]; w = N / sum(N)
    obs = rowSums(sweep(d$y[[i]], 2, d$rate_per / N, "*") * matrix(w, nrow(d$y[[i]]), 3, byrow = TRUE), na.rm = FALSE)
    mod = rowSums(sweep(f$mu[[i]], 2, d$rate_per / N, "*") * matrix(w, nrow(d$y[[i]]), 3, byrow = TRUE))
    data.frame(country = d$countries[ic], season = d$seasons[d$cs_season[i] + 1L],
               week = seq_len(nrow(d$y[[i]])), observed = obs, model = mod)
  })
  dd = do.call(rbind, rows)
  ggplot(dd, aes(week)) +
    geom_point(aes(y = observed), size = 0.5, colour = "grey25", na.rm = TRUE) +
    geom_line(aes(y = model), colour = "#C1541E", linewidth = 0.65) +
    facet_grid(country ~ season, scales = "free_y") +
    labs(title = "Every country-season: the model against the data",
         subtitle = paste("Points are observed influenza-positive ILI consultations per 100 000, pooled over age groups by population weight.",
                          "The line is the model at the joint optimum. One shared transmissibility per season and one shared reporting",
                          "deviation per season have to work for all countries at once, so a systematically low or high line is informative,",
                          "not a failure of that country's own fit.", sep = "\n"),
         x = "week of the season (from 1 August)", y = "ILI+ per 100 000") +
    .jm_theme(9) + theme(axis.text = element_text(size = 6.5), panel.spacing = unit(0.25, "lines"))
}

# ---- |-2. fit vs data by age, one country ----
plot_jm_fit_country = function(fit, cc){
  dd = jm_tidy_fit(fit) %>% filter(country == cc) %>% mutate(group = factor(group, levels = .jm_grp))
  ggplot(dd, aes(week, value, colour = group)) +
    geom_point(data = dd %>% filter(what == "observed"), size = 0.6, alpha = 0.8, na.rm = TRUE) +
    geom_line(data = dd %>% filter(what == "model"), linewidth = 0.7) +
    facet_grid(group ~ season, scales = "free_y") +
    scale_colour_manual(values = .jm_gcol, guide = "none") +
    labs(title = paste0(cc, ": the three age groups, season by season"),
         subtitle = paste("Points observed, lines the model. The age profile is set by the contact matrix, one global elderly susceptibility",
                          "and this country's two age reporting offsets; nothing here varies by season, so all season-to-season",
                          "movement comes from the shared transmissibility, the shared reporting deviation and this wave's seed.", sep = "\n"),
         x = "week of the season", y = "ILI+ per 100 000 of the age group") +
    .jm_theme(9) + theme(axis.text = element_text(size = 6.5))
}

# ---- |-3. the two shared season parameters ----
plot_jm_seasons = function(fit){
  s = jm_summary_season(fit); d = fit$d
  s$season = factor(s$season, levels = d$seasons)
  p1 = ggplot(s, aes(season, R0, group = 1)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5 * exp(-1.96 * d$pr_R0_sd), ymax = 1.5 * exp(1.96 * d$pr_R0_sd),
             fill = "#2B5D8A", alpha = 0.10) +
    geom_hline(yintercept = 1.5, linetype = "dashed", colour = "grey50") +
    geom_line(colour = "#2B5D8A", linewidth = 0.8) + geom_point(size = 3, colour = "#2B5D8A") +
    geom_text(aes(label = sprintf("%.2f", R0)), vjust = -1.1, size = 3, colour = "grey20") +
    labs(title = "Transmissibility of each season's virus, shared by all countries",
         subtitle = paste("R0 is how many people one case would infect in a fully susceptible population, given the age mixing.",
                          "It is ONE number per season for the whole of Europe, so it is carried by every wave in that season.",
                          "Shaded band is the prior's 95% range around 1.5; the dashed line is its centre. Movement away from the",
                          "centre is the data speaking, since the prior pulls towards 1.5.", sep = "\n"),
         x = NULL, y = expression(R[0]~"of the season")) +
    .jm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  p2 = ggplot(s, aes(season, reporting_mult, group = 1)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(colour = "#C1541E", linewidth = 0.8) + geom_point(size = 3, colour = "#C1541E") +
    geom_text(aes(label = sprintf("%.2f", reporting_mult)), vjust = -1.1, size = 3, colour = "grey20") +
    scale_y_log10() +
    labs(title = "How visible each season was, per infection, shared by all countries",
         subtitle = paste("The season observation deviation: positive consultations per infection that season, relative to normal.",
                          "1.3 means a third more consultations for the same number of infections. It mixes how symptomatic the",
                          "season's strain was with any season size the mechanism cannot produce. Constrained to average one,",
                          "which is what makes it separable from each country's own reporting level.", sep = "\n"),
         x = NULL, y = "multiplier on reporting (log scale)") +
    .jm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  list(R0 = p1, deviation = p2)
}

# ---- |-4. the country parameters ----
plot_jm_countries = function(fit){
  s = jm_summary_country(fit); d = fit$d
  p1 = ggplot(s, aes(S0, reorder(country, S0))) +
    annotate("rect", xmin = plogis(d$pr_S0_mean - 1.96 * d$pr_S0_sd), xmax = plogis(d$pr_S0_mean + 1.96 * d$pr_S0_sd),
             ymin = -Inf, ymax = Inf, fill = "#2B5D8A", alpha = 0.08) +
    geom_vline(xintercept = plogis(d$pr_S0_mean), linetype = "dashed", colour = "grey50") +
    geom_segment(aes(x = 0, xend = S0, yend = reorder(country, S0)), colour = "grey80") +
    geom_point(size = 3.2, colour = "#2B5D8A") +
    geom_text(aes(label = sprintf("%.3f", S0)), hjust = -0.35, size = 3, colour = "grey20") +
    scale_x_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
    labs(title = "Susceptible fraction at the start of a season, by country",
         subtitle = paste("S0 is the share of the population that could be infected on 1 August, held the same across that country's",
                          "seasons. It is the project's target quantity. It is identified here only because transmissibility is shared",
                          "across countries: fitted country by country, the two trade off almost perfectly. Shaded band is the prior's",
                          "95% range. Read the RANKING with more confidence than the absolute level.", sep = "\n"),
         x = "S0", y = NULL) + .jm_theme()
  age = s %>% select(country, young = rel_young, elderly = rel_elderly) %>%
    pivot_longer(c(young, elderly), names_to = "group", values_to = "rel") %>%
    mutate(group = factor(group, levels = c("young", "elderly")))
  p2 = ggplot(age, aes(rel, reorder(country, rel), colour = group)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3) + scale_x_log10() +
    scale_colour_manual(values = .jm_gcol[c("young", "elderly")], name = NULL) +
    labs(title = "How much more visible a child or an elderly person is, per infection",
         subtitle = paste("Age reporting offsets, relative to adults, one pair per country. Above one means that group generates more",
                          "positive consultations per infection than an adult does: they consult more readily, or the sentinel network",
                          "sees them more. This is surveillance, not biology, which is exactly why it varies by country while the",
                          "elderly susceptibility below is one number for everyone.", sep = "\n"),
         x = "reporting relative to adults (log scale)", y = NULL) + .jm_theme()
  p3 = ggplot(s, aes(c_adult * 100, reorder(country, c_adult))) +
    geom_segment(aes(x = 0, xend = c_adult * 100, yend = reorder(country, c_adult)), colour = "grey80") +
    geom_point(size = 3.2, colour = "#1B7F5C") +
    geom_text(aes(label = sprintf("%.1f%%", c_adult * 100)), hjust = -0.3, size = 3, colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Reporting: adult infections that become a counted positive consultation",
         subtitle = paste("The chance that one adult infection turns into an influenza-positive ILI consultation in that country's",
                          "surveillance. A country with a low value is not having fewer infections, it is seeing fewer of them, so",
                          "raw curves are not comparable across countries but these numbers make them so.", sep = "\n"),
         x = "per cent of adult infections reported", y = NULL) + .jm_theme()
  list(S0 = p1, age_offsets = p2, reporting = p3)
}

# ---- |-5. arrival time: what the seeds are really doing ----
plot_jm_arrival = function(fit){
  d = fit$d; p = jm_unpack(fit$theta, d)
  rows = do.call(rbind, lapply(seq_len(d$n_country), function(ic){
    ics = d$cs_of_country[[ic]] + 1L
    data.frame(country = d$countries[ic], season = d$seasons[d$cs_season[ics] + 1L],
               log10_I0 = log10(unname(p$country[[ic]]$I0)))
  }))
  rows$season = factor(rows$season, levels = d$seasons)
  ggplot(rows, aes(season, country, fill = log10_I0)) + geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.1f", log10_I0)), size = 2.7, colour = "grey15") +
    scale_fill_distiller(palette = "PuOr", direction = 1, name = expression(log[10]~I[0])) +
    labs(title = "Seed size, which is really the arrival time of each wave",
         subtitle = paste("The infected fraction the model plants on 1 August. It has no epidemiological meaning of its own: a seed ten",
                          "times smaller simply makes that wave arrive about two weeks later. It is the only parameter that can move a",
                          "wave sideways without changing its shape, which is why it is free for every country-season. Darker orange",
                          "= smaller seed = later wave.", sep = "\n"),
         x = NULL, y = NULL) +
    .jm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ---- |-6. modelled attack rates by age ----
plot_jm_attack = function(fit){
  d = fit$d; f = jm_fitted_cpp(fit$theta, d)
  att = as.data.frame(f$attack); names(att) = d$groups
  att$country = d$countries[d$cs_country + 1L]; att$season = d$seasons[d$cs_season + 1L]
  long = att %>% pivot_longer(all_of(d$groups), names_to = "group", values_to = "attack") %>%
    mutate(group = factor(group, levels = .jm_grp), season = factor(season, levels = d$seasons))
  ggplot(long, aes(season, attack, colour = group, group = interaction(country, group))) +
    geom_line(alpha = 0.35, linewidth = 0.5) +
    stat_summary(aes(group = group), fun = median, geom = "line", linewidth = 1.4) +
    scale_colour_manual(values = .jm_gcol, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "Who actually gets infected: modelled attack rate by age group",
         subtitle = paste("Share of each age group infected over the season. Thin lines are countries, thick lines the median across them.",
                          "This is the one output that is free of reporting: it comes from the contact matrix, the shared elderly",
                          "susceptibility and vaccination, not from how many consultations were counted. Compare with prospective",
                          "cohort evidence rather than with surveillance curves.", sep = "\n"),
         x = NULL, y = "attack rate") +
    .jm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ---- |-7. model adequacy: is the dispersion explainable as measurement noise? ----
plot_jm_adequacy = function(fit){
  a = jm_adequacy(fit)
  long = a %>% select(country, model = cv_fitted, data = cv_data) %>%
    pivot_longer(c(model, data), names_to = "what", values_to = "cv")
  ggplot(long, aes(cv, reorder(country, cv), colour = what)) +
    geom_line(aes(group = country), colour = "grey80", linewidth = 1) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = c(model = "#C1541E", data = "#2B5D8A"),
                        labels = c(data = "the data's own week-to-week scatter", model = "what the fit needed"), name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "How much noise the fit needed, against how much the data actually have",
         subtitle = paste("Blue is the observed scatter of the weekly series around its own 3-week average: a floor that no mean could",
                          "beat, so it is pure measurement noise. Orange is the noise the fitted dispersion implies. The gap is the",
                          "deterministic curve failing to follow the wave, absorbed as if it were measurement error. A large gap is the",
                          "case for putting the Kalman filter back, and it is the main thing to watch in this version.", sep = "\n"),
         x = "coefficient of variation of a weekly count", y = NULL) + .jm_theme()
}

# ---- |-8. identifiability: what the data pin and what the priors hold up ----
plot_jm_identifiability = function(id){
  f = id$family %>% mutate(family = reorder(family, contraction_med))
  ggplot(f, aes(contraction_med, family)) +
    geom_segment(aes(x = 0, xend = contraction_med, yend = family), colour = "grey80", linewidth = 1) +
    geom_point(aes(size = n), colour = "#2B5D8A") +
    geom_text(aes(label = sprintf("%.2f", contraction_med)), hjust = -0.5, size = 3, colour = "grey20") +
    scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
    scale_size_continuous(range = c(2, 5.5), name = "parameters") +
    labs(title = "For each kind of parameter: did the data decide, or the prior?",
         subtitle = paste("Prior-to-posterior contraction. Zero means the fit gave back the prior unchanged, so that number carries no",
                          "information from the surveillance data. One means the data determined it and the prior was irrelevant.",
                          "Anything below about 0.3 should be reported as an assumption rather than a result.", sep = "\n"),
         x = "contraction (0 = all prior, 1 = all data)", y = NULL) + .jm_theme()
}

# ---- |-write everything ----
save_jm_report = function(fit, id = NULL, dir = "output/joint_model"){
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  d = fit$d
  ggsave(file.path(dir, "fit_overview.png"), plot_jm_fit_overview(fit),
         width = 2.05 * d$n_season + 2, height = 1.35 * d$n_country + 2.6, dpi = 115, limitsize = FALSE)
  for (cc in d$countries)
    ggsave(file.path(dir, sprintf("fit_%s.png", cc)), plot_jm_fit_country(fit, cc),
           width = 2.05 * d$n_season + 2, height = 7.2, dpi = 110, limitsize = FALSE)
  ps = plot_jm_seasons(fit)
  ggsave(file.path(dir, "season_R0.png"), ps$R0, width = 9.5, height = 5.4, dpi = 115)
  ggsave(file.path(dir, "season_deviation.png"), ps$deviation, width = 9.5, height = 5.4, dpi = 115)
  pc = plot_jm_countries(fit)
  ggsave(file.path(dir, "country_S0.png"), pc$S0, width = 9.5, height = 5.6, dpi = 115)
  ggsave(file.path(dir, "country_age_offsets.png"), pc$age_offsets, width = 9.5, height = 5.4, dpi = 115)
  ggsave(file.path(dir, "country_reporting.png"), pc$reporting, width = 9.5, height = 5.2, dpi = 115)
  ggsave(file.path(dir, "arrival_time.png"), plot_jm_arrival(fit), width = 9.5, height = 5.4, dpi = 115)
  ggsave(file.path(dir, "attack_rates.png"), plot_jm_attack(fit), width = 9.5, height = 5.6, dpi = 115)
  ggsave(file.path(dir, "adequacy.png"), plot_jm_adequacy(fit), width = 9.5, height = 5.4, dpi = 115)
  if (!is.null(id)) ggsave(file.path(dir, "identifiability.png"), plot_jm_identifiability(id),
                           width = 9.5, height = 5.2, dpi = 115)
  write.csv(jm_summary_season(fit), file.path(dir, "summary_season.csv"), row.names = FALSE)
  write.csv(jm_summary_country(fit), file.path(dir, "summary_country.csv"), row.names = FALSE)
  write.csv(jm_adequacy(fit), file.path(dir, "adequacy.csv"), row.names = FALSE)
  if (!is.null(id)) write.csv(id$table, file.path(dir, "identifiability.csv"), row.names = FALSE)
  invisible(dir)
}
