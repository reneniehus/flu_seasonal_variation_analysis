# joint_report.R -- THE DEFAULT FIGURE SET for the joint model.
#
# The set is numbered and answers three questions in order, so a reader can walk it front to back:
#   HOW THE MODEL WORKS       01 what varies where, 02 how a wave's arrival is set
#   WHAT DATA IT FITS         03 every country-season, 04 one country by age, 05 the noise budget
#   WHAT IT LEARNS            06 season transmissibility, 07 season visibility, 08 country
#                             susceptibility, 09 reporting level, 10 age reporting, 11 attack rates
#   WHETHER TO BELIEVE IT     12 data or prior, 13 recovery of a known truth
# Every figure carries its own interpretation in the subtitle: what the parameter means, and what the
# pattern in front of you is saying. They are meant to be readable without the code open.
# save_jm_report() writes the whole set plus a MANIFEST.md listing them in this order.

suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

.jm_theme = function(base = 10){
  theme_minimal(base_size = base) +
    theme(plot.title = element_text(face = "bold", size = base + 3),
          plot.subtitle = element_text(size = base - 0.5, colour = "grey30", lineheight = 1.25),
          plot.caption = element_text(size = base - 1.5, colour = "grey45", hjust = 0),
          strip.text = element_text(face = "bold", size = base - 1),
          panel.grid.minor = element_blank(), legend.position = "top")
}
.jm_grp  = c("young", "medium", "elderly")
.jm_gcol = c(young = "#1B7F5C", medium = "#2B5D8A", elderly = "#C1541E")
.jm_blue = "#2B5D8A"; .jm_orange = "#C1541E"; .jm_green = "#1B7F5C"
# pull the interval rows for a set of parameter names, in that order
.jm_iv = function(iv, names_wanted){
  if (is.null(iv)) return(NULL)
  i = match(names_wanted, iv$parameter)
  data.frame(estimate = iv$estimate[i], lower = iv$lower[i], upper = iv$upper[i])
}
.jm_ivnote = function(iv) if (is.null(iv)) "" else
  "\nBars are 95% intervals from the curvature of the fitted surface."

# ================= HOW THE MODEL WORKS =================

# ---- |-01 what varies where ----
plot_jm_design = function(fit){
  d = fit$d; S = d$n_season; C = d$n_country
  spec = tibble::tribble(
    ~group,        ~parameter,                  ~season, ~country, ~age, ~n,                       ~note,
    "dynamics",    "R0  transmissibility",      TRUE,  FALSE, FALSE, S,           "one number per season for all of Europe",
    "dynamics",    "S0  susceptibility",        FALSE, TRUE,  FALSE, C,           "one per country, same across its seasons",
    "dynamics",    "sigma  elderly suscept.",   FALSE, FALSE, TRUE,  1,           "one number for everyone",
    "dynamics",    "I0  seed / arrival",        TRUE,  TRUE,  FALSE, d$n_cs,      "free for every wave: sets when it arrives",
    "observation", "c  reporting level",        FALSE, TRUE,  FALSE, C,           "one per surveillance system",
    "observation", "delta  season visibility",  TRUE,  FALSE, FALSE, S,           "shared, constrained to average one",
    "observation", "off  age reporting",        FALSE, TRUE,  TRUE,  2 * C,       "adults the reference",
    "observation", "b  off-season baseline",    FALSE, TRUE,  FALSE, sum(d$n_src),"one per data source present",
    "observation", "phi  dispersion",           FALSE, TRUE,  FALSE, C,           "one per country",
    "fixed",       "gamma  infectious period",  FALSE, FALSE, FALSE, 0,           "3.6 days, from the literature",
    "fixed",       "vaccine effects (3)",       FALSE, FALSE, TRUE,  0,           "fixed, 65+ pulse on 1 October",
    "fixed",       "contact matrix",            FALSE, TRUE,  TRUE,  0,           "fixed, rescaled to spectral radius 1")
  long = spec %>%
    pivot_longer(c(season, country, age), names_to = "dim", values_to = "varies") %>%
    mutate(dim = factor(dim, levels = c("season", "country", "age"),
                        labels = c("by season", "by country", "by age")),
           parameter = factor(parameter, levels = rev(spec$parameter)),
           state = ifelse(spec$n[match(as.character(parameter), spec$parameter)] == 0, "fixed",
                   ifelse(varies, "varies", "same")))
  lab = spec %>% mutate(parameter = factor(parameter, levels = rev(spec$parameter)),
                        txt = ifelse(n == 0, "fixed", paste0("x", n)))
  ggplot(long, aes(dim, parameter)) +
    geom_tile(aes(fill = state), colour = "white", linewidth = 1.1) +
    geom_text(data = lab, aes(x = 3.95, y = parameter, label = txt), hjust = 0,
              size = 3.1, colour = "grey25", fontface = "bold") +
    geom_text(data = lab, aes(x = 4.65, y = parameter, label = note), hjust = 0,
              size = 2.9, colour = "grey40") +
    scale_fill_manual(values = c(varies = .jm_blue, same = "grey88", fixed = "grey70"),
                      labels = c(varies = "varies along this dimension",
                                 same = "assumed the same along it", fixed = "not fitted at all"),
                      name = NULL) +
    scale_x_discrete(expand = expansion(add = c(0.5, 9))) +
    labs(title = "What varies where: the whole model in one picture",
         subtitle = paste("Each row is a parameter; the three squares say whether it is allowed to differ between seasons,",
                          "between countries, or between age groups. The count on the right is how many numbers that row",
                          "contributes. The design is the science: sharing transmissibility across countries is what makes",
                          "susceptibility identifiable, and constraining season visibility to average one is what keeps it",
                          "separable from each country's reporting level.", sep = "\n"),
         x = NULL, y = NULL) +
    .jm_theme() + theme(panel.grid = element_blank(), axis.text.y = element_text(size = 9.5))
}

# ---- |-02 how a wave's arrival is set ----
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
    labs(title = "How the model decides when a wave arrives",
         subtitle = paste("The seed: the infected fraction planted on 1 August. It has no epidemiological meaning of its own.",
                          "A seed ten times smaller simply makes that wave arrive about two weeks later, and it is the only",
                          "parameter that can move a wave sideways without changing its shape. That is why it is free for",
                          "every country-season: waves genuinely arrive at different times and nothing else can express it.",
                          "Darker orange = smaller seed = later wave.", sep = "\n"),
         x = NULL, y = NULL) +
    .jm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ================= WHAT DATA IT FITS =================

# ---- |-03 every country-season ----
plot_jm_fit_overview = function(fit){
  d = fit$d; f = jm_fitted_cpp(fit$theta, d)
  dd = do.call(rbind, lapply(seq_len(d$n_cs), function(i){
    ic = d$cs_country[i] + 1L; N = d$N[[ic]]; w = matrix(N / sum(N), nrow(d$y[[i]]), 3, byrow = TRUE)
    data.frame(country = d$countries[ic], season = d$seasons[d$cs_season[i] + 1L],
               week = seq_len(nrow(d$y[[i]])),
               observed = rowSums(sweep(d$y[[i]], 2, d$rate_per / N, "*") * w),
               model    = rowSums(sweep(f$mu[[i]], 2, d$rate_per / N, "*") * w))
  }))
  ggplot(dd, aes(week)) +
    geom_point(aes(y = observed), size = 0.5, colour = "grey25", na.rm = TRUE) +
    geom_line(aes(y = model), colour = .jm_orange, linewidth = 0.65) +
    facet_grid(country ~ season, scales = "free_y") +
    labs(title = "Every country-season: the model against the data",
         subtitle = paste("Points are observed influenza-positive ILI consultations per 100 000, pooled over age groups by",
                          "population weight; the line is the model at the joint optimum. One transmissibility and one",
                          "visibility number per season must serve every country at once, so the line is not free to chase",
                          "each panel. A panel that misses is therefore information about the shared assumption, not just",
                          "about that country.", sep = "\n"),
         x = "week of the season (from 1 August)", y = "ILI+ per 100 000") +
    .jm_theme(9) + theme(axis.text = element_text(size = 6.5), panel.spacing = unit(0.25, "lines"))
}

# ---- |-04 one country, by age ----
plot_jm_fit_country = function(fit, cc){
  dd = jm_tidy_fit(fit) %>% filter(country == cc) %>% mutate(group = factor(group, levels = .jm_grp))
  ggplot(dd, aes(week, value, colour = group)) +
    geom_point(data = dd %>% filter(what == "observed"), size = 0.6, alpha = 0.8, na.rm = TRUE) +
    geom_line(data = dd %>% filter(what == "model"), linewidth = 0.7) +
    facet_grid(group ~ season, scales = "free_y") +
    scale_colour_manual(values = .jm_gcol, guide = "none") +
    labs(title = paste0(cc, ": the three age groups, season by season"),
         subtitle = paste("Points observed, lines the model. The age profile comes from the contact matrix, the one global",
                          "elderly susceptibility and this country's two age reporting offsets. None of those varies by",
                          "season, so every season-to-season change you see is produced by the shared transmissibility,",
                          "the shared season visibility and this wave's seed.", sep = "\n"),
         x = "week of the season", y = "ILI+ per 100 000 of the age group") +
    .jm_theme(9) + theme(axis.text = element_text(size = 6.5))
}

# ---- |-05 the noise budget ----
plot_jm_adequacy = function(fit){
  a = jm_adequacy(fit)
  long = a %>% select(country, model = cv_fitted, data = cv_data) %>%
    pivot_longer(c(model, data), names_to = "what", values_to = "cv")
  ggplot(long, aes(cv, reorder(country, cv), colour = what)) +
    geom_line(aes(group = country), colour = "grey80", linewidth = 1) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = c(model = .jm_orange, data = .jm_blue),
                        labels = c(data = "the data's own week-to-week scatter", model = "what the fit needed"),
                        name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "The noise budget: how much noise the fit needed against how much the data have",
         subtitle = paste("Blue is the observed scatter of each weekly series around its own 3-week average, a floor no mean",
                          "could beat, so it is pure measurement noise. Orange is the noise the fitted dispersion implies.",
                          "The gap is the deterministic curve failing to follow the wave, written off as if it were",
                          "measurement error. This is the model's honest limitation and the trigger for putting the Kalman",
                          "filter back: if the filter closes this gap it has earned its place, and if not the missing",
                          "structure is something else.", sep = "\n"),
         x = "coefficient of variation of a weekly count", y = NULL) + .jm_theme()
}

# ================= WHAT IT LEARNS =================

# ---- |-06 season transmissibility ----
plot_jm_season_R0 = function(fit, iv = NULL){
  d = fit$d; s = jm_summary_season(fit)
  s$season = factor(s$season, levels = d$seasons)
  ci = .jm_iv(iv, paste0("log_R0_", d$seasons)); if (!is.null(ci)) s = cbind(s, ci[, c("lower", "upper")])
  g = ggplot(s, aes(season, R0, group = 1)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5 * exp(-1.96 * d$pr_R0_sd),
             ymax = 1.5 * exp(1.96 * d$pr_R0_sd), fill = .jm_blue, alpha = 0.10) +
    geom_hline(yintercept = 1.5, linetype = "dashed", colour = "grey50")
  if (!is.null(ci)) g = g + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, colour = .jm_blue)
  g + geom_line(colour = .jm_blue, linewidth = 0.8) + geom_point(size = 3, colour = .jm_blue) +
    geom_text(aes(label = sprintf("%.2f", R0)), vjust = -1.4, size = 3, colour = "grey20") +
    labs(title = "What it learns, 1: how transmissible each season's virus was",
         subtitle = paste0("R0 is how many people one case would infect in a fully susceptible population, given the age mixing.",
                          "\nIt is ONE number per season for the whole of Europe, carried by every wave in that season. The shaded",
                          "\nband is the prior's 95% range around 1.5 and the dashed line its centre, so movement away from the",
                          "\ncentre is the data speaking against the prior's pull.", .jm_ivnote(iv)),
         x = NULL, y = expression(R[0]~"of the season")) +
    .jm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ---- |-07 season visibility ----
plot_jm_season_visibility = function(fit, iv = NULL){
  d = fit$d; s = jm_summary_season(fit)
  s$season = factor(s$season, levels = d$seasons)
  ci = .jm_iv(iv, paste0("delta_", d$seasons))
  if (!is.null(ci)){ s$lower = exp(ci$lower); s$upper = exp(ci$upper) }
  g = ggplot(s, aes(season, reporting_mult, group = 1)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50")
  if (!is.null(ci)) g = g + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, colour = .jm_orange)
  g + geom_line(colour = .jm_orange, linewidth = 0.8) + geom_point(size = 3, colour = .jm_orange) +
    geom_text(aes(label = sprintf("%.2f", reporting_mult)), vjust = -1.4, size = 3, colour = "grey20") +
    scale_y_log10() +
    labs(title = "What it learns, 2: how visible each season was, per infection",
         subtitle = paste0("Season visibility: positive consultations per infection that season, relative to normal. 1.3 means a",
                          "\nthird more consultations for the same number of infections. It mixes how symptomatic the season's strain",
                          "\nwas with any season size the mechanism cannot produce, so it is partly a residual and must be read with",
                          "\nthe noise budget in mind. Constrained to average one, which is what makes it separable from each",
                          "\ncountry's own reporting level.", .jm_ivnote(iv)),
         x = NULL, y = "multiplier on reporting (log scale)") +
    .jm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ---- |-08 country susceptibility ----
plot_jm_country_S0 = function(fit, iv = NULL){
  d = fit$d; s = jm_summary_country(fit)
  ci = .jm_iv(iv, paste0(d$countries, ":logit_S0")); if (!is.null(ci)) s = cbind(s, ci[, c("lower", "upper")])
  g = ggplot(s, aes(S0, reorder(country, S0))) +
    annotate("rect", xmin = plogis(d$pr_S0_mean - 1.96 * d$pr_S0_sd),
             xmax = plogis(d$pr_S0_mean + 1.96 * d$pr_S0_sd), ymin = -Inf, ymax = Inf,
             fill = .jm_blue, alpha = 0.08) +
    geom_vline(xintercept = plogis(d$pr_S0_mean), linetype = "dashed", colour = "grey50")
  if (!is.null(ci)) g = g + geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.22, colour = .jm_blue)
  g + geom_point(size = 3.2, colour = .jm_blue) +
    geom_text(aes(label = sprintf("%.3f", S0)), hjust = -0.35, size = 3, colour = "grey20") +
    scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
    labs(title = "What it learns, 3: how susceptible each country was at the season start",
         subtitle = paste0("S0 is the share of the population that could be infected on 1 August, held the same across that",
                          "\ncountry's seasons. It is the project's target quantity, and it is identifiable here only because",
                          "\ntransmissibility is shared across countries: fitted one country at a time, the two trade off almost",
                          "\nperfectly and the answer comes from the prior. Shaded band is the prior's 95% range. Trust the",
                          "\nRANKING more than the absolute level.", .jm_ivnote(iv)),
         x = "S0", y = NULL) + .jm_theme()
}

# ---- |-09 reporting level ----
plot_jm_country_reporting = function(fit, iv = NULL){
  d = fit$d; s = jm_summary_country(fit)
  ci = .jm_iv(iv, paste0(d$countries, ":log_c"))
  if (!is.null(ci)){ s$lower = 100 * ci$lower; s$upper = 100 * ci$upper }
  g = ggplot(s, aes(c_adult * 100, reorder(country, c_adult)))
  if (!is.null(ci)) g = g + geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.22, colour = .jm_green)
  g + geom_point(size = 3.2, colour = .jm_green) +
    geom_text(aes(label = sprintf("%.1f%%", c_adult * 100)), hjust = -0.35, size = 3, colour = "grey20") +
    scale_x_log10(expand = expansion(mult = c(0.08, 0.22))) +
    labs(title = "What it learns, 4: what fraction of adult infections is actually counted",
         subtitle = paste0("The chance that one adult infection becomes an influenza-positive ILI consultation in that country's",
                          "\nsurveillance. A country low on this axis is not having fewer infections, it is seeing fewer of them.",
                          "\nThis is why raw ILI+ curves cannot be compared between countries and these numbers are what make",
                          "\nthem comparable. Log scale: the spread is more than tenfold.", .jm_ivnote(iv)),
         x = "per cent of adult infections reported (log scale)", y = NULL) + .jm_theme()
}

# ---- |-10 age reporting ----
plot_jm_age_offsets = function(fit, iv = NULL){
  d = fit$d; s = jm_summary_country(fit)
  age = s %>% select(country, young = rel_young, elderly = rel_elderly) %>%
    pivot_longer(c(young, elderly), names_to = "group", values_to = "rel") %>%
    mutate(group = factor(group, levels = c("young", "elderly")))
  ci = .jm_iv(iv, c(paste0(d$countries, ":off_young"), paste0(d$countries, ":off_eld")))
  if (!is.null(ci)){ age$lower = ci$lower; age$upper = ci$upper }
  g = ggplot(age, aes(rel, reorder(country, rel), colour = group)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50")
  if (!is.null(ci)) g = g + geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.5,
                                           position = position_dodge(width = 0.5))
  g + geom_point(size = 3, position = position_dodge(width = 0.5)) + scale_x_log10() +
    scale_colour_manual(values = .jm_gcol[c("young", "elderly")], name = NULL) +
    labs(title = "What it learns, 5: how much more visible a child or an elderly person is, per infection",
         subtitle = paste0("Age reporting offsets relative to adults, one pair per country. Above one means that group generates",
                          "\nmore positive consultations per infection than an adult does: they consult more readily, or the",
                          "\nsentinel network sees them more. This is surveillance, not biology, which is exactly why it is allowed",
                          "\nto differ by country while the elderly's susceptibility per contact is one number for all of Europe.",
                          .jm_ivnote(iv)),
         x = "reporting relative to adults (log scale)", y = NULL) + .jm_theme()
}

# ---- |-11 attack rates ----
plot_jm_attack = function(fit){
  d = fit$d; f = jm_fitted_cpp(fit$theta, d)
  att = as.data.frame(f$attack); names(att) = d$groups
  att$country = d$countries[d$cs_country + 1L]; att$season = d$seasons[d$cs_season + 1L]
  long = att %>% pivot_longer(all_of(d$groups), names_to = "group", values_to = "attack") %>%
    mutate(group = factor(group, levels = .jm_grp), season = factor(season, levels = d$seasons))
  p = jm_unpack(fit$theta, d)
  ggplot(long, aes(season, attack, colour = group, group = interaction(country, group))) +
    geom_line(alpha = 0.32, linewidth = 0.5) +
    stat_summary(aes(group = group), fun = median, geom = "line", linewidth = 1.5) +
    scale_colour_manual(values = .jm_gcol, name = NULL) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = "What it learns, 6: who actually gets infected",
         subtitle = paste0("Modelled share of each age group infected over the season. Thin lines are countries, thick lines the",
                          "\nmedian across them. This is the one output entirely free of reporting: it comes from the contact",
                          "\nmatrix, the global elderly susceptibility (fitted at ", sprintf("%.2f", p$sigma_eld),
                          "x an adult's) and vaccination, not from how",
                          "\nmany consultations were counted. So compare it with prospective cohort evidence, never with",
                          "\nsurveillance curves."),
         x = NULL, y = "attack rate") +
    .jm_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ================= WHETHER TO BELIEVE IT =================

# ---- |-12 data or prior ----
plot_jm_identifiability = function(id){
  f = id$family %>% mutate(family = reorder(family, contraction_med))
  ggplot(f, aes(contraction_med, family)) +
    geom_segment(aes(x = 0, xend = contraction_med, yend = family), colour = "grey80", linewidth = 1) +
    geom_point(aes(size = n), colour = .jm_blue) +
    geom_text(aes(label = sprintf("%.2f", contraction_med)), hjust = -0.5, size = 3, colour = "grey20") +
    scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
    scale_size_continuous(range = c(2, 5.5), name = "parameters") +
    labs(title = "Can you believe it, 1: did the data decide, or the prior?",
         subtitle = paste("Prior-to-posterior contraction. Zero means the fit handed back the prior unchanged, so that number",
                          "carries no information from the surveillance data at all. One means the data determined it and the",
                          "prior was irrelevant. Anything below about 0.3 should be reported as an assumption rather than a",
                          "result. Dot size is how many parameters of that kind there are.", sep = "\n"),
         x = "contraction (0 = all prior, 1 = all data)", y = NULL) + .jm_theme()
}

# ---- |-13 recovery of a known truth ----
plot_jm_recovery = function(rec, d, summ = NULL){
  if (is.null(summ)) summ = jm_recovery_summary(rec, d)
  cmp = summ$comparison %>%
    filter(family %in% c("R0 (season, shared)", "season deviation (shared)", "S0 (country)",
                         "reporting c (country)", "elderly susceptibility (global)"))
  cov_txt = summ$by_family %>% filter(!is.na(coverage)) %>%
    summarise(m = median(coverage)) %>% pull(m)
  # A replicate in which one country fell into the flat-line optimum throws an estimate far off scale
  # and would squash every other panel flat. Such points are WINSORISED FOR DISPLAY ONLY, drawn as
  # open triangles at the panel edge and counted in the subtitle, so they are visible rather than
  # hidden and the informative range stays readable.
  cmp = cmp %>% group_by(family) %>%
    mutate(rng = diff(range(truth)),
           cap_hi = max(truth) + 0.35 * ifelse(rng > 0, rng, abs(max(truth))),
           cap_lo = min(truth) - 0.35 * ifelse(rng > 0, rng, abs(max(truth))),
           off = estimate > cap_hi | estimate < cap_lo,
           shown = pmin(pmax(estimate, cap_lo), cap_hi)) %>% ungroup()
  n_off = sum(cmp$off)
  ggplot(cmp, aes(truth, shown)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey45", linetype = "dashed") +
    geom_point(data = cmp %>% filter(!off), alpha = 0.55, size = 1.7, colour = .jm_blue) +
    geom_point(data = cmp %>% filter(off), shape = 2, size = 2.2, stroke = 0.8, colour = .jm_orange) +
    facet_wrap(~ family, scales = "free") +
    labs(title = "Can you believe it, 2: recovering a known truth",
         subtitle = paste0("Data were simulated from the model itself, at a known parameter set, onto the REAL design: the same",
                          "\ncountries, seasons, observed weeks, missing cells and populations as the actual panel. Then the whole",
                          "\npipeline was refitted from scratch. Each point is one parameter in one replicate; the dashed line is",
                          "\nperfect recovery. Season transmissibility rank correlation ",
                          sprintf("%.2f", summ$rank_shared$spearman_med[summ$rank_shared$block == "R0 by season"][1]),
                          ", susceptibility ranking ", sprintf("%.2f", summ$rank_S0_spearman),
                          ",\n95% interval coverage ", sprintf("%.0f%%", 100 * cov_txt),
                          ". Orange triangles are ", n_off, " estimate(s) off scale, where one country",
                          "\nin one replicate fell into the flat-line optimum; they are drawn at the panel edge, not dropped.",
                          "\nAnything the fit cannot recover from its own simulation cannot be trusted from real data either."),
         x = "true value", y = "estimated value") + .jm_theme()
}

# ================= write the default set =================
save_jm_report = function(fit, id = NULL, iv = NULL, rec = NULL, dir = "output/joint_model"){
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  d = fit$d; man = character(0)
  put = function(file, plot, w, h, what){
    ggsave(file.path(dir, file), plot, width = w, height = h, dpi = 115, limitsize = FALSE)
    man <<- c(man, sprintf("| `%s` | %s |", file, what))
  }
  put("01_design_what_varies_where.png", plot_jm_design(fit), 12.5, 6.2, "how the model works: which parameters vary by season, by country, by age, and which are fixed")
  put("02_arrival_time.png", plot_jm_arrival(fit), 10, 5.8, "how the model sets when each wave arrives (the seed)")
  put("03_fit_overview.png", plot_jm_fit_overview(fit), 2.05 * d$n_season + 2, 1.35 * d$n_country + 2.8,
      "what data it fits: every country-season, observed against modelled")
  for (cc in d$countries)
    ggsave(file.path(dir, sprintf("04_fit_%s.png", cc)), plot_jm_fit_country(fit, cc),
           width = 2.05 * d$n_season + 2, height = 7.2, dpi = 110, limitsize = FALSE)
  man = c(man, sprintf("| `04_fit_<country>.png` | the same by age group, one file per country (%s) |",
                       paste(d$countries, collapse = ", ")))
  put("05_noise_budget.png", plot_jm_adequacy(fit), 10, 5.8, "how well it fits, honestly: the noise the fit needed against the noise the data have")
  put("06_season_R0.png", plot_jm_season_R0(fit, iv), 10, 6.0, "what it learns: transmissibility of each season's virus, shared across countries")
  put("07_season_visibility.png", plot_jm_season_visibility(fit, iv), 10, 6.0, "what it learns: how visible each season was per infection, shared across countries")
  put("08_country_S0.png", plot_jm_country_S0(fit, iv), 10, 6.0, "what it learns: susceptibility at the season start, by country")
  put("09_country_reporting.png", plot_jm_country_reporting(fit, iv), 10, 5.6, "what it learns: fraction of adult infections that is counted, by country")
  put("10_age_reporting.png", plot_jm_age_offsets(fit, iv), 10, 5.8, "what it learns: how visible children and the elderly are per infection, by country")
  put("11_attack_rates.png", plot_jm_attack(fit), 10, 6.0, "what it learns: modelled attack rate by age group, the one reporting-free output")
  if (!is.null(id)) put("12_data_or_prior.png", plot_jm_identifiability(id), 10, 5.4,
                        "whether to believe it: how much each parameter owes to the data rather than its prior")
  if (!is.null(rec)) put("13_recovery.png", plot_jm_recovery(rec, d), 10, 7.0,
                         "whether to believe it: recovery of a known truth simulated from the model onto the real design")
  write.csv(jm_summary_season(fit), file.path(dir, "summary_season.csv"), row.names = FALSE)
  write.csv(jm_summary_country(fit), file.path(dir, "summary_country.csv"), row.names = FALSE)
  write.csv(jm_adequacy(fit), file.path(dir, "noise_budget.csv"), row.names = FALSE)
  if (!is.null(iv)) write.csv(iv, file.path(dir, "intervals.csv"), row.names = FALSE)
  if (!is.null(id)) write.csv(id$table, file.path(dir, "identifiability.csv"), row.names = FALSE)
  writeLines(c("# The default figure set", "",
               "Walk them in order: how the model works, what data it fits, what it learns, whether to believe it.",
               "Every figure carries its own interpretation in the subtitle.", "",
               "| file | what it shows |", "|---|---|", man, "",
               sprintf("Regenerate with `Rscript code/07_joint_model/run_joint_model.R` (%d countries, %d seasons, %d parameters).",
                       d$n_country, d$n_season, d$n_par)),
             file.path(dir, "MANIFEST.md"))
  invisible(dir)
}
