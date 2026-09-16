# joint_report.R -- THE DEFAULT FIGURE SET for the joint model.
#
# The set is numbered and answers five questions in order, so a reader can walk it front to back:
#   WHAT DATA THERE IS        01 the panel, 02 the four features that dictate the model
#   HOW THE MODEL WORKS       03 what each parameter does to a wave, 04 what varies where,
#                             05 how a wave's arrival is set
#   WHAT DATA IT FITS         06 every country-season, 07 one country by age, 08 the noise budget
#   WHAT IT LEARNS            09 season transmissibility, 10 season visibility, 11 country
#                             susceptibility, 12 reporting level, 13 age reporting, 14 attack rates
#   WHETHER TO BELIEVE IT     15 data or prior, 16 recovery of a known truth
# Every figure carries its own interpretation in the subtitle: what the quantity means, and what the
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
  # A name that is not in the interval table yields an NA row, which ggplot drops SILENTLY -- the bar
  # just is not there and nothing says so. Warn, so a renamed parameter cannot quietly remove the
  # uncertainty from a figure. The returned frame keeps the requested ORDER, which is what callers
  # rely on when they bind it to a plotting frame.
  if (anyNA(i)) warning(sprintf("no interval for %d parameter(s), bars omitted: %s",
                                sum(is.na(i)), paste(head(names_wanted[is.na(i)], 5), collapse = ", ")),
                        call. = FALSE)
  data.frame(estimate = iv$estimate[i], lower = iv$lower[i], upper = iv$upper[i])
}
.jm_ivnote = function(iv) if (is.null(iv)) "" else
  "\nBars are 95% intervals from the curvature of the fitted surface."

# ================= WHAT DATA THERE IS =================
#
# The panel comes before the model in the set deliberately: every design choice the model makes is
# forced by a feature of this data, and a reader who has seen the data first recognises each choice
# instead of taking it on trust.

# a tidy long frame of the OBSERVED series, pooled over age by population weight -- the same quantity
# figure 09 plots the model against, but with no model in it
.jm_obs = function(d){
  do.call(rbind, lapply(seq_len(d$n_cs), function(i){
    ic = d$cs_country[i] + 1L; N = d$N[[ic]]
    per = d$rate_per / N
    y = d$y[[i]]
    data.frame(country = d$countries[ic], season = d$seasons[d$cs_season[i] + 1L],
               week = seq_len(nrow(y)),
               # population-weighted ILI+ per 100 000, i.e. sum(counts) / sum(population)
               value = rowSums(sweep(y, 2, per, "*") * matrix(N / sum(N), nrow(y), length(N), byrow = TRUE)),
               stringsAsFactors = FALSE)
  }))
}

# ---- |-01 the panel: what data there is ----
plot_jm_data_panel = function(fit){
  d = fit$d
  obs = .jm_obs(d)
  cs = obs %>% group_by(country, season) %>%
    summarise(weeks = sum(is.finite(value)), peak = suppressWarnings(max(value, na.rm = TRUE)),
              .groups = "drop")
  src = data.frame(country = d$countries[d$cs_country + 1L],
                   season = d$seasons[d$cs_season + 1L],
                   source = vapply(seq_len(d$n_cs), function(i)
                     d$sources[[d$cs_country[i] + 1L]][d$cs_src[i] + 1L], character(1)),
                   stringsAsFactors = FALSE)
  cs = left_join(cs, src, by = c("country", "season"))
  # every country x season the design COULD have had, so absent cells are visible as absent
  grid = expand.grid(country = d$countries, season = d$seasons, stringsAsFactors = FALSE)
  ex = d$excluded_ambiguous
  exkey = if (!is.null(ex) && nrow(ex)) paste(ex$country, ex$season) else character(0)
  g = left_join(grid, cs, by = c("country", "season")) %>%
    mutate(season = factor(season, levels = d$seasons),
           country = factor(country, levels = rev(d$countries)),
           # an EXCLUDED cell is not an absent one: the country reported that season, we chose not to
           # fit it. Showing both as blank would hide a decision behind a data gap.
           state = ifelse(paste(country, season) %in% exkey, "excluded",
                          ifelse(is.na(weeks), "no data", source)))
  exd = if (length(exkey)) g %>% filter(state == "excluded") else NULL
  ggplot(g, aes(season, country)) +
    geom_tile(aes(fill = state), colour = "white", linewidth = 0.9) +
    geom_text(aes(label = ifelse(is.na(weeks), "", as.character(weeks))), size = 2.9,
              colour = "grey15") +
    { if (!is.null(exd)) geom_point(data = exd, shape = 4, size = 6, stroke = 1.4,
                                    colour = "#B02020") } +
    scale_fill_manual(values = c(ERVISS = "#7FA8C9", RespiCompass = "#BFD8B0",
                                 excluded = "#F2C9C9", `no data` = "grey93"),
                      na.value = "grey93", name = NULL,
                      breaks = c("ERVISS", "RespiCompass", "excluded", "no data"),
                      labels = c("ERVISS", "RespiCompass", "excluded (see below)", "no data")) +
    labs(title = "The panel: 85 country-seasons of weekly influenza-positive ILI",
         subtitle = paste("One tile per country and season; the number is how many weeks carry an observation, and the",
                          "colour which surveillance source that season came from. Blank means the country reported",
                          "nothing usable that season. The red cross is the one country-season excluded because its",
                          "influenza positivity cannot be computed for 14 weeks (see MODEL.md). Week counts differ",
                          "because reporting stops at different points in the year -- which is why the attack rate is",
                          "integrated to a fixed horizon rather than to the end of each series.", sep = "\n"),
         caption = sprintf("%d countries x %d seasons; %d country-seasons fitted, %s observed age-week cells. RespiCompass covers up to 2023/24, ERVISS from 2024/25.",
                           d$n_country, d$n_season, d$n_cs,
                           format(sum(vapply(d$y, function(m) sum(is.finite(m)), numeric(1))), big.mark = ",")),
         x = NULL, y = NULL) +
    .jm_theme() + theme(panel.grid = element_blank(),
                        axis.text.x = element_text(angle = 30, hjust = 1))
}

# ---- |-02 the four features of the data that dictate the model ----
# Every one of the model's sharing choices is forced by one of these. Putting them in one figure is
# the shortest honest answer to "why is the model built this way".
plot_jm_data_features = function(fit){
  d = fit$d
  obs = .jm_obs(d)
  pk = obs %>% group_by(country, season) %>%
    summarise(peak = suppressWarnings(max(value, na.rm = TRUE)),
              peak_week = week[which.max(replace(value, is.na(value), -Inf))], .groups = "drop") %>%
    filter(is.finite(peak), peak > 0)

  # (a) the between-country spread in observed level, and that it is a COUNTRY property
  a = pk %>% mutate(country = reorder(country, peak, median))
  pa = ggplot(a, aes(peak, country)) +
    geom_line(aes(group = country), colour = "grey85", linewidth = 2.6) +
    geom_point(colour = .jm_blue, size = 1.9, alpha = 0.85) +
    scale_x_log10() +
    labs(title = "1. Observed level is a COUNTRY property, not an epidemic one",
         subtitle = paste0("Peak weekly ILI+ per 100 000, one point per season. The spread between countries is about ",
                           sprintf("%.0f", max(a$peak) / min(a$peak)),
                           "-fold\nand a country keeps its place across seasons, so it cannot be how many people were",
                           " infected.\nIt is how many infections that country's surveillance SEES -> one reporting level per country."),
         x = "peak ILI+ per 100 000 (log scale)", y = NULL)

  # (b) seasons move TOGETHER across countries -> shared season parameters
  b = pk %>% group_by(country) %>% mutate(rel = peak / median(peak)) %>% ungroup() %>%
    mutate(season = factor(season, levels = d$seasons))
  pb = ggplot(b, aes(season, rel)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_line(aes(group = country), colour = "grey80", linewidth = 0.5) +
    stat_summary(fun = median, geom = "line", aes(group = 1), colour = .jm_orange, linewidth = 1.6) +
    stat_summary(fun = median, geom = "point", colour = .jm_orange, size = 2.6) +
    scale_y_log10() +
    labs(title = "2. Seasons rise and fall TOGETHER across Europe",
         subtitle = paste("Each country's peak divided by its own median peak, so the country level is removed. Thin",
                          "lines are countries, thick the median. They move in step -- a big season is big almost",
                          "everywhere -- which is what licenses ONE transmissibility and ONE visibility per season", sep = "\n"),
         x = NULL, y = "peak relative to that country's median")

  # (c) a third of the observations are exactly zero -> negative binomial, not Gaussian
  cnt = unlist(lapply(d$y, function(m) m[is.finite(m)]))
  cz = data.frame(count = cnt) %>%
    mutate(bin = cut(count, c(-1, 0, 1, 3, 10, 30, 100, 300, Inf),
                     labels = c("0", "1", "2-3", "4-10", "11-30", "31-100", "101-300", ">300"))) %>%
    count(bin)
  cz$pct_lab = ifelse(cz$bin == "0", sprintf("%.0f%%", 100 * cz$n / sum(cz$n)), "")
  pc = ggplot(cz, aes(bin, n)) +
    geom_col(fill = .jm_green, width = 0.75) +
    geom_text(aes(label = pct_lab), vjust = -0.5, size = 3.4, fontface = "bold", colour = .jm_green) +
    labs(title = "3. A third of the observations are exactly zero",
         subtitle = paste("Weekly influenza-positive consultations, all age groups and country-seasons pooled. A",
                          "Gaussian around a small mean puts a sixth of its mass below zero and cannot represent",
                          "this spike, so the observation model is negative-binomial on counts.", sep = "\n"),
         x = "weekly count in one age group", y = "observations")

  # (d) waves arrive at different times -> a free seed per country-season
  pd = ggplot(pk, aes(peak_week)) +
    geom_histogram(binwidth = 1, fill = .jm_blue, colour = "white", linewidth = 0.2) +
    labs(title = "4. Waves arrive weeks apart, even in the same season",
         subtitle = paste0("Week of the observed peak, counted from 1 August, over all ", nrow(pk),
                           " country-seasons: a spread of ",
                           diff(range(pk$peak_week)), " weeks.\nNothing else in the model can move a wave sideways",
                           " without changing its shape, so the seed size\nis left free for every country-season."),
         x = "week of the observed peak (from 1 August)", y = "country-seasons")

  th = .jm_theme(9) + theme(plot.title = element_text(face = "bold", size = 11),
                            plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.2))
  patchwork::wrap_plots(pa + th, pb + th, pc + th, pd + th, ncol = 2) +
    patchwork::plot_annotation(
      title = "Four features of the data, and the model choice each one forces",
      subtitle = "Read this before the model. Every sharing decision in the design answers one of these; none of them is a modelling preference.",
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, colour = "grey30")))
}

# ================= HOW THE MODEL WORKS =================

# ---- |-03 what each parameter DOES to a wave ----
# The design table (figure 07) says where each parameter varies; this says what it MEANS, by showing
# the curve the model actually produces when you move it. Every curve here comes from
# jm_fitted_cpp at a perturbed parameter vector -- the real C++ model, not a redrawing of it -- so
# the figure cannot drift from the thing it describes.
#
# The last panel is the point of the whole design: two pairs of curves that lie on top of each other.
# Transmissibility and susceptibility enter the rise rate as a product, and reporting level and
# season visibility enter the observation as a product, so within one country-season each pair is
# indistinguishable. Sharing R0 across countries and constraining the visibilities to average one is
# what breaks the two ties.
plot_jm_mechanism = function(fit, ref = NULL){
  d = fit$d; th = fit$theta; S = d$n_season
  # a reference wave with a full grid and a clear peak, and NOT in the last season, whose visibility
  # is the constrained one rather than a free slot
  cand = which(d$cs_season < S - 1L & d$n_weeks >= 45)
  if (!length(cand)) cand = which(d$cs_season < S - 1L)
  pk = vapply(cand, function(i) suppressWarnings(max(d$y[[i]][, 2], na.rm = TRUE)), numeric(1))
  i0 = if (is.null(ref)) cand[which.max(pk)] else ref
  ic = d$cs_country[i0] + 1L; s_ix = d$cs_season[i0] + 1L
  cc = d$countries[ic]; ss = d$seasons[s_ix]
  base = d$off_country[ic]; nsrc = d$n_src[ic]
  nm = jm_par_names(d)
  # index of each slot we perturb
  j = list(R0 = s_ix, delta = if (s_ix <= S - 1L) S + s_ix else NA_integer_,
           S0 = base + 1L, c = base + 2L, off_eld = base + 4L, phi = base + 5L,
           I0 = base + 5L + nsrc + d$cs_pos[i0] + 1L)
  N = d$N[[ic]]; per = d$rate_per / N
  mu_of = function(theta, grp = 2L){
    m = jm_fitted_cpp(theta, d)$mu[[i0]]
    data.frame(week = seq_len(nrow(m)), value = m[, grp] * per[grp])
  }
  bump = function(slot, by){ t2 = th; t2[slot] = t2[slot] + by; t2 }

  # each panel: the reference curve plus one lower and one higher setting of a single parameter
  panels = list(
    list(key = "R0  transmissibility of the season", lo = bump(j$R0, log(0.93)), hi = bump(j$R0, log(1.07)),
         lab = c("-7%", "+7%"),
         note = "Steeper AND taller, and it peaks earlier. A 7% move is enough to double the peak, which is how little of this the data need to see."),
    list(key = "S0  susceptibility of the country", lo = bump(j$S0, -0.25), hi = bump(j$S0, 0.25),
         lab = c("lower", "higher"),
         note = "Does the SAME thing as R0: the two enter the rise rate as a product. See the bottom panel."),
    list(key = "I0  seed, i.e. when the wave arrives", lo = bump(j$I0, log(0.01)), hi = bump(j$I0, log(100)),
         lab = c("/100", "x100"),
         note = "Moves the wave sideways and changes nothing else. The only handle on timing."),
    list(key = "c  reporting level of the country", lo = bump(j$c, log(0.6)), hi = bump(j$c, log(1.6)),
         lab = c("-40%", "+60%"),
         note = "Scales the whole curve. Same infections, more or fewer of them counted."),
    list(key = "delta  visibility of the season", lo = bump(j$delta, log(0.6)), hi = bump(j$delta, log(1.6)),
         lab = c("-40%", "+60%"),
         note = "Also scales the curve -- but shared across countries, so it cannot absorb c."),
    list(key = "off_eld  age reporting, 65+", lo = bump(j$off_eld, -1), hi = bump(j$off_eld, 1),
         lab = c("half", "double"), grp = 3L,
         note = "Scales ONE age group's curve only (65+ shown). Surveillance, not biology.")
  )
  # colour by DIRECTION, not by the label: every panel then reads the same way -- blue is the lower
  # setting and orange the higher one, whatever the units of that particular parameter happen to be.
  rows = do.call(rbind, lapply(panels, function(p){
    g = if (is.null(p$grp)) 2L else p$grp
    if (any(is.na(unlist(p[c("lo", "hi")])))) return(NULL)
    rbind(cbind(mu_of(th, g),   dir = "as fitted", key = p$key),
          cbind(mu_of(p$lo, g), dir = "lower",     key = p$key),
          cbind(mu_of(p$hi, g), dir = "higher",    key = p$key))
  }))
  # wrap to the panel width, or the note runs off the right edge and the sentence is lost
  notes = data.frame(key = vapply(panels, `[[`, character(1), "key"),
                     note = vapply(panels, function(p)
                       paste(strwrap(sprintf("%s  (%s / %s)", p$note, p$lab[1], p$lab[2]), width = 46),
                             collapse = "\n"), character(1)),
                     stringsAsFactors = FALSE)
  rows$key = factor(rows$key, levels = notes$key)
  notes$key = factor(notes$key, levels = notes$key)
  notes$x = 1; notes$y = Inf
  rows$dir = factor(rows$dir, levels = c("lower", "as fitted", "higher"))
  pmain = ggplot(rows, aes(week, value, group = dir)) +
    geom_line(data = rows %>% filter(dir != "as fitted"), aes(colour = dir), linewidth = 0.75) +
    geom_line(data = rows %>% filter(dir == "as fitted"), colour = "grey25", linewidth = 1.3) +
    geom_text(data = notes, aes(x = x, y = y, label = note), hjust = 0, vjust = 1.15,
              size = 2.8, colour = "grey30", lineheight = 1.1, inherit.aes = FALSE) +
    facet_wrap(~ key, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = c(lower = .jm_blue, higher = .jm_orange), name = NULL,
                        labels = c(lower = "the lower setting", higher = "the higher setting")) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.55))) +
    labs(x = "week of the season (from 1 August)",
         y = sprintf("ILI+ per 100 000 (%s %s)", cc, ss)) +
    .jm_theme(9) + theme(strip.text = element_text(face = "bold", size = 8.5),
                         plot.margin = margin(4, 8, 4, 4))

  # the two exact trade-offs, each drawn as a pair that coincides
  stopifnot(!is.na(j$delta))
  t_R0S0 = bump(j$R0, log(1.25))                   # +25% R0 ...
  t_R0S0[j$S0] = qlogis(plogis(th[j$S0]) / 1.25)   # ... with S0 divided by 1.25: product held
  t_cdel = bump(j$c, log(1.5)); t_cdel[j$delta] = th[j$delta] + log(1 / 1.5)
  tr = rbind(cbind(mu_of(th),     set = "as fitted",                      pair = "R0 x S0 held constant"),
             cbind(mu_of(t_R0S0), set = "R0 +25%, S0 /1.25",              pair = "R0 x S0 held constant"),
             cbind(mu_of(th),     set = "as fitted",                      pair = "c x delta held constant"),
             cbind(mu_of(t_cdel), set = "c +50%, visibility /1.5",        pair = "c x delta held constant"))
  ptr = ggplot(tr %>% filter(value > 0.5), aes(week, value, colour = set, linetype = set)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~ pair, ncol = 2) +
    scale_y_log10() +
    scale_colour_manual(values = c("as fitted" = "grey30", "R0 +25%, S0 /1.25" = "#C1541E",
                                   "c +50%, visibility /1.5" = "#C1541E"), name = NULL) +
    scale_linetype_manual(values = c("as fitted" = "solid", "R0 +25%, S0 /1.25" = "22",
                                     "c +50%, visibility /1.5" = "22"), name = NULL) +
    labs(title = "Why the sharing in the design is not optional",
         subtitle = paste("Each panel holds a PRODUCT fixed and moves both of its factors; log scale, so the rise is a",
                          "straight line. RIGHT: reporting level and season visibility give exactly the same curve --",
                          "they are the same number to the data. LEFT: transmissibility and susceptibility give the same",
                          "RISE and separate only at the peak, because susceptibility also sets how many people there",
                          "are left to infect. So one wave pins their product and almost nothing else -- which is why",
                          "susceptibility needs transmissibility shared across countries before it means anything.", sep = "\n"),
         x = "week of the season", y = "ILI+ per 100 000 (log scale)") +
    .jm_theme(9) + theme(legend.position = "top")

  patchwork::wrap_plots(pmain, ptr, ncol = 1, heights = c(2.1, 1)) +
    patchwork::plot_annotation(
      title = "What each parameter does to a wave",
      subtitle = paste0("One country-season (", cc, " ", ss,
                        "), adults unless stated. Grey is the fit; coloured lines move ONE parameter and leave the rest alone.",
                        "\nEvery curve is produced by the model's own C++, so this figure cannot drift from the model it explains."),
      theme = theme(plot.title = element_text(face = "bold", size = 15),
                    plot.subtitle = element_text(size = 10, colour = "grey30", lineheight = 1.25)))
}

# ---- |-04 what varies where ----
# The spec is a SEPARATE function so a test can hold it against the real layout: the `n` column must
# add up to d$n_par exactly. It drifted once already -- the season-visibility row claimed one number
# per season when the sum-to-zero constraint leaves only S-1 of them free -- and a figure that
# miscounts the model is worse than no figure, because it is the one a reader trusts for the design.
jm_design_spec = function(d){
  S = d$n_season; C = d$n_country
  tibble::tribble(
    ~group,        ~parameter,                  ~season, ~country, ~age, ~n,                       ~note,
    "dynamics",    "R0  transmissibility",      TRUE,  FALSE, FALSE, S,           "one number per season for all of Europe",
    "dynamics",    "S0  susceptibility",        FALSE, TRUE,  FALSE, C,           "one per country, same across its seasons",
    "dynamics",    "sigma  elderly suscept.",   FALSE, FALSE, TRUE,  1,           "one number for everyone",
    "dynamics",    "I0  seed / arrival",        TRUE,  TRUE,  FALSE, d$n_cs,      "free for every wave: sets when it arrives",
    "observation", "c  reporting level",        FALSE, TRUE,  FALSE, C,           "one per surveillance system",
    "observation", "delta  season visibility",  TRUE,  FALSE, FALSE, S - 1L,      "one per season, the last set by the average-one constraint",
    "observation", "off  age reporting",        FALSE, TRUE,  TRUE,  2 * C,       "adults the reference",
    "observation", "b  off-season baseline",    FALSE, TRUE,  FALSE, sum(d$n_src),"one per data source present",
    "observation", "phi  dispersion",           FALSE, TRUE,  FALSE, C,           "one per country",
    "fixed",       "gamma  infectious period",  FALSE, FALSE, FALSE, 0L,          "3.6 days, from the literature",
    "fixed",       "vaccine effects (3)",       FALSE, FALSE, TRUE,  0L,          "fixed, 65+ pulse on 1 October",
    "fixed",       "contact matrix",            FALSE, TRUE,  TRUE,  0L,          "fixed, rescaled to spectral radius 1")
}

plot_jm_design = function(fit){
  d = fit$d; spec = jm_design_spec(d)
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
         caption = sprintf("The counts add to %d fitted numbers, estimated from %d country-seasons of weekly age-specific data across %d countries and %d seasons.",
                           sum(spec$n), d$n_cs, d$n_country, d$n_season),
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
  # the value label sits to the right of the whole INTERVAL, not of the point -- at the point it is
  # overprinted by the upper cap as soon as intervals are drawn. Added BEFORE ggplot(), which captures
  # its data by value: a column attached afterwards is not in the plot and the aesthetic fails.
  s$lab_at = if (!is.null(ci)) s$upper else s$S0
  g = ggplot(s, aes(S0, reorder(country, S0))) +
    annotate("rect", xmin = plogis(d$pr_S0_mean - 1.96 * d$pr_S0_sd),
             xmax = plogis(d$pr_S0_mean + 1.96 * d$pr_S0_sd), ymin = -Inf, ymax = Inf,
             fill = .jm_blue, alpha = 0.08) +
    geom_vline(xintercept = plogis(d$pr_S0_mean), linetype = "dashed", colour = "grey50")
  if (!is.null(ci)) g = g + geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.22, colour = .jm_blue)
  g + geom_point(size = 3.2, colour = .jm_blue) +
    geom_text(aes(x = lab_at, label = sprintf("%.3f", S0)), hjust = -0.25, size = 3, colour = "grey20") +
    scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
    labs(title = "What it learns, 3: how susceptible each country was at the season start",
         subtitle = paste0("S0 is the share of the population that could be infected on 1 August, held the same across that",
                          "\ncountry's seasons. It is the project's target quantity, and it is identifiable here only because",
                          "\ntransmissibility is shared across countries: fitted one country at a time, the two trade off almost",
                          "\nperfectly and the answer comes from the prior. That cuts both ways, and it is the caveat to carry:",
                          "\nif countries genuinely differ in transmissibility, the difference has nowhere to go but into S0. On",
                          "\nsimulated data where countries' R0 really did differ by 10%, this RANKING fell from 0.97 to 0.09.",
                          "\nSo read it as conditional on shared transmissibility, and trust the ranking over the absolute level.",
                          .jm_ivnote(iv)),
         x = "S0", y = NULL) + .jm_theme()
}

# ---- |-09 reporting level ----
plot_jm_country_reporting = function(fit, iv = NULL){
  d = fit$d; s = jm_summary_country(fit)
  ci = .jm_iv(iv, paste0(d$countries, ":log_c"))
  if (!is.null(ci)){ s$lower = 100 * ci$lower; s$upper = 100 * ci$upper }
  s$lab_at = if (!is.null(ci)) s$upper else s$c_adult * 100      # right of the interval, not the point
  g = ggplot(s, aes(c_adult * 100, reorder(country, c_adult)))
  if (!is.null(ci)) g = g + geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.22, colour = .jm_green)
  g + geom_point(size = 3.2, colour = .jm_green) +
    geom_text(aes(x = lab_at, label = sprintf("%.1f%%", c_adult * 100)), hjust = -0.25, size = 3,
              colour = "grey20") +
    scale_x_log10(expand = expansion(mult = c(0.08, 0.30))) +
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
  # Build the parameter name FROM EACH ROW. pivot_longer interleaves (country1 young, country1
  # elderly, country2 young, ...) while a blocked request returns (all young, then all elderly), so
  # binding the two attached every bar to the wrong country and group. Constructing the name per row
  # cannot go out of order however the frame is reshaped; the test asserts each bar brackets its point.
  age$par = paste0(age$country, ifelse(age$group == "young", ":off_young", ":off_eld"))
  ci = .jm_iv(iv, age$par)
  if (!is.null(ci)){ age$lower = ci$lower; age$upper = ci$upper }
  g = ggplot(age, aes(rel, reorder(country, rel), colour = group)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50")
  if (!is.null(ci)) g = g + geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0, linewidth = 0.5,
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
         subtitle = paste0("Modelled share of each age group infected over a full ", d$attack_weeks,
                          "-week season. Thin lines are countries, thick lines",
                          "\nthe median across them. This is the one output entirely free of reporting: it comes from the contact",
                          "\nmatrix, the global elderly susceptibility (fitted at ", sprintf("%.2f", p$sigma_eld),
                          "x an adult's) and vaccination, not from how",
                          "\nmany consultations were counted. So compare it with prospective cohort evidence, never with",
                          "\nsurveillance curves. The dynamics run to the same horizon for every country-season, so these are",
                          "\ncomparable with each other: read over only the weeks each country happened to report, five cells",
                          "\nmoved by more than 5% and one by 15% for no epidemiological reason at all."),
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
  # run_joint_recovery.R saves THREE studies in one object (base / local / prior / driver). Being
  # handed that wrapper instead of one study's `rec` is the easy mistake; say so instead of failing
  # inside a dplyr verb with an unrecognisable message.
  if (!is.data.frame(rec$comparison))
    stop("plot_jm_recovery() needs one study, not the whole saved object: pass e.g. ",
         "readRDS('output/joint_model/joint_recovery.rds')$local$rec")
  if (is.null(summ)) summ = jm_recovery_summary(rec, d)
  cmp = summ$comparison %>%
    filter(family %in% c("R0 (season, shared)", "season deviation (shared)", "S0 (country)",
                         "reporting c (country)", "elderly susceptibility (global)"))
  cov_txt = summ$by_family %>% filter(!is.na(coverage)) %>%
    summarise(m = median(coverage)) %>% pull(m)
  # A replicate in which one country landed in a bad local optimum throws an estimate far off scale
  # and would squash every other panel flat. Such points are WINSORISED FOR DISPLAY ONLY, drawn as
  # open triangles at the panel edge and counted in the subtitle, so they are visible rather than
  # hidden and the informative range stays readable. The cap is 35% of the family's own truth range
  # beyond its extremes, falling back to 35% of the level itself where a family has a single truth
  # value (the global elderly susceptibility), so that band is never zero-width.
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
                          ". Orange triangles are ", n_off, " estimate(s) that fell outside",
                          "\ntheir panel's range; they are drawn at the edge rather than dropped, so a failure cannot hide.",
                          "\nAnything the fit cannot recover from its own simulation cannot be trusted from real data either."),
         x = "true value", y = "estimated value") + .jm_theme()
}

# ================= write the default set =================
save_jm_report = function(fit, id = NULL, iv = NULL, rec = NULL, dir = "output/joint_model"){
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  # Omitting `iv` is legitimate -- a quick look without waiting minutes for the Hessian -- but it must
  # be announced, because figures 09-13 then carry NO uncertainty and look no different for it. The
  # pipeline wrote the entire set that way for two days without anything saying so.
  if (is.null(iv))
    warning("no intervals supplied: figures 09-13 will be written WITHOUT uncertainty bars. ",
            "Pass iv = jm_intervals(fit) for the publishable set.", call. = FALSE)
  d = fit$d; man = character(0)
  put = function(file, plot, w, h, what){
    ggsave(file.path(dir, file), plot, width = w, height = h, dpi = 115, limitsize = FALSE)
    man <<- c(man, sprintf("| `%s` | %s |", file, what))
  }
  put("01_data_panel.png", plot_jm_data_panel(fit), 10.5, 6.6,
      "what data there is: every country-season, its source, its observed weeks, and what was excluded")
  put("02_data_features.png", plot_jm_data_features(fit), 14, 10.5,
      "the four features of the data that dictate the model's design -- read before the model")
  put("03_mechanism.png", plot_jm_mechanism(fit), 13, 11,
      "how the model works: what each parameter does to a wave, and the two products the design has to break")
  put("04_design_what_varies_where.png", plot_jm_design(fit), 12.5, 6.2,
      "which parameters vary by season, by country, by age, and which are fixed")
  put("05_arrival_time.png", plot_jm_arrival(fit), 10, 5.8,
      "how the model sets when each wave arrives (the seed)")
  put("06_fit_overview.png", plot_jm_fit_overview(fit), 2.05 * d$n_season + 2, 1.35 * d$n_country + 2.8,
      "what data it fits: every country-season, observed against modelled")
  for (cc in d$countries)
    ggsave(file.path(dir, sprintf("07_fit_%s.png", cc)), plot_jm_fit_country(fit, cc),
           width = 2.05 * d$n_season + 2, height = 7.2, dpi = 110, limitsize = FALSE)
  man = c(man, sprintf("| `07_fit_<country>.png` | the same by age group, one file per country (%s) |",
                       paste(d$countries, collapse = ", ")))
  put("08_noise_budget.png", plot_jm_adequacy(fit), 10, 5.8,
      "how well it fits, honestly: the noise the fit needed against the noise the data have")
  put("09_season_R0.png", plot_jm_season_R0(fit, iv), 10, 6.0,
      "what it learns: transmissibility of each season's virus, shared across countries")
  put("10_season_visibility.png", plot_jm_season_visibility(fit, iv), 10, 6.0,
      "what it learns: how visible each season was per infection, shared across countries")
  put("11_country_S0.png", plot_jm_country_S0(fit, iv), 10, 6.4,
      "what it learns: susceptibility at the season start, by country")
  put("12_country_reporting.png", plot_jm_country_reporting(fit, iv), 10, 5.8,
      "what it learns: fraction of adult infections that is counted, by country")
  put("13_age_reporting.png", plot_jm_age_offsets(fit, iv), 10, 5.8,
      "what it learns: how visible children and the elderly are per infection, by country")
  put("14_attack_rates.png", plot_jm_attack(fit), 10, 6.0,
      "what it learns: modelled attack rate by age group, the one reporting-free output")
  if (!is.null(id)) put("15_data_or_prior.png", plot_jm_identifiability(id), 10, 5.4,
                        "whether to believe it: how much each parameter owes to the data rather than its prior")
  # Figure 16 is normally written by run_joint_recovery.R, which is a separate (much longer) run. So
  # when this function is called without `rec` an EXISTING 16 is left on disk untouched: list it
  # anyway, with its date, rather than leaving a figure present but unmentioned and silently older
  # than the rest of the set.
  f16 = file.path(dir, "16_recovery.png")
  if (!is.null(rec)){
    put("16_recovery.png", plot_jm_recovery(rec, d), 10, 7.0,
        "whether to believe it: recovery of a known truth simulated from the model onto the real design")
  } else if (file.exists(f16)){
    man = c(man, sprintf(paste0("| `16_recovery.png` | whether to believe it: recovery of a known truth. ",
                                "NOT regenerated by this run -- written by `run_joint_recovery.R`, last on %s |"),
                         format(file.mtime(f16), "%Y-%m-%d %H:%M")))
  } else {
    man = c(man, "| `16_recovery.png` | whether to believe it: recovery of a known truth. NOT YET RUN -- produce it with `Rscript code/07_joint_model/run_joint_recovery.R` |")
  }
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
