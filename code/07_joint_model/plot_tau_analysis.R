# plot_tau_analysis.R -- figure 18: what makes a wave fat, and what the data can tell apart.
# Called by run_tau_analysis.R with the fits, profiles, recovery and comparisons it computed.
suppressMessages({library(ggplot2); library(patchwork); library(dplyr)})

plot_tau_analysis = function(f0, fS, fC, prof, prof_c, rec, rec_more = NULL, cmp, theory){
  blue = "#2C6E9B"; orange = "#D0671F"; green = "#5E8C31"; purple = "#8E5BA8"; grey = "grey35"
  th_small = theme_minimal(10) + theme(plot.title = element_text(face = "bold", size = 11),
                                       plot.subtitle = element_text(size = 8.6, colour = "grey30", lineheight = 1.15),
                                       legend.position = "top", legend.text = element_text(size = 8))
  d = fS$d; th = fS$theta; S = d$n_season; per_of = function(i) d$rate_per / sum(d$N[[d$cs_country[i] + 1L]])

  # (a) R0 and S0 are ONE number to the data: (k R0, S0/k, I0/k, k c) reproduces every expected count
  i_ref = which(d$cs_label == "PL 2017/2018"); if (!length(i_ref)) i_ref = 1L
  ic = d$cs_country[i_ref] + 1L; base = d$off_country[ic]; s = d$cs_season[i_ref] + 1L
  xs = c(th[seq_len(S - 1)], -sum(th[seq_len(S - 1)])); S0 = plogis(th[base + 1] + xs[s])
  jI = base + 5 + d$n_src[ic] + d$cs_pos[i_ref] + 1L; k = 1.2
  t_ex = th; t_ex[base + 1] = qlogis(S0 / k) - xs[s]; t_ex[base + 2] = th[base + 2] + log(k); t_ex[jI] = th[jI] - log(k)
  d_k = d; d_k$R0_fixed = d$R0_fixed * k
  curve = function(theta, dd) rowSums(jm_fitted_cpp(theta, dd)$mu[[i_ref]]) * per_of(i_ref)
  m_fit = curve(th, d); m_ex = curve(t_ex, d_k); m_raw = curve(th, d_k)
  dev = max(vapply(seq_len(d$n_cs), function(i){       # the same exchange in every wave, for the caption
    ic2 = d$cs_country[i] + 1L; b2 = d$off_country[ic2]; s2 = d$cs_season[i] + 1L; S02 = plogis(th[b2 + 1] + xs[s2])
    if (S02 / k >= 1) return(0)
    t2 = th; t2[b2 + 1] = qlogis(S02 / k) - xs[s2]; t2[b2 + 2] = th[b2 + 2] + log(k)
    j2 = b2 + 5 + d$n_src[ic2] + d$cs_pos[i] + 1L; t2[j2] = th[j2] - log(k)
    max(abs(jm_fitted_cpp(t2, d_k)$mu[[i]] / jm_fitted_cpp(th, d)$mu[[i]] - 1)) }, numeric(1)))
  da = rbind(data.frame(week = seq_along(m_fit), v = m_fit, set = sprintf("R0 %.1f, S0 %.2f (as fitted)", d$R0_fixed, S0)),
             data.frame(week = seq_along(m_ex), v = m_ex, set = sprintf("R0 %.1f, S0 %.2f, seed and c rescaled", d$R0_fixed * k, S0 / k)),
             data.frame(week = seq_along(m_raw), v = m_raw, set = sprintf("R0 %.1f, S0 %.2f (only R0 moved)", d$R0_fixed * k, S0)))
  da$set = factor(da$set, levels = unique(da$set))
  pa = ggplot(da %>% filter(v > 0.3), aes(week, v, colour = set, linetype = set)) + geom_line(linewidth = 1) +
    scale_y_log10() + scale_colour_manual(values = c(grey, orange, blue), name = NULL) +
    scale_linetype_manual(values = c("solid", "22", "solid"), name = NULL) +
    guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
    labs(title = "(a) R0 and S0 are one number to the data",
         subtitle = sprintf(paste("A wave's shape depends on R0 x S0 only; S0 also scales its height, which reporting",
                                  "absorbs. Moving R0 up by 20%% with S0, the seed and c rescaled gives the SAME curve",
                                  "(dashed on grey): max relative difference %.0e over all %d waves. Fixing R0 is",
                                  "therefore a judgement; its value reaches the data only through S0 <= 1.", sep = "\n"), dev, d$n_cs),
         x = sprintf("week (%s)", d$cs_label[i_ref]), y = "ILI+ per 100 000 (log)") + th_small

  # (b) S0 and tau in theory: the rise-width plane of the model's own noiseless waves
  ss = theory %>% filter(is.finite(rise), is.finite(fwhm))
  lab_pts = ss %>% filter(tau == 0, abs(S0 - round(S0 / 0.05) * 0.05) < 1e-9, S0 >= 0.75)
  pb = ggplot(ss, aes(rise, fwhm, colour = factor(tau))) + geom_path(linewidth = 1) +
    geom_point(data = lab_pts, colour = grey, size = 1.6) +
    geom_text(data = lab_pts, aes(label = sprintf("S0 %.2f", S0)), colour = grey, size = 2.6, hjust = -0.15, vjust = -0.4) +
    scale_colour_manual(values = c(`0` = grey, `7` = green, `14` = orange, `21` = purple), name = "tau (days)") +
    coord_cartesian(ylim = c(4.5, 17)) +
    labs(title = "(b) S0 moves a wave along the line, tau lifts it off",
         subtitle = paste("Each line is the model's own wave (one design, noise-free) as S0 varies at a fixed spread.",
                          "Lower S0: slower rise AND wider. More spread: wider at almost the same rise. Up to a",
                          "week of spread stays within half a week of the no-spread line -- a small S0 change",
                          "mimics it; from two weeks on the wave leaves the line and tau becomes visible.", sep = "\n"),
         x = "rise rate (per week, from 10% to 60% of the peak)", y = "width at half maximum (weeks)") + th_small

  # (c) one tau for all: the profile
  pr = prof %>% arrange(tau) %>% mutate(dll = loglik - max(loglik))
  tau_hat = exp(fS$theta[["log_tau"]])
  at = function(v) pr$dll[which.min(abs(pr$tau - v))]     # the held values are near, not on, the grid
  pc = ggplot(pr, aes(tau, dll)) +
    geom_hline(yintercept = -qchisq(0.95, 1) / 2, linetype = "22", colour = "grey55") +
    geom_line(colour = blue, linewidth = 1) + geom_point(colour = blue, size = 2) +
    annotate("point", x = tau_hat, y = 0, colour = orange, size = 3) +
    annotate("text", x = tau_hat, y = 0, label = sprintf("  fitted %.1f d", tau_hat), hjust = 0, vjust = -0.6, size = 3, colour = orange) +
    labs(title = "(c) One tau for every country: a few days at most",
         subtitle = sprintf(paste("Profile: tau held, the other 181 parameters refitted. Flat up to about two days, then",
                                  "falling: %.1f nats lost at one week, %.0f at two, %.0f at three; dotted: the 95%% cut.",
                                  "The fitted 2.1 days gains %+.2f nats over none. Not the S0 ceiling: with R0 pinned",
                                  "at 1.7 the shared tau again settles near two days.", sep = "\n"),
                            -at(7), -at(14), -at(21), fS$loglik - f0$loglik),
         x = "tau, days (held)", y = "log-likelihood vs the best") + th_small

  # (d) one tau per country: estimate, how sharply its own data pin it, and country size
  pC = jm_unpack(fC$theta, fC$d)
  pcs = prof_c %>% group_by(country) %>% mutate(dll = loglik - max(loglik)) %>% ungroup()
  gain0 = pcs %>% group_by(country) %>% summarise(lose_at_0 = -dll[tau == 0], .groups = "drop")
  dd = data.frame(country = fC$d$countries, pop = vapply(fC$d$N, sum, 0) / 1e6,
                  tau = vapply(pC$country, `[[`, 0, "tau")) %>% left_join(gain0, by = "country")
  k2 = cmp[cmp$contrast == "tau by country vs shared tau", ]
  pd = ggplot(dd, aes(pop, tau)) + geom_point(aes(size = lose_at_0), colour = orange, alpha = 0.8) +
    geom_text(aes(label = country), size = 3, vjust = -1.1) + scale_x_log10() +
    scale_y_log10(expand = expansion(mult = c(0.06, 0.14))) +
    scale_size_area(max_size = 9, name = "nats lost if\nforced to 0") +
    labs(title = "(d) One tau per country: real, but not geography",
         subtitle = sprintf(paste("+%.0f nats over one shared tau, but %d of %d waves agree (sign test p = %.2f) and Estonia",
                                  "carries half of it. Country size does not predict tau (rank correlation %.2f): Estonia's",
                                  "36 days lays a broad mean through a very noisy series. tau is a blunt sensor of how fat a",
                                  "country's waves are -- spread, subtype mixing, reporting -- not a map of its cities.", sep = "\n"),
                            k2$total, k2$favour_B, k2$n_waves, k2$sign_p, cor(dd$tau, dd$pop, method = "spearman")),
         x = "population (millions, log)", y = "fitted tau (days, log)") + th_small + theme(legend.position = "right")

  # (e) the trade-off, country by country: hold tau, refit the country, watch its S0
  pcs = pcs %>% filter(tau <= 30)
  best = pcs %>% group_by(country) %>% slice_max(loglik, n = 1) %>% ungroup()
  spread_wanted = best$country[best$tau >= 7]
  pe = ggplot(pcs, aes(tau, S0, group = country)) +
    geom_line(aes(colour = country %in% spread_wanted), show.legend = FALSE) +
    geom_point(aes(shape = dll > -qchisq(0.95, 1) / 2), size = 1.5, colour = "grey45") +
    geom_point(data = best, colour = orange, size = 2.8) +
    geom_text(data = best %>% filter(country %in% spread_wanted) %>% group_by(tau) %>% arrange(S0) %>%
                mutate(vj = ifelse(row_number() %% 2 == 1, 1.9, -0.9)) %>% ungroup(),
              aes(label = country, vjust = vj), size = 3, colour = orange) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "within 1.92 nats of that country's best",
                                                                    `FALSE` = "outside"), name = NULL) +
    scale_colour_manual(values = c(`TRUE` = blue, `FALSE` = "grey75")) +
    scale_x_continuous(trans = "sqrt", breaks = c(0, 2, 7, 14, 21, 30)) +
    labs(title = "(e) More spread asks for a sharper local wave: S0 goes up",
         subtitle = sprintf(paste("One line per country: its tau held, its own parameters refitted, its S0 read off; orange",
                                  "dot = its best tau. Blue: the %d countries whose data want a spread (labelled). To keep",
                                  "the rise while the spread widens the wave, the local epidemic must be faster, so S0",
                                  "climbs with tau -- and the country ranking of S0 moves (0.80 against no spread).", sep = "\n"),
                            length(spread_wanted)),
         x = "tau held (days, square-root scale)", y = "country S0 level") + th_small

  # (f) recovery: is a known tau recovered from data like ours?
  rows = list()
  for (nm in grep("^shared_", names(rec), value = TRUE)){
    cmpr = rec[[nm]]$rec$comparison; r = cmpr[cmpr$parameter == "log_tau", ]
    rows[[nm]] = data.frame(truth = exp(r$truth_theta), est = exp(r$est_theta), what = "one shared tau")
  }
  if (!is.null(rec$country_hetero)){
    fc = rec$country_hetero$fit_c; dcc = fc$d
    tt = ifelse(rec$country_hetero$big, 14, 3)
    rows$c = data.frame(truth = tt, est = vapply(jm_unpack(fc$theta, dcc)$country, `[[`, 0, "tau"),
                        what = "per country (big 14 d, small 3 d)")
  }
  if (!is.null(rec_more))
    rows$more = data.frame(truth = rec_more$tau_true, est = rec_more$tau_est, what = "one shared tau")
  rr = do.call(rbind, rows)
  sh_rows = rr[rr$what == "one shared tau", ]; sh_rows$truth = round(sh_rows$truth, 1)
  pf = ggplot(rr, aes(truth, est, colour = what)) + geom_abline(linetype = "22", colour = "grey55") +
    geom_jitter(width = 0.03, height = 0, size = 2.4, alpha = 0.85) + scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c(blue, orange), name = NULL) +
    labs(title = "(f) Recovery: simulate with a known tau, refit",
         subtitle = sprintf(paste("Simulated from the fitted model on the real design, tau set, refitted.",
                                  "A shared week comes back at %s days, two weeks at %s:",
                                  "part is read as lower S0 -- the likelihood's own ridge, not the",
                                  "priors. Per country (orange) big and small separate but shrink",
                                  "together, and the S0 ranking drops to 0.69.", sep = "\n"),
                            paste(sprintf("%.1f", sort(sh_rows$est[sh_rows$truth == 7])), collapse = "/"),
                            paste(sprintf("%.1f", sort(sh_rows$est[sh_rows$truth == 14])), collapse = "/")),
         x = "true tau (days, log)", y = "estimated tau (days, log)") + th_small + guides(colour = guide_legend(ncol = 1))

  (pa | pb | pc) / (pd | pe | pf) +
    plot_annotation(title = "What makes a wave fat -- R0, S0 or spread -- and what the data can tell apart",
                    subtitle = paste("R0 and S0 are exactly one number (a). S0 and the spatial spread tau are separable in principle, through the two tails of a wave (b).",
                                     "The data reject a common spread (c). A per-country spread fits better but behaves as a blunt fatness sensor, not as geography (d),",
                                     "and trades against S0 within each country (e); on simulated data, spreads up to a week are partly read as lower S0 (f).", sep = "\n"),
                    theme = theme(plot.title = element_text(face = "bold", size = 15),
                                  plot.subtitle = element_text(size = 10, colour = "grey30", lineheight = 1.2)))
}
