# comp_model_report.R -- the figures to EYEBALL a compartmental-model fit (ggplot2)
#
# Three questions, three figure families (ASSUMPTIONS.md sections 6-9 give the vocabulary):
#   1. How good are the fits?   plot_cm_fit(fit)        age x season grid: observed ILI+ (points),
#      the EKF one-step-ahead filtered mean (solid; what the likelihood scores), the deterministic
#      SIR at the fitted parameters (dashed; the mechanistic mean), the vaccination pulse (dotted).
#      Rates per 100 000 of the age group so panels are comparable. Each panel carries its cor.
#      plot_cm_wiggle(fit)     how much the filter had to move the latent infected fraction away
#      from the deterministic path each week (the process noise the fit actually used).
#   2. Parameters vs priors?     plot_cm_params(fit)     every fitted parameter as a point with its
#      Laplace 95% interval on the fitted (transformed) scale, over the prior 95% band; the
#      season R0_s sit in a row so the pull-to-1.5 prior is read against the estimates at a glance.
#   3. Season and country variation?  plot_cm_seasons(fits)  across countries: R0_s by season (one
#      line per country, the shared-season hypothesis = lines move together), country S0_c with
#      intervals, attack rates by age x season (tiles), R_eff = R0_s * S0_c.
# Each function returns a ggplot (or patchwork) object; save_cm_report() writes them to output/comp_model/.

suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})

.cm_long = function(fit){
  K = length(fit$seasons); A = length(fit$groups)
  do.call(rbind, lapply(seq_len(K), function(s){
    per100k = fit$settings$rate_per / fit$N                     # counts -> rate per 100 000 of the group
    nw = nrow(fit$y[[s]]); row = function(what, m) data.frame(season = fit$seasons[s], week = seq_len(nw), group = rep(fit$groups, each = nw),
                                                              what = what, value = as.numeric(sweep(m, 2, per100k, "*")))
    rbind(row("observed", fit$y[[s]]),
          row("EKF filtered (one-step-ahead)", fit$mu_filt[[s]]),
          row("deterministic SIR", fit$mu_det[[s]]),
          if (!is.null(fit$stage1)) row("stage-1 deterministic fit", fit$stage1$mu[[s]]))
  })) %>% mutate(group = factor(group, levels = fit$groups))
}

# ---- |-1. fit quality: age x season grid ----
plot_cm_fit = function(fit){
  d = .cm_long(fit)
  cors = summarise_comp_fit(fit) %>% transmute(season, label = sprintf("R0=%.2f  cor=%.2f", R0, cor))
  vax_week = ceiling(62 / 7)
  ggplot() +
    geom_vline(xintercept = vax_week, linetype = "dotted", colour = "grey60") +
    geom_line(data = d %>% filter(what != "observed"), aes(week, value, colour = what, linetype = what), linewidth = 0.7) +
    geom_point(data = d %>% filter(what == "observed"), aes(week, value), size = 0.9, colour = "grey25", alpha = 0.8) +
    geom_text(data = cors, aes(x = 1, y = Inf, label = label), hjust = 0, vjust = 1.4, size = 2.7, colour = "grey30") +
    facet_grid(group ~ season, scales = "free_y") +
    scale_colour_manual(values = c("EKF filtered (one-step-ahead)" = "#d95f02", "deterministic SIR" = "#1b9e77", "stage-1 deterministic fit" = "#7570b3"), name = NULL) +
    scale_linetype_manual(values = c("EKF filtered (one-step-ahead)" = "solid", "deterministic SIR" = "dashed", "stage-1 deterministic fit" = "dotted"), name = NULL) +
    labs(title = sprintf("%s -- compartmental model fit: S0 = %.2f (shared across seasons), c = %s, b = %.1f, phi = %.1f, q = %.3f",
                         fit$country, fit$params$S0, paste(signif(fit$params$c_age, 3), collapse = "/"), mean(fit$params$b), fit$params$phi, fit$params$q),
         subtitle = "points = observed ILI+ per 100 000 of the age group | orange = EKF one-step-ahead mean (what the likelihood scores) | green dashed = deterministic SIR at the fitted parameters | dotted = 65+ vaccination pulse (1 Oct)",
         x = "season week (1 = week of 1 Aug)", y = "ILI+ per 100 000") +
    theme_minimal(base_size = 10) + theme(legend.position = "top", plot.subtitle = element_text(size = 8))
}

# ---- |-1b. the wiggle: where the filter left the deterministic path ----
# log2 ratio of the EKF one-step-ahead mean to the deterministic SIR mean, per week, age and season:
# 0 = the mechanistic model alone explains the week; +/-1 = the filter had to double / halve it.
# A season that is only fitted THROUGH the noise (large sustained ratios) is a season the SIR with
# the fitted S0_c and R0_s does not capture -- exactly what the process noise is meant to expose.
plot_cm_wiggle = function(fit){
  d = .cm_long(fit) %>% filter(what != "observed") %>%
    pivot_wider(names_from = what, values_from = value) %>%
    mutate(ratio = log2(pmax(`EKF filtered (one-step-ahead)`, 1e-9) / pmax(`deterministic SIR`, 1e-9)))
  ggplot(d, aes(week, group, fill = pmin(pmax(ratio, -2), 2))) + geom_tile() + facet_wrap(~ season, ncol = 4) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "log2(filtered /\ndeterministic)") +
    labs(title = sprintf("%s -- how much process noise each week needed (q = %.3f%s)", fit$country, fit$params$q, if (!is.null(fit$settings$q_fixed)) ", fixed" else ""),
         subtitle = "red = the filter ran above the deterministic SIR, blue = below; sustained colour = a season the mechanistic model alone does not capture",
         x = "season week", y = NULL) + theme_minimal(base_size = 10) + theme(plot.subtitle = element_text(size = 8))
}

# ---- |-1c. age misfit: observed / modelled cumulative ILI+ per age group and season ----
# Under age-invariant reporting (F1/F2) the model's age profile comes from the contacts alone; this
# tile shows, per season, how far each age group's observed total sits from the deterministic model's
# total. A stable pattern across seasons (e.g. elderly always ~2x) is an age-reporting or
# age-susceptibility signal the current assumptions cannot absorb.
plot_cm_age_misfit = function(fit){
  d = do.call(rbind, lapply(seq_along(fit$seasons), function(s){
    y = fit$y[[s]]; m = fit$mu_det[[s]]; ok = is.finite(y)
    data.frame(season = fit$seasons[s], group = fit$groups,
               ratio = vapply(seq_len(ncol(y)), function(a) sum(y[ok[, a], a]) / sum(m[ok[, a], a]), numeric(1)))
  })) %>% mutate(group = factor(group, levels = fit$groups))
  ggplot(d, aes(season, group, fill = pmin(pmax(log2(ratio), -2), 2))) + geom_tile() +
    geom_text(aes(label = sprintf("%.2fx", ratio)), size = 3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "log2(obs /\nmodel)") +
    labs(title = sprintf("%s -- observed / deterministic-model cumulative ILI+ by age group and season", fit$country),
         subtitle = "1x = the age profile implied by the contact structure and shared reporting matches; a stable colour per row = a systematic age effect",
         x = NULL, y = NULL) + theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.subtitle = element_text(size = 8))
}

# ---- |-2. parameters vs priors (transformed scale, Laplace 95% intervals) ----
plot_cm_params = function(fit){
  s = fit$settings; th = fit$theta; se = fit$se; nm = names(th)
  prior_mean = ifelse(grepl("^log_R0_", nm), log(s$R0_reference), NA_real_)
  prior_sd   = ifelse(grepl("^log_R0_", nm), s$prior_logR0_sd, NA_real_)
  prior_mean[nm == "logit_S0"] = s$prior_logitS0["mean"]; prior_sd[nm == "logit_S0"] = s$prior_logitS0["sd"]
  prior_mean[nm == "log_phi"]  = s$prior_logphi["mean"];  prior_sd[nm == "log_phi"]  = s$prior_logphi["sd"]
  prior_mean[nm == "log_q"]    = s$prior_logq["mean"];    prior_sd[nm == "log_q"]    = s$prior_logq["sd"]
  prior_mean[grepl("^log_I0_", nm)] = log(s$I0_fraction); prior_sd[grepl("^log_I0_", nm)] = s$prior_logI0_sd
  prior_mean[grepl("^logc_dev_", nm)] = 0; prior_sd[grepl("^logc_dev_", nm)] = s$prior_logc_season_sd
  prior_mean[grepl("^log2c_", nm)] = 0;    prior_sd[grepl("^log2c_", nm)] = s$prior_logc_age_sd
  d = data.frame(param = nm, est = th, lo = th - 1.96 * se, hi = th + 1.96 * se, prior_mean, prior_sd, stringsAsFactors = FALSE) %>%
    mutate(family = case_when(grepl("^log_R0_", param) ~ "season R0 (log)", grepl("^log_I0_", param) ~ "season seed I0 (log; arrival time)",
                              param == "logit_S0" ~ "S0 (logit)",
                              grepl("^logc_dev_", param) ~ "season reporting deviation (log)",
                              param %in% c("log_c", "log_b") | grepl("^log2c_|^log_b_", param) ~ "reporting (log; age offsets log2)", TRUE ~ "noise / dispersion (log)"),
           label = sub("^log_R0_", "", param), label = sub("^logit_", "", label), label = sub("^log_", "", label),
           label = factor(label, levels = rev(unique(label))))
  ggplot(d, aes(y = label)) +
    geom_linerange(data = d %>% filter(is.finite(prior_sd)), aes(y = label, xmin = prior_mean - 1.96 * prior_sd, xmax = prior_mean + 1.96 * prior_sd),
                   colour = "#8da0cb", linewidth = 6, alpha = 0.35) +
    geom_vline(data = d %>% filter(is.finite(prior_mean)) %>% distinct(family, prior_mean), aes(xintercept = prior_mean),
               colour = "#8da0cb", linetype = "dashed") +
    geom_pointrange(aes(x = est, xmin = lo, xmax = hi), colour = "#d95f02") +
    facet_wrap(~ family, scales = "free", ncol = 3) +
    labs(title = sprintf("%s -- fitted parameters (points, Laplace 95%% intervals) against their priors (blue band = prior 95%%)", fit$country),
         subtitle = "on the fitted (log / logit) scale; c and b carry no prior. Season R0 read against the log(1.5) +/- 0.10 pull.",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) + theme(plot.subtitle = element_text(size = 8))
}

# ---- |-3. across countries: season R0, country S0, attack rates ----
plot_cm_seasons = function(fits){
  summ = do.call(rbind, lapply(fits, summarise_comp_fit))
  se_S0 = vapply(fits, function(f) f$se["logit_S0"], numeric(1))
  S0 = data.frame(country = vapply(fits, `[[`, "", "country"), S0 = vapply(fits, function(f) f$params$S0, numeric(1)),
                  lo = plogis(vapply(fits, function(f) f$theta["logit_S0"], numeric(1)) - 1.96 * se_S0),
                  hi = plogis(vapply(fits, function(f) f$theta["logit_S0"], numeric(1)) + 1.96 * se_S0))
  p1 = ggplot(summ, aes(season, R0, group = country, colour = country)) +
    geom_hline(yintercept = 1.5, colour = "grey60", linetype = "dashed") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5 * exp(-1.96 * fits[[1]]$settings$prior_logR0_sd),
             ymax = 1.5 * exp(1.96 * fits[[1]]$settings$prior_logR0_sd), fill = "#8da0cb", alpha = 0.25) +
    geom_line() + geom_point(size = 2) +
    labs(title = "Season transmissibility R0_s per country (per-country stage; the joint stage will share one value per season)",
         subtitle = "blue band = prior 95% around 1.5. Lines moving together across countries = evidence for a shared season factor.",
         x = NULL, y = "fitted R0_s") + theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.subtitle = element_text(size = 8))
  p2 = ggplot(S0, aes(S0, reorder(country, S0))) + geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#d95f02") +
    labs(title = "Country susceptibility S0_c (shared across seasons)", subtitle = "Laplace 95% intervals; relative comparison only", x = "S0", y = NULL) +
    theme_minimal(base_size = 10) + theme(plot.subtitle = element_text(size = 8))
  att = summ %>% select(country, season, young = attack_young, medium = attack_medium, elderly = attack_elderly) %>%
    pivot_longer(c(young, medium, elderly), names_to = "group", values_to = "attack") %>%
    mutate(group = factor(group, levels = c("young", "medium", "elderly")))
  p3 = ggplot(att, aes(season, group, fill = attack)) + geom_tile() + facet_wrap(~ country, ncol = 3) +
    geom_text(aes(label = sprintf("%.0f%%", 100 * attack)), size = 2.5) +
    scale_fill_viridis_c(name = "attack\nrate", option = "C", limits = c(0, NA)) +
    labs(title = "Modelled attack rate (cumulative infections / group) by age group and season", x = NULL, y = NULL) +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p4 = ggplot(summ, aes(season, R_eff, group = country, colour = country)) + geom_hline(yintercept = 1, colour = "grey60") +
    geom_line() + geom_point(size = 2) +
    labs(title = "Effective reproduction number at season start, R0_s x S0_c", x = NULL, y = "R_eff") +
    theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  list(R0 = p1, S0 = p2, attack = p3, R_eff = p4)
}

# ---- |-write the per-country and cross-country figures ----
save_cm_report = function(fits, dir = "output/comp_model"){
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (fit in fits){
    K = length(fit$seasons)
    ggsave(file.path(dir, sprintf("fit_%s.png", fit$country)), plot_cm_fit(fit), width = 2.6 * K + 2, height = 7.5, dpi = 110)
    ggsave(file.path(dir, sprintf("params_%s.png", fit$country)), plot_cm_params(fit), width = 13, height = 5.5, dpi = 110)
    ggsave(file.path(dir, sprintf("wiggle_%s.png", fit$country)), plot_cm_wiggle(fit), width = 11, height = 4.5, dpi = 110)
    ggsave(file.path(dir, sprintf("agemisfit_%s.png", fit$country)), plot_cm_age_misfit(fit), width = 8, height = 3.2, dpi = 110)
  }
  if (length(fits) > 1){
    ps = plot_cm_seasons(fits)
    ggsave(file.path(dir, "seasons_R0.png"), ps$R0, width = 9, height = 5, dpi = 110)
    ggsave(file.path(dir, "countries_S0.png"), ps$S0, width = 6, height = 4, dpi = 110)
    ggsave(file.path(dir, "attack_rates.png"), ps$attack, width = 11, height = 2.2 * ceiling(length(fits) / 3) + 1.5, dpi = 110)
    ggsave(file.path(dir, "R_eff.png"), ps$R_eff, width = 9, height = 5, dpi = 110)
  }
  write.csv(do.call(rbind, lapply(fits, summarise_comp_fit)), file.path(dir, "comp_model_summary.csv"), row.names = FALSE)
  invisible(dir)
}
