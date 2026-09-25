# plot_model_comparison.R -- ONE figure for "where does the variation live", read from the outputs of
# run_model_comparison.R. Three panels:
#   (a) the 2x2 of models, each cell its log-likelihood relative to the working model;
#   (b) per country-season, which mechanism the wave prefers for the SEASON effect (R0 minus S0, with
#       the country effect held on S0 in both) -- a diverging heatmap, so the answer is visible wave
#       by wave rather than only in the total;
#   (c) the wave that discriminates most, with both fitted curves over the data, so the shape
#       difference the likelihood is scoring can be seen.
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

plot_model_comparison = function(cmp, fits, out = NULL){
  tab = cmp$table; cc = cmp$pairwise
  blue = "#2B5D8A"; orange = "#C1541E"

  # (a) the 2x2
  g = tab %>% filter(model %in% c("S0/S0", "R0/R0", "R0/S0", "S0/R0")) %>%
    mutate(dll = loglik - loglik[model == "S0/S0"],
           season_on = factor(paste("season effect on", season_on), levels = paste("season effect on", c("S0", "R0"))),
           country_on = factor(paste("country effect on", country_on), levels = paste("country effect on", c("S0", "R0"))),
           lab = sprintf("%s\n%+.1f nats\n%s", model, dll, ifelse(converged, "converged", "NOT converged")))
  pa = ggplot(g, aes(country_on, season_on, fill = dll)) +
    geom_tile(colour = "white", linewidth = 2) +
    geom_text(aes(label = lab), size = 3.3, lineheight = 0.95,
              colour = ifelse(abs(g$dll) > 8, "white", "grey15")) +
    scale_fill_gradient2(low = blue, mid = "grey92", high = orange, midpoint = 0,
                         name = "log-likelihood\nvs the working model") +
    scale_y_discrete(limits = rev) +
    labs(title = "(a) Four models, one parameter count",
         subtitle = paste("Each quantity that carries no effect is pinned (R0 = 1.5, S0 = 0.75). Same data, same likelihood,",
                          "same count, so the numbers compare directly. Positive favours that model over the working one.", sep = "\n"),
         x = NULL, y = NULL) +
    theme_minimal(10) + theme(panel.grid = element_blank(), legend.position = "right",
                              plot.title = element_text(face = "bold"))

  # (b) the per-cell heatmap
  seasons = unique(cc$season[order(cc$season)])
  cc$season = factor(cc$season, levels = seasons)
  cty_order = cc %>% group_by(country) %>% summarise(t = sum(dll), .groups = "drop") %>% arrange(t) %>% pull(country)
  cc$country = factor(cc$country, levels = cty_order)
  lim = max(abs(cc$dll))
  pb = ggplot(cc, aes(season, country, fill = dll)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%+.1f", dll)), size = 2.6, colour = "grey15") +
    scale_fill_gradient2(low = blue, mid = "grey95", high = orange, midpoint = 0, limits = c(-lim, lim),
                         name = "log-likelihood,\nR0 minus S0") +
    labs(title = sprintf("(b) Wave by wave: does the SEASON effect prefer R0 or S0?  total %+.1f nats, %d of %d waves favour R0",
                         sum(cc$dll), sum(cc$dll > 0), nrow(cc)),
         subtitle = paste("The same fit with the season effect moved from S0 to R0 (country effect on S0 in both). Orange: the wave is",
                          "better explained by a faster rise WITHOUT a bigger pool; blue: by a bigger pool as well. Countries",
                          "ordered by their total.", sep = "\n"),
         x = NULL, y = NULL) +
    theme_minimal(10) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
                              plot.title = element_text(face = "bold"))

  # (c) the most discriminating wave, both fits over the data
  k = which.max(abs(cc$dll)); ic_name = as.character(cc$country[k]); s_name = as.character(cc$season[k])
  curves = do.call(rbind, lapply(c("S0/S0", "R0/S0"), function(m){
    fit = fits[[m]]; d = fit$d
    i = which(d$countries[d$cs_country + 1L] == ic_name & d$seasons[d$cs_season + 1L] == s_name)
    ic = d$cs_country[i] + 1L; N = d$N[[ic]]; per = d$rate_per / N
    mu = jm_fitted_cpp(fit$theta, d)$mu[[i]]
    data.frame(model = m, week = seq_len(nrow(mu)),
               value = rowSums(sweep(mu, 2, per, "*") * matrix(N / sum(N), nrow(mu), 3, byrow = TRUE)))
  }))
  d = fits[["S0/S0"]]$d
  i = which(d$countries[d$cs_country + 1L] == ic_name & d$seasons[d$cs_season + 1L] == s_name)
  ic = d$cs_country[i] + 1L; N = d$N[[ic]]; y = d$y[[i]]
  obs = data.frame(week = seq_len(nrow(y)),
                   value = rowSums(sweep(y, 2, d$rate_per / N, "*") * matrix(N / sum(N), nrow(y), 3, byrow = TRUE)))
  pc = ggplot() +
    geom_point(data = obs, aes(week, value), colour = "grey25", size = 1.4, na.rm = TRUE) +
    geom_line(data = curves, aes(week, value, colour = model), linewidth = 1) +
    scale_colour_manual(values = c("S0/S0" = blue, "R0/S0" = orange),
                        labels = c("S0/S0" = "season effect on S0 (working model)", "R0/S0" = "season effect on R0"), name = NULL) +
    labs(title = sprintf("(c) The wave that discriminates most: %s %s (%+.1f nats for R0)", ic_name, s_name, cc$dll[k]),
         subtitle = paste("Points are the observed ILI+ per 100 000, lines the two fitted means. The two mechanisms reach the peak",
                          "with a different rise-to-decay shape: an S0 effect that makes the rise faster also deepens the pool,",
                          "so the decay differs. This is the shape the likelihood is scoring across all 85 waves.", sep = "\n"),
         x = "week of the season (from 1 August)", y = "ILI+ per 100 000") +
    theme_minimal(10) + theme(legend.position = "top", plot.title = element_text(face = "bold"))

  fig = (pa | pc) / pb + plot_layout(heights = c(1, 1.25)) +
    plot_annotation(
      title = "Where does the variation live: in the susceptible pool, or in transmissibility?",
      subtitle = "An effect on S0 changes how fast a wave rises AND how many are left to infect; an effect on R0 changes the speed only. The data can tell these apart.",
      theme = theme(plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 10, colour = "grey30")))
  if (!is.null(out)) ggsave(out, fig, width = 15, height = 12.5, dpi = 115)
  fig
}

if (sys.nframe() == 0){
  setwd(here::here())
  suppressMessages(source("code/01_main_supporting/setup.R"))
  source("code/01_main_supporting/stitch_iliplus.R"); source("code/01_main_supporting/sir_core.R")
  source("code/06_comp_model/contact_matrix.R"); source("code/06_comp_model/comp_model_settings.R")
  source("code/06_comp_model/comp_model_core.R"); source("code/06_comp_model/comp_model_data.R")
  source("code/07_joint_model/joint_model.R"); jm_load_cpp()
  cmp = readRDS("output/joint_model/compare/comparison.rds")
  fits = list(`S0/S0` = readRDS("output/joint_model/compare/fit_S0_S0.rds"),
              `R0/S0` = readRDS("output/joint_model/compare/fit_R0_S0.rds"))
  plot_model_comparison(cmp, fits, out = "output/joint_model/17_model_comparison.png")
  cat("wrote output/joint_model/17_model_comparison.png\n")
}
