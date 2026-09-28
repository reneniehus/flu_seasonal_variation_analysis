# report_figures.R -- the figures of the printable report (documentation/joint_model_report.pdf),
# drawn from the working fit for A4. The default figure set (joint_report.R) explains itself inside
# each figure; here the report's captions carry that, so titles and explanations are dropped, the type
# is Source Sans 3 at print size (falling back to the system sans), and every figure is drawn on a
# 9-inch canvas that prints at the 170 mm text width. One figure is new: what informs the children's
# modifier, the surveillance data or the cohort evidence.
#   Rscript code/07_joint_model/report/report_figures.R      (then build_report.sh)
setwd(here::here())
suppressMessages(source("code/01_main_supporting/setup.R"))
source("code/07_joint_model/joint_model.R"); source("code/07_joint_model/joint_recovery.R")
source("code/07_joint_model/joint_report.R"); jm_load_cpp()
suppressMessages(library(patchwork))
o = readRDS("output/joint_model/joint_fit.rds"); fit = o$fit; d = fit$d; iv = o$intervals
out = "output/joint_model/report"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
FONT = if ("Source Sans 3" %in% systemfonts::system_fonts()$family) "Source Sans 3" else "sans"

# every plot function draws with .jm_theme(): redefined here for print, one family, type 15% larger
# on the 9-inch canvas so that it prints at about 8.5 pt
.jm_theme = function(base = 10){
  theme_minimal(base_size = 1.15 * base, base_family = FONT) +
    theme(plot.title = element_text(face = "bold", size = rel(1.0)),
          strip.text = element_text(face = "bold", size = rel(0.85)),
          panel.grid.minor = element_blank(), legend.position = "top")
}
bare = function(p) p + labs(title = NULL, subtitle = NULL, caption = NULL)
unlabelled = function(p){ p$layers = Filter(function(l) !inherits(l$geom, "GeomText"), p$layers); p }  # the tables carry the values
ages = c(young = "children 0-14", medium = "adults 15-64", elderly = "65+")
save = function(p, name, h, w = 9)
  ggsave(file.path(out, name), p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")

# ---- the data ----
save(bare(plot_jm_data_panel(fit)), "fig01_data_panel.png", 4.3)
feat = (plot_jm_data_features(fit) & labs(subtitle = NULL)) + plot_annotation(title = NULL, subtitle = NULL)
feat[[2]] = feat[[2]] + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save(feat, "fig02_data_features.png", 6.6)

# ---- the model ----
mech = plot_jm_mechanism(fit)
short = c("x_s  season effect on susceptibility" = "season effect on S0",
          "S0_c  the country's susceptibility level" = "country level of S0",
          "I0  seed, i.e. when the wave arrives" = "seed (arrival)",
          "c  reporting level of the country" = "country reporting level",
          "delta  visibility of the season" = "season visibility",
          "off_eld  age reporting, 65+" = "65+ reporting (65+ shown)",
          "kappa  children's extra susceptibility" = "children's S0 modifier (0-14 shown)",
          "tau  spread of the wave across the country" = "spatial spread (fixed at 0)")
pm = mech[[1]]
pm$layers = Filter(function(l) !inherits(l$geom, "GeomText"), pm$layers)   # the caption explains each panel
pm = pm + facet_wrap(~ key, scales = "free_y", ncol = 4, labeller = as_labeller(short)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08)))
save(pm, "fig03_mechanism.png", 4.2)
save(bare(mech[[2]]), "fig04_ties.png", 3.0)
save(bare(plot_jm_design(fit)), "fig05_design.png", 4.4)

# ---- what informs the children's modifier: the surveillance data or the cohort ----
prof = readRDS("output/joint_model/kappa/profile_kappa.rds")$table
kf = iv[iv$parameter == "log_kappa_young", ]
f0 = readRDS("output/joint_model/kappa/fit_kappa_0.rds")
att = jm_fitted_cpp(fit$theta, d)$attack
fitted = data.frame(k = kf$estimate, lo = kf$lower, hi = kf$upper,
                    dll = fit$loglik - f0$loglik, ratio = median(att[, 1] / att[, 2]))
pr = exp(d$pr_kappa_mean + c(-1.96, 1.96) * d$pr_kappa_sd)     # the prior's 95% range
pa = ggplot(prof, aes(kappa, d_loglik)) +
  annotate("rect", xmin = pr[1], xmax = pr[2], ymin = -Inf, ymax = Inf, fill = "grey88") +
  annotate("text", x = mean(pr), y = Inf, label = "cohort prior, 95%", vjust = 1.6, size = 3.1, colour = "grey35", family = FONT) +
  geom_line(colour = .jm_blue, linewidth = 0.7) + geom_point(colour = .jm_blue, size = 2.2) +
  geom_segment(data = fitted, aes(x = lo, xend = hi, y = dll, yend = dll), colour = .jm_orange, linewidth = 0.9) +
  geom_point(data = fitted, aes(k, dll), colour = .jm_orange, size = 3) +
  annotate("text", x = fitted$hi + 0.08, y = fitted$dll, label = "fitted", hjust = 0, size = 3.1, colour = .jm_orange, family = FONT) +
  labs(title = "(a) The surveillance data: fit against no modifier",
       x = "children's S0 modifier, kappa", y = "log-likelihood difference (nats)") + .jm_theme(10)
pb = ggplot(prof, aes(kappa, young_adult)) +
  annotate("rect", xmin = pr[1], xmax = pr[2], ymin = -Inf, ymax = Inf, fill = "grey88") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.65, ymax = 1.80, fill = .jm_green, alpha = 0.22) +
  annotate("text", x = 0.02, y = 1.84, label = "PHIRST cohort, 1.65-1.80", hjust = 0, vjust = 0, size = 3.1, colour = "grey25", family = FONT) +
  geom_line(colour = .jm_blue, linewidth = 0.7) + geom_point(colour = .jm_blue, size = 2.2) +
  geom_segment(data = fitted, aes(x = lo, xend = hi, y = ratio, yend = ratio), colour = .jm_orange, linewidth = 0.9) +
  geom_point(data = fitted, aes(k, ratio), colour = .jm_orange, size = 3) +
  annotate("text", x = fitted$hi + 0.08, y = fitted$ratio, label = "fitted", hjust = 0, size = 3.1, colour = .jm_orange, family = FONT) +
  labs(title = "(b) The cohort: children's infection relative to adults'",
       x = "children's S0 modifier, kappa", y = "children / adults, attack rate") + .jm_theme(10)
save(pa + pb, "fig06_kappa.png", 3.3)
save(bare(plot_jm_identifiability(o$id)), "fig07_data_or_prior.png", 3.2)

# ---- what it learned ----
tagged = function(p) (p + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")) &
  theme(plot.tag = element_text(family = FONT, face = "bold", size = 11))
save(tagged((unlabelled(bare(plot_jm_season_S0(fit, iv))) + labs(y = "adult S0, typical country")) /
            (unlabelled(bare(plot_jm_season_visibility(fit, iv))) + labs(y = "visibility (log scale)"))),
     "fig08_seasons.png", 5.6)
save(tagged((unlabelled(bare(plot_jm_country_S0(fit, iv))) + labs(x = "adult S0")) |
            unlabelled(bare(plot_jm_country_reporting(fit, iv)))),
     "fig09_countries.png", 3.9)
pa = bare(plot_jm_attack(fit)) + scale_colour_manual(values = .jm_gcol, labels = ages, name = NULL)
pb = bare(plot_jm_age_offsets(fit, iv)) + scale_colour_manual(values = .jm_gcol[c("young", "elderly")], labels = ages, name = NULL)
save(tagged(pa | pb), "fig10_age.png", 4.0)

# ---- how well it fits ----
save(bare(plot_jm_fit_overview(fit)) + theme(axis.text = element_text(size = 6.2), strip.text = element_text(size = 7.5)),
     "fig11_fit_overview.png", 10.0)
save(bare(plot_jm_fit_country(fit, "DK")) + facet_grid(group ~ season, scales = "free_y", labeller = labeller(group = ages)) +
       theme(axis.text = element_text(size = 6.8), strip.text = element_text(size = 7.8)),
     "fig12_fit_DK.png", 4.3)
save(bare(plot_jm_adequacy(fit)), "fig13_noise_budget.png", 3.4)
save(bare(plot_jm_arrival(fit)), "fig14_seeds.png", 3.8)

# ---- whether to trust it ----
rec = readRDS("output/joint_model/joint_recovery.rds")
save(bare(plot_jm_recovery(rec$local$rec, d, rec$local$summary)), "fig15_recovery.png", 5.2)
cat("report figures written to", out, "\n")
