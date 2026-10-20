# Coverage map: the best influenza-hospitalisation indicator for each country and season.
# Colour carries the kind of evidence (three families); the letter in each tile names the exact class,
# so identity never rests on colour alone. Pale tiles hold fewer periods than a full season needs.
#
# Usage: Rscript code/08_flu_hospitalisation/plot_coverage.R <data dir>
#   reads  <data dir>/coverage_by_country_season.csv
#   writes <data dir>/figures/coverage_by_country_season.png
suppressMessages({library(ggplot2); library(dplyr)})

dir  <- commandArgs(TRUE)[1]
d    <- read.csv(file.path(dir, "coverage_by_country_season.csv"))
d    <- d %>% filter(season >= "2010/2011", season <= "2025/2026")
ink  <- c(primary = "#0b0b0b", secondary = "#52514e", surface = "#fcfcfb")
fam  <- c("Laboratory-based: admissions or SARI+ (A, B)"      = "#2a78d6",
          "Laboratory-confirmed: ICU or severe cases only (F)" = "#eb6834",
          "Diagnosis code or clinical diagnosis (C, D, E, G)"  = "#1baf7a")

d <- d %>%
  mutate(letter = substr(best_class, 1, 1),
         family = factor(case_when(letter %in% c("A", "B") ~ names(fam)[1],
                                   letter == "F"           ~ names(fam)[2],
                                   TRUE                    ~ names(fam)[3]), levels = names(fam)),
         full   = periods >= ifelse(letter == "E", 5, 20),
         season = sub("/20", "/", season))
ord <- d %>% group_by(country) %>% summarise(n = sum(periods)) %>% arrange(n) %>% pull(country)
d$country <- factor(d$country, levels = ord)

g <- ggplot(d, aes(season, country)) +
  geom_tile(aes(fill = family, alpha = full), colour = ink[["surface"]], linewidth = 0.7) +
  geom_text(aes(label = letter), size = 2.7, colour = ink[["primary"]]) +
  scale_fill_manual(values = fam, name = NULL, drop = FALSE) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  labs(x = "Season", y = NULL,
       title = "Influenza hospitalisation data by country and season: ERVISS and national sources",
       subtitle = paste0("Each tile shows the best indicator available. Pale tiles: fewer than 20 weeks of it ",
                         "(5 months for E).\nSeasons 2019/20 to 2022/23 are affected by COVID-19. ",
                         "Sweden's monthly series starts in 1997/98 (not shown)."),
       caption = paste0("A laboratory-confirmed admissions  |  B SARI rate x laboratory influenza positivity (SARI+)  |  ",
                        "C SARI coded as influenza\nD hospitalisations after an emergency visit  |  ",
                        "E inpatient diagnoses J09-J11 (monthly)  |  F ICU or severe cases only  |  ",
                        "G hospital referrals, suspected influenza")) +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = ink[["surface"]], colour = NA),
        panel.grid = element_blank(),
        text = element_text(colour = ink[["primary"]]),
        axis.text = element_text(colour = ink[["secondary"]]),
        axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = ink[["secondary"]]),
        plot.caption = element_text(colour = ink[["secondary"]], hjust = 0),
        plot.caption.position = "plot", plot.title.position = "plot",
        legend.position = "top", legend.justification = "left", legend.direction = "vertical")

dir.create(file.path(dir, "figures"), showWarnings = FALSE)
ggsave(file.path(dir, "figures", "coverage_by_country_season.png"), g, width = 10, height = 9, dpi = 150)
cat("coverage map written\n")
