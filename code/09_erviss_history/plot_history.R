# The merged ERVISS primary-care series over time, one row per country, 2014/15 to 2025/26.
#   ILI rates: each row scaled to its own 99th-percentile week (units differ between countries), so a
#   single erroneous week (Slovakia 2015-W51: 155 763 per 100 000) cannot wash out the rest of its row.
#   Influenza positivity: the share of tests positive, on one scale for all countries, from the
#   stream the ILI+ panel uses (non-sentinel for HR, FI, IS, LV, MT, RO; sentinel elsewhere), and
#   only in weeks with at least 5 tests, so single-test weeks do not show as 0% or 100%.
# A blank means no data. Sequential blue: light = low, dark = high.
#
# Usage: Rscript code/09_erviss_history/plot_history.R <repository data dir>
#   writes <data>/primary_care_ili_rates/figures/ili_rates_by_country.png
#          <data>/primary_care_flu_positivity/figures/positivity_by_country.png
suppressMessages({library(ggplot2); library(dplyr); library(ISOweek)})

dir    <- commandArgs(TRUE)[1]
window <- as.Date(c("2014-08-01", "2026-08-01"))
ink    <- c(primary = "#0b0b0b", secondary = "#52514e", grid = "#e1e0d9", surface = "#fcfcfb")
ramp   <- c("#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b")   # sequential blue, steps 100-700
nonsentinel <- c("Croatia", "Finland", "Iceland", "Latvia", "Malta", "Romania")

strips <- function(d, fill_name, limits, title, subtitle, file) {
  d <- d %>% filter(start >= window[1], start + 7 <= window[2]) %>% mutate(country = factor(country, levels = rev(sort(unique(country)))))
  g <- ggplot(d, aes(x = start + 3.5, y = country, fill = fill)) +
    geom_tile(width = 7, height = 0.78) +
    scale_fill_gradientn(colours = ramp, limits = limits, oob = scales::squish, name = fill_name,
                         guide = guide_colourbar(barwidth = 12, barheight = 0.5, title.position = "top")) +
    scale_x_date(limits = window, expand = c(0, 0), breaks = as.Date(paste0(2015:2026, "-02-01")),
                 labels = paste0(14:25, "/", 15:26), minor_breaks = as.Date(paste0(2014:2026, "-08-01"))) +
    labs(x = "Season", y = NULL, title = title, subtitle = subtitle,
         caption = "Source: ERVISS (ECDC), merged across its releases of 2024-10-18, 2025-06-20, 2026-05-22 and 2026-10-09.") +
    theme_minimal(base_size = 9) +
    theme(plot.background = element_rect(fill = ink[["surface"]], colour = NA),
          text = element_text(colour = ink[["primary"]]), axis.text = element_text(colour = ink[["secondary"]]),
          panel.grid.major = element_blank(), panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_line(colour = ink[["grid"]], linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(colour = ink[["secondary"]]),
          plot.caption = element_text(colour = ink[["secondary"]], hjust = 0),
          plot.caption.position = "plot", plot.title.position = "plot",
          legend.position = "top", legend.justification = "left")
  dir.create(dirname(file), showWarnings = FALSE)
  ggsave(file, g, width = 11, height = 1.8 + 0.2 * n_distinct(d$country), dpi = 150)
  cat("written:", file, "\n")
}

ili <- read.csv(file.path(dir, "primary_care_ili_rates", "erviss_iliari_merged.csv")) %>%
  filter(indicator == "ILIconsultationrate", age == "total") %>%
  transmute(country = countryname, start = ISOweek2date(paste0(yearweek, "-1")), value) %>%
  group_by(country) %>% mutate(fill = value / quantile(value, 0.99)) %>% ungroup()
strips(ili, "Share of the country's 99th-percentile week (higher values darkest)", c(0, 1),
       "Primary care ILI rates, all ages, by country",
       "Each row is scaled to its own 99th-percentile week (units differ: see the README). Blank means no data.",
       file.path(dir, "primary_care_ili_rates", "figures", "ili_rates_by_country.png"))

read_tests <- function(f) read.csv(file.path(dir, "primary_care_flu_positivity", f)) %>%
  filter(pathogentype == "Influenza", age == "total", indicator %in% c("tests", "detections")) %>%
  select(countryname, yearweek, indicator, value) %>%
  tidyr::pivot_wider(names_from = indicator, values_from = value)
pos <- bind_rows(read_tests("erviss_detections_sentinel_merged.csv") %>% filter(!countryname %in% nonsentinel),
                 read_tests("erviss_detections_nonsentinel_merged.csv") %>% filter(countryname %in% nonsentinel)) %>%
  filter(countryname != "EU/EEA", tests >= 5) %>%
  transmute(country = ifelse(countryname %in% nonsentinel, paste0(countryname, " (non-sentinel)"), countryname),
            start = ISOweek2date(paste0(yearweek, "-1")), fill = 100 * coalesce(detections, 0) / tests)
strips(pos, "Influenza positivity (%), 60% and above darkest", c(0, 60),
       "Influenza test positivity in primary care, by country",
       paste0("Share of tests positive, in weeks with at least 5 tests; a week without a detections row counts as 0.\n",
              "Non-sentinel tests where the ILI+ panel uses them. Blank means no data or fewer than 5 tests."),
       file.path(dir, "primary_care_flu_positivity", "figures", "positivity_by_country.png"))
