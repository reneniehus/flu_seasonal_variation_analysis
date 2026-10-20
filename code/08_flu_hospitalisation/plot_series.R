# Every influenza-hospitalisation series over time, one row per series, grouped by country.
# Units differ between series, so each row is scaled to its own maximum over the window (a sequential
# blue ramp: light = low, dark = that series' peak). A blank stretch means no data. Rows of the same
# country stacked together show whether national and ERVISS series agree on timing.
#
# Usage: Rscript code/08_flu_hospitalisation/plot_series.R <data dir>
#   reads  <data dir>/erviss_flu_hospital_weekly.csv and national_flu_hospital_series.csv
#   writes <data dir>/figures/series_by_country.png
suppressMessages({library(ggplot2); library(dplyr); library(ISOweek)})

dir    <- commandArgs(TRUE)[1]
window <- as.Date(c("2014-08-01", "2026-08-01"))   # seasons 2014/15 to 2025/26
ink    <- c(primary = "#0b0b0b", secondary = "#52514e", grid = "#e1e0d9", surface = "#fcfcfb")
names_iso <- c(AT = "Austria", DE = "Germany", DK = "Denmark", ES = "Spain", FR = "France", PL = "Poland", SE = "Sweden")
icu_like  <- c("Finland", "France", "Italy", "Netherlands", "Portugal", "Sweden")   # ERVISS admissions = ICU there

# ERVISS: admissions, ICU and SARI+ (all ages); where admissions are the ICU series, keep the ICU row only
e <- read.csv(file.path(dir, "erviss_flu_hospital_weekly.csv")) %>%
  filter(age == "total", subtype == "total", country != "EU/EEA",
         indicator %in% c("hosp_admissions", "icu_admissions", "sari_plus"),
         !(indicator == "hosp_admissions" & country %in% icu_like)) %>%
  mutate(series = case_when(indicator == "hosp_admissions" ~ "admissions, ERVISS [A]",
                            indicator == "icu_admissions"  ~ "ICU admissions, ERVISS [F]",
                            country == "Germany"           ~ "SARI coded as influenza, ERVISS [C]",
                            TRUE                           ~ "SARI+, ERVISS [B]"),
         start = ISOweek2date(paste0(yearweek, "-1")), end = start + 7) %>%
  select(country, series, start, end, value)

# national series: weekly, monthly (SE) and quarter-monthly (PL, days 1-7, 8-15, 16-22, 23-end)
n <- read.csv(file.path(dir, "national_flu_hospital_series.csv"))
first_of_next_month <- function(d) as.Date(format(as.Date(format(d, "%Y-%m-01")) + 32, "%Y-%m-01"))
n$start <- as.Date(NA); n$end <- as.Date(NA)
wk <- n$period_type == "week"
n$start[wk] <- ISOweek2date(paste0(n$period[wk], "-1")); n$end[wk] <- n$start[wk] + 7
mo <- n$period_type == "month"
n$start[mo] <- as.Date(paste0(n$period[mo], "-01")); n$end[mo] <- first_of_next_month(n$start[mo])
qm <- n$period_type == "quarter-month"
q  <- substr(n$period[qm], 8, 8); ym <- substr(n$period[qm], 1, 7)
n$start[qm] <- as.Date(paste0(ym, "-", c(A = "01", B = "08", C = "16", D = "23")[q]))
n$end[qm]   <- as.Date(ifelse(q == "D", as.character(first_of_next_month(n$start[qm])),
                              paste0(ym, "-", c(A = "08", B = "16", C = "23")[q])))
n <- n %>%
  mutate(country = names_iso[country],
         series = recode(indicator,
                         lab_confirmed_admissions = "lab-confirmed admissions, SSI [A]",
                         SARI_influenza_rate = "SARI+, SiVIRA [B]",
                         SARI_influenza_coded_rate = "SARI coded as influenza, RKI [C]",
                         SARI_influenza_coded_admissions = "SARI coded as influenza, SARI-Dashboard [C]",
                         ED_to_hospital_influenza = "after emergency visit, OSCOUR [D]",
                         `inpatient_episodes_ICD_J09-J11` = "J09-J11 inpatients (monthly), Socialstyrelsen [E]",
                         hospital_referrals_influenza_or_suspected = "referrals, suspected influenza, PZH [G]")) %>%
  select(country, series, start, end, value)

d <- bind_rows(e, n) %>%
  filter(start >= window[1], end <= window[2], !is.na(value)) %>%
  group_by(country, series) %>%
  filter(n() >= 20, max(value) > 0) %>%
  mutate(scaled = value / max(value)) %>%
  ungroup() %>%
  mutate(class = substr(series, nchar(series) - 1, nchar(series) - 1)) %>%
  arrange(country, class, series)
d$series <- factor(d$series, levels = rev(unique(d$series)))
n_rows   <- n_distinct(paste(d$country, d$series))

ramp <- c("#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b")   # sequential blue, steps 100-700
g <- ggplot(d, aes(x = start + (end - start) / 2, y = series)) +
  geom_tile(aes(width = as.numeric(end - start), fill = scaled), height = 0.78) +
  facet_grid(country ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(colours = ramp, limits = c(0, 1), name = "Share of the series' maximum",
                       guide = guide_colourbar(barwidth = 12, barheight = 0.5, title.position = "top")) +
  scale_x_date(limits = window, expand = c(0, 0),
               breaks = as.Date(paste0(2015:2026, "-02-01")),
               labels = paste0(14:25, "/", 15:26),
               minor_breaks = as.Date(paste0(2014:2026, "-08-01"))) +
  labs(x = "Season", y = NULL,
       title = "Influenza hospitalisation indicators over time, one row per series",
       subtitle = paste0("Each row is scaled to its own maximum over 2014/15 to 2025/26; blank means no data.\n",
                         "The letter gives the indicator class, as in the coverage map. ",
                         "Sweden's monthly series reaches back to 1998 (cut here)."),
       caption = "Sources: ERVISS (ECDC), SSI, SARI-Dashboard Austria, RKI, ISCIII SiVIRA, Sante publique France, Socialstyrelsen, NIZP-PZH.") +
  theme_minimal(base_size = 9) +
  theme(plot.background = element_rect(fill = ink[["surface"]], colour = NA),
        text = element_text(colour = ink[["primary"]]),
        axis.text = element_text(colour = ink[["secondary"]]),
        panel.grid.major = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_line(colour = ink[["grid"]], linewidth = 0.3),
        panel.spacing.y = unit(3, "pt"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", colour = ink[["primary"]]),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = ink[["secondary"]]),
        plot.caption = element_text(colour = ink[["secondary"]], hjust = 0),
        plot.caption.position = "plot", plot.title.position = "plot",
        legend.position = "top", legend.justification = "left")

dir.create(file.path(dir, "figures"), showWarnings = FALSE)
ggsave(file.path(dir, "figures", "series_by_country.png"), g, width = 11, height = 1.6 + 0.19 * n_rows, dpi = 150)
cat("series figure written:", n_rows, "series\n")
