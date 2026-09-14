# stitch_iliplus.R
#
# THE single implementation of the combined ILI+ panel-assembly rules (see
# documentation/decisions.md). Before this file existed the same stitch was copy-pasted in
# build_slim_panel.R, precovid_predict_postcovid.R and (approximately) data_availability.R,
# and the copies had started to diverge -- every consumer now calls stitch_iliplus_panel().
#
# The rules it encodes:
#   - ERVISS ILI+ = ILI rate x positivity, SENTINEL except NON-sentinel for MT/IS/HR/RO/LV/FI;
#     per-100-consultations countries CY/LU/MT scaled x1000 onto the per-100 000 basis;
#   - stitch PER WEEK: RespiCompass where present, else ERVISS aligned to the RespiCompass scale
#     by a per-country factor (median RespiCompass/ERVISS over the 2023/24 overlap weeks). In
#     practice RespiCompass covers <=2023/24 and ERVISS 2024/25+, but within the 2023/24 overlap
#     season aligned-ERVISS weeks DO fill RespiCompass gaps (14 country-seasons mix sources);
#     the per-season `source` label is the FIRST contributing source = the dominant early-season
#     one (RespiCompass for all current mixed seasons);
#   - single-source overrides: NO,ES = RespiCompass only; SK,LV = ERVISS only (native scale);
#   - optionally exclude the four COVID seasons; require >= min_wk strictly-POSITIVE observed
#     weeks (a season of zeros carries no wave); contiguous weekly grid from season week 1.
#
# Requires the tidyverse (source setup.R first). models_in comes from output/models_in.rds.

stitch_covid_seasons <- c("2019/2020", "2020/2021", "2021/2022", "2022/2023")

# ---- |-the stitch rules as data (one definition, shared by every stitch function below) ----
.stitch_rules <- list(
  nonsentinel = c("MT", "IS", "HR", "RO", "LV", "FI"),   # ERVISS ILI+ uses NON-sentinel positivity
  per_1000    = c("CY", "LU", "MT"),                     # per-100-consultations -> x1000 to per-100 000
  resp_only   = c("NO", "ES"),                           # single source = RespiCompass
  erviss_only = c("SK", "LV")                            # single source = ERVISS
)

# ---- |-week-level stitched values for one or more age groups (the shared core) ----
# Applies the source rules per week and age group; the per-country alignment factor is ALWAYS
# estimated on the age TOTAL (the panel's definition) and applied to every band, so age-specific
# series stay on the panel's scale. Returns country_short, season, date, season_week, agegroup,
# value, source (all requested agegroups; no season filtering, no week grid).
.stitch_iliplus_weeks <- function(models_in, agegroups = "age_total", overlap_season = "2023/2024"){
  r <- .stitch_rules
  ili_plus <- models_in$data_timeseries_long %>%
    filter(indicator=="ili_plus", pathogen=="Influenza", agegroup %in% union("age_total", agegroups)) %>%
    select(stream, country_short, season, date, season_week, agegroup, value)

  erviss <- ili_plus %>% filter(stream %in% c("ili_plus_sentinel","ili_plus_nonsentinel")) %>%
    mutate(chosen_stream = ifelse(country_short %in% r$nonsentinel, "ili_plus_nonsentinel", "ili_plus_sentinel")) %>%
    filter(stream==chosen_stream) %>%
    mutate(value = value * ifelse(country_short %in% r$per_1000, 1000, 1)) %>%
    transmute(country_short, season, date, season_week, agegroup, erviss = value)
  respicompass <- ili_plus %>% filter(stream=="ili_plus_respicompass") %>%
    transmute(country_short, season, date, season_week, agegroup, respicompass = value)

  combined <- full_join(erviss, respicompass, by=c("country_short","season","date","season_week","agegroup"))

  # per-country alignment factor from the overlap season, on the TOTAL (median RespiCompass / ERVISS over weeks where BOTH streams are finite and > 0)
  align_factors <- combined %>% filter(agegroup=="age_total", season==overlap_season, is.finite(erviss), erviss>0, is.finite(respicompass), respicompass>0) %>%
    group_by(country_short) %>% summarise(align_factor = median(respicompass/erviss), .groups="drop")

  combined %>% left_join(align_factors, by="country_short") %>%
    mutate(align_factor = ifelse(is.na(align_factor), 1, align_factor),
           value  = case_when(country_short %in% r$resp_only   ~ respicompass,
                              country_short %in% r$erviss_only ~ erviss,             # native ERVISS scale
                              !is.na(respicompass) ~ respicompass,                   # default: RespiCompass where present
                              TRUE        ~ erviss * align_factor),                  # default ERVISS era, aligned to RespiCompass
           source = case_when(country_short %in% r$resp_only   ~ "RespiCompass",
                              country_short %in% r$erviss_only ~ "ERVISS",
                              !is.na(respicompass) ~ "RespiCompass", TRUE ~ "ERVISS")) %>%
    filter(agegroup %in% agegroups) %>%
    select(country_short, season, date, season_week, agegroup, value, source, align_factor)
}

stitch_iliplus_panel <- function(models_in, exclude_covid = TRUE, min_wk = 15,
                                 overlap_season = "2023/2024"){
  combined <- .stitch_iliplus_weeks(models_in, "age_total", overlap_season) %>%
    select(country_short, season, date, season_week, value, source)
  if (exclude_covid) combined <- combined %>% filter(!season %in% stitch_covid_seasons)
  combined <- combined %>% filter(is.finite(value))

  # keep country-seasons with >= min_wk positive observed weeks; contiguous weekly grid from week 1
  combined %>%
    group_by(country_short, season) %>%
    filter(sum(is.finite(value) & value>0) >= min_wk) %>%
    group_modify(function(df, key){
      season_source <- df$source[which(!is.na(df$source))][1]   # first contributing source = the dominant early-season one (overlap seasons can mix sources; see header)
      last_week     <- max(df$season_week[is.finite(df$value)])
      tibble(season_week = 1:last_week) %>%
        left_join(df %>% select(season_week, date, value), by="season_week") %>%
        arrange(season_week) %>%
        transmute(week = season_week, season_week, date, value, source = season_source)
    }) %>% ungroup() %>%
    arrange(country_short, season, week)
}

# ---- |-age-specific stitched series for the committed panel's country-seasons ----
# For the compartmental model (code/06_comp_model): the same per-week source rules and the same
# total-based alignment factor as the panel, evaluated for the four ILI age bands, restricted to
# the country-seasons that are IN the committed panel (inclusion decided on the total, as the panel
# does), and laid on the panel's contiguous weekly grid (weeks absent from a band are NA).
# Returns country_short, season, week, agegroup, value (rate per 100 000 of the band), source.
stitch_iliplus_by_age <- function(models_in, panel = read.csv("data/slim_flu_iliplus.csv", stringsAsFactors=FALSE),
                                  agegroups = c("age_00_04","age_05_14","age_15_64","age_65_99")){
  weeks <- .stitch_iliplus_weeks(models_in, agegroups) %>% filter(is.finite(value)) %>%
    select(country_short, season, week = season_week, agegroup, value)
  grid <- panel %>% distinct(country_short, season, week, source)          # the panel's grid + season source label
  tidyr::crossing(grid, agegroup = agegroups) %>%
    left_join(weeks, by = c("country_short","season","week","agegroup")) %>%
    arrange(country_short, season, agegroup, week)
}
