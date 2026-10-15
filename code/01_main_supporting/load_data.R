# 1: small helpers shared by the data-loading functions (incl. the cache helper)
# 2: defining functions that load each data stream
# 3: define a mother function that calls each data stream function
#
# Conventions used by every loader below:
#   data             - the growing list of data streams; each loader adds one named slot and returns it
#   params           - run settings (see settings_version0.R); threaded through so any loader MIGHT use it
#   regenerate       - F: reuse the cached output/<name>.Rdata if present; T: rebuild it from scratch
#   new_from_online  - T: (re)fetch from the internet and refresh the local snapshot;
#                      F: read the local snapshot, self-bootstrapping (fetch on the fly) if it is missing

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Helpers shared across the data-loading functions ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-cache helper: load a stream from cache, or (re)build and cache it ----
# Every data stream shares the same caching dance, so it lives here once instead of being
# copy-pasted into each loader. Pass a build_fn() that constructs the object from scratch:
#   - if output/<name>.Rdata exists and regenerate==F  -> return the cached object
#   - otherwise                                         -> run build_fn(), cache it, return it
# The cache is read into a throw-away environment so we recover the object whatever variable
# name it was saved under (keeps us compatible with older caches).
load_or_build = function(name, build_fn, regenerate=FALSE){
  cache_path = here("output", paste0(name, ".Rdata"))
  if (file.exists(cache_path) & !regenerate) {
    pr=paste0("Loading '",name,"' from cache (output/",name,".Rdata) ... \n"); cat(green(pr))
    cache_env = new.env()
    load(cache_path, envir=cache_env)
    return(get(ls(cache_env)[1], envir=cache_env))
  }
  obj = build_fn()                  # stream-specific construction
  save(obj, file=cache_path)        # refresh the on-disk cache
  return(obj)
}

# ---- |-build a RespiCompass raw-file URL for the active hub round ----
# the round folder changes each season, so keep it in one place (params$respicompass_round)
# rather than repeating the full path at every read_csv() call.
respicompass_file = function(relpath, params=NULL){
  round = if (!is.null(params$respicompass_round)) params$respicompass_round else "2024-2025_round_1"
  paste0("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/",
         "refs/heads/main/Previous_Rounds/", round, "/", relpath)
}

# ---- |-recode the ERVISS age labels to our age-group labels ----
recode_age = function(age){
  case_when(
    age == "0-4"   ~ "age_00_04",
    age == "5-14"  ~ "age_05_14",
    age == "15-64" ~ "age_15_64",
    age == "65+"   ~ "age_65_99",
    age == "total" ~ "age_total",
    age == "unk"   ~ "age_unk",
    .default = age          # carry through any age label we have not mapped yet
  )
}

# ---- |-standardise one raw ERVISS table ----
# every ERVISS file shares the core columns: survtype, countryname, yearweek, indicator, age, value
# here we (i) add an ISO-week (Wednesday) date, (ii) recode age, (iii) add the short country code,
# then either slim down to the rate columns or keep every column.
standardise_erviss = function(df, schema){
  df = df %>%
    mutate(date          = ISOweek2date(paste0(yearweek,"-3"))) %>%   # mid-week (Wednesday) date
    mutate(age           = recode_age(age)) %>%                       # standard age-group labels (in place)
    mutate(country_short = EU_short(countryname))                     # ISO2-style short country code
  if (schema == "rates") {
    # slim consultation/hospitalisation-rate schema (no pathogen breakdown) used by the SIR-type models
    df = df %>% select(country_short, date, target=indicator, agegroup=age, value)
  }
  # schema == "detailed": keep every original column (pathogen / pathogentype / pathogensubtype /
  # datasource / variant / detectableprevalence are all carried through unchanged)
  return(df)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Defining data-loading functions for each data stream ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load_data_epi = function(data=list(), params=NULL, regenerate=FALSE, new_from_online=TRUE){

  build_epi = function(){

    # ---- |-ERVISS file registry ----
    # one row per dataset published in the ECDC ERVISS data folder:
    #   https://github.com/EU-ECDC/Respiratory_viruses_weekly_data/tree/main/data
    # to add a new ERVISS stream, simply add a row here -- the loop below does the rest.
    #   name     : key under data$epi$<name> (used by the downstream code)
    #   file     : exact CSV file name in the ERVISS /data folder
    #   snapshot : local snapshot file under data/ (used when new_from_online==F); streams that serve a
    #              data need of the project one-pager sit in that need's folder, next to its README
    #   schema   : "rates"    -> slim consultation/hospitalisation-rate table
    #              "detailed" -> keep all columns (pathogen / variant / datasource kept)
    erviss_url = "https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/"
    erviss_registry = tribble(
      ~name,                         ~file,                                       ~snapshot,                                                   ~schema,
      "erviss_ili_ari",              "ILIARIRates.csv",                           "primary_care_ili_rates/erviss_iliari.csv",                  "rates",
      "erviss_sari_rates",           "SARIRates.csv",                             "sari_rates/erviss_sari_rates.csv",                          "rates",
      "erviss_typing_sentinel",      "sentinelTestsDetectionsPositivity.csv",     "primary_care_flu_positivity/erviss_detections_sentinel.csv",    "detailed",
      "erviss_typing_nonsentinel",   "nonSentinelTestsDetections.csv",            "primary_care_flu_positivity/erviss_detections_nonsentinel.csv", "detailed",
      "erviss_typing_sari",          "SARITestsDetectionsPositivity.csv",         "sari_flu_positivity/erviss_detections_sari.csv",            "detailed",
      "erviss_flu_type_subtype",     "activityFluTypeSubtype.csv",                "erviss_flu_type_subtype.csv",                               "detailed",
      "erviss_severity_nonsentinel", "nonSentinelSeverity.csv",                   "flu_hospitalisation_indicators/erviss_severity_nonsentinel.csv", "detailed",
      "erviss_sequencing",           "sequencingVolumeDetectablePrevalence.csv",  "erviss_sequencing.csv",                                     "detailed",
      "erviss_variants",             "variants.csv",                              "erviss_variants.csv",                                       "detailed"
    )

    # ---- |-load + standardise every ERVISS stream ----
    epi = list( date_epilist_created = today() )    # initiate the epi list with a creation time-stamp
    for (i in seq_len(nrow(erviss_registry))) {     # i = 1
      reg_i         = erviss_registry[i,]
      snapshot_path = here("data", reg_i$snapshot)
      # fetch online when asked to, or when the local snapshot is missing (self-bootstrapping)
      fetch_online  = new_from_online==T | !file.exists(snapshot_path)
      if (fetch_online) {
        pr=paste0("Loading ERVISS '",reg_i$name,"' from github ... \n"); cat(green(pr))
        raw_i = read_csv(file=paste0(erviss_url, reg_i$file), show_col_types = FALSE)
        raw_i %>% write_csv(file=snapshot_path)     # keep / refresh the local snapshot
      } else {
        pr=paste0("Loading ERVISS '",reg_i$name,"' from disk ... \n"); cat(green(pr))
        raw_i = read_csv(file=snapshot_path, show_col_types = FALSE)
      }
      # add date, recode age, add country code (and slim down if it is a rate table)
      check_required_columns(raw_i, erviss_required_cols(reg_i$schema),
                             paste0("ERVISS '", reg_i$name, "' (", reg_i$file, ")"))  # fail loudly on upstream column changes
      epi[[reg_i$name]] = standardise_erviss(raw_i, schema=reg_i$schema)
    }

    # ---- |-RespiCompass ili_plus (not an ERVISS file, kept alongside the epi streams) ----
    iliplus_snapshot = here("data/data_respicompass_iliplus.csv")
    if (new_from_online==T | !file.exists(iliplus_snapshot)) {
      pr=paste("Loading RespiCompass ili_plus from github ... \n"); cat(green(pr))
      data_respicompass_iliplus = read_csv(file=respicompass_file("target-data/influenza/ili_plus.csv", params), show_col_types = FALSE)
      data_respicompass_iliplus %>% write_csv(file=iliplus_snapshot)
    } else {
      pr=paste("Loading RespiCompass ili_plus from disk ... \n"); cat(green(pr))
      data_respicompass_iliplus = read_csv(file=iliplus_snapshot,show_col_types = FALSE)
    }
    epi$respicompass_iliplus = data_respicompass_iliplus %>%
      mutate(date=ISOweek2date(paste0(yearweek,"-3"))) %>%
      mutate(age=recode_age(age)) %>%
      mutate(countrycode=EU_short(location_name),target="ili_plus") %>%
      select(
        country_short=countrycode,
        date,
        target,
        agegroup=age,
        value=value
      )

    return(epi)
  }

  data$epi = load_or_build("epi", build_epi, regenerate=regenerate)
  return(data)
}

load_data_vax = function(data=list(), params=NULL, regenerate=FALSE, new_from_online=TRUE){

  build_vax = function(){
    # fetch online when asked to, or when the local snapshots are missing (self-bootstrapping)
    fetch_online = new_from_online==T | !file.exists(here("data/vax_flu_scenarios.csv"))
    if (fetch_online) {
      pr=paste("Loading vaccination data from github ... \n"); cat(green(pr))
      data_vax          = read_csv(respicompass_file("auxiliary-data/influenza/vaccination/influenza_vax_scenarios.csv", params), show_col_types = FALSE)
      data_vax_hist     = read_csv(respicompass_file("auxiliary-data/influenza/vaccination/vaccine_coverage_65plus.csv", params), show_col_types = FALSE)
      data_vax_hist_all = read_csv(respicompass_file("auxiliary-data/influenza/vaccination/vaccine_coverage_all.csv",    params), show_col_types = FALSE)
      # keep / refresh the local snapshots
      data_vax          %>% write_csv(file=here("data/vax_flu_scenarios.csv"))
      data_vax_hist     %>% write_csv(file=here("data/vax_flu_history.csv"))
      data_vax_hist_all %>% write_csv(file=here("data/vax_flu_history_all.csv"))
    } else {
      pr=paste("Loading vaccination data from disk ... \n"); cat(green(pr))
      data_vax          = read_csv(file=here("data/vax_flu_scenarios.csv"),   show_col_types = FALSE)
      data_vax_hist     = read_csv(file=here("data/vax_flu_history.csv"),     show_col_types = FALSE)
      data_vax_hist_all = read_csv(file=here("data/vax_flu_history_all.csv"), show_col_types = FALSE)
    }

    # ":" is the Eurostat/ECDC "not available" placeholder -> map to NA before the numeric cast,
    # so the coercion is deliberate (no "NAs introduced by coercion" warning).
    # data_vax_history     = 65+ coverage only; data_vax_history_all = every reported target group.
    list(
      data_vax             = data_vax %>% mutate(vaccine_coverage=vaccine_coverage/100) %>% pivot_wider(names_from = "scenario", values_from = vaccine_coverage),
      data_vax_history     = data_vax_hist     %>% mutate(vaccine_coverage=as.numeric(na_if(vaccine_coverage,":"))/100) %>% mutate(season = str_replace(season, "-", "/")),
      data_vax_history_all = data_vax_hist_all %>% mutate(vaccine_coverage=as.numeric(na_if(vaccine_coverage,":"))/100) %>% mutate(season = str_replace(season, "-", "/"))
    )
  }

  data$vax = load_or_build("vax", build_vax, regenerate=regenerate)
  return(data)
}

# contact data is LOCAL-ONLY: the Prem et al. synthetic contact matrices ship with the repo
# (data/MUestimates_all_locations_{1,2}.xlsx), so there is no online mode here.
load_data_contact = function(data=list(), params=NULL, regenerate=FALSE){

  build_contact = function(){
    # the country list comes from the already-loaded RespiCompass helpers (run that loader first)
    if (is.null(data$helpers_respicompass)) {
      stop("load_data_contact() needs data$helpers_respicompass -- run load_data_helpers_respicompass() first.")
    }
    xlocations = data$helpers_respicompass$iso2_code

    # RespiCompass and the Prem workbooks disagree on some country names; map before the sheet lookup
    # (without this, the mismatch fails silently and the country is dropped from models_in$contacts)
    sheet_overrides = c("Czechia" = "Czech Republic")

    xdata = list()
    for (country_i in xlocations$location_name){
      sheet_i = if (country_i %in% names(sheet_overrides)) sheet_overrides[[country_i]] else country_i
      contacts = 0
      # the matrices are split across two workbooks (each holds ~half the countries); a country is a
      # sheet in exactly one of them, so we try both and keep whichever read succeeds
      try({
        contacts = read_excel(here("data/MUestimates_all_locations_1.xlsx"), sheet = sheet_i, col_names = F, .name_repair = "unique_quiet", skip = 1)
      }, silent = T)
      try({
        contacts = read_excel(here("data/MUestimates_all_locations_2.xlsx"), sheet = sheet_i, col_names = F, .name_repair = "unique_quiet")
      }, silent = T)
      xdata[[country_i]] = contacts
    }
    # declare, rather than silently swallow, the countries with no matrix in either workbook
    # (left at the scalar-0 sentinel; transform_contacts() skips them)
    missing_i = names(xdata)[vapply(xdata, function(x) length(x) == 1, logical(1))]
    if (length(missing_i)) cat(yellow(paste0("load_data_contact: no Prem contact matrix for: ",
                                             paste(missing_i, collapse=", "), "\n")))
    return(xdata)
  }

  data$contact = load_or_build("contact", build_contact, regenerate=regenerate)
  return(data)
}

load_data_helpers_respicompass = function(data=list(), params=NULL, regenerate=FALSE, new_from_online=TRUE){

  build_helpers = function(){
    # fetch online when asked to, or when either snapshot is missing (self-bootstrapping)
    fetch_online = new_from_online==T |
      !file.exists(here("output/respicompass_locations.csv")) |
      !file.exists(here("output/respicompass_weeks.csv"))
    if (fetch_online) {
      pr=paste("Loading RespiCompass helpers from github ... \n"); cat(green(pr))
      xlocations = read_csv(respicompass_file("supporting-files/locations_iso2_codes.csv", params), show_col_types = F)
      xweeks     = read_csv(respicompass_file("supporting-files/iso_weeks.csv",            params), show_col_types = F)
      # keep / refresh the local snapshots
      xlocations %>% write_csv(file=here("output/respicompass_locations.csv"))
      xweeks     %>% write_csv(file=here("output/respicompass_weeks.csv"))
    } else {
      pr=paste("Loading RespiCompass helpers from disk ... \n"); cat(green(pr))
      xlocations = read_csv(here("output/respicompass_locations.csv"), show_col_types = F)
      xweeks     = read_csv(here("output/respicompass_weeks.csv"),     show_col_types = F)
    }

    list(
      iso2_code = xlocations,
      iso_weeks = xweeks
    )
  }

  data$helpers_respicompass = load_or_build("respicompass_helpers", build_helpers, regenerate=regenerate)
  return(data)
}

load_data_demography_ECDC = function(data=list(), params=NULL, regenerate=FALSE, new_from_online=FALSE){

  build_demography_ECDC = function(){
    if (new_from_online==T) {
      # The live ECDC SQL database (out.DM_Population_ByCountryEU) is reachable only inside the
      # ECDC network; its client is not part of this repository. Use the committed snapshot.
      stop("Live ECDC demography DB access is not available here; set params$use_ecdc_db=FALSE ",
           "to use the committed data/population_pyramid.fst snapshot.")
    }
    # ---- read the committed snapshot (works without the ECDC network) ----
    pr=paste("Loading demography data from disk ... \n"); cat(green(pr))
    mdat = read_fst(path=here("data/population_pyramid.fst")) %>% as_tibble()

    list(
      population_pyramid = mdat
    )
  }

  data$demography_ECDC = load_or_build("demography", build_demography_ECDC, regenerate=regenerate)
  return(data)
}

load_data_demography_respicast = function(data=list(), params=NULL, regenerate=FALSE, new_from_online=TRUE){

  build_demography_respicast = function(){
    # fetch online when asked to, or when the snapshot is missing (self-bootstrapping)
    fetch_online = new_from_online==T | !file.exists(here("output/population_pyramid_respicast.csv"))
    if (fetch_online) {
      # RespiCompass ships one population file per country (aggregated + fine age bands); loop and stack.
      if (is.null(data$helpers_respicompass)) {
        stop("load_data_demography_respicast() needs data$helpers_respicompass -- run load_data_helpers_respicompass() first.")
      }
      country_v = data$helpers_respicompass$iso2_code$location_name
      # read one country's population file (aggregated age bands, or fine single-year bands)
      load_one_country = function(country_i, fine=FALSE){
        pr=paste("> Loading respicast pop for:",country_i, if (fine) "(fine)" else "(aggregated)", "... \n"); cat(green(pr))
        relpath = paste0("auxiliary-data/miscellaneous/population/", country_i, if (fine) ".csv" else "_aggr.csv")
        read_csv(respicompass_file(relpath, params), show_col_types = FALSE) %>% mutate(country = country_i)
      }
      pop_df      = map_dfr(country_v, load_one_country, fine=FALSE) %>% select(country, age_group, population)
      pop_fine_df = map_dfr(country_v, load_one_country, fine=TRUE)  %>% select(country, age_group, population)
      # keep / refresh the local snapshots
      pop_df      %>% write_csv(here("output/population_pyramid_respicast.csv"))
      pop_fine_df %>% write_csv(here("output/population_pyramid_fine_respicast.csv"))
    } else {
      pr=paste("Loading respicast demography data from disk ... \n"); cat(green(pr))
      pop_df      = read_csv(here("output/population_pyramid_respicast.csv"),      show_col_types = F)
      pop_fine_df = read_csv(here("output/population_pyramid_fine_respicast.csv"), show_col_types = F)
    }

    list(
      population_pyramid      = pop_df,
      population_pyramid_fine = pop_fine_df
    )
  }

  data$demography_respicast = load_or_build("demography_respicast", build_demography_respicast, regenerate=regenerate)
  return(data)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Mother function: calling the data-loading functions for each data stream ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Order matters: helpers_respicompass produces the country list that the contact and
# respicast-demography streams depend on, so it runs before them.
load_data = function( params=NULL , regenerate=FALSE , new_from_online=TRUE ){

  data = list() # reset data list

  data = load_data_epi( data=data, params=params, regenerate=regenerate, new_from_online=new_from_online )

  data = load_data_vax( data=data, params=params, regenerate=regenerate, new_from_online=new_from_online )

  data = load_data_helpers_respicompass( data=data, params=params, regenerate=regenerate, new_from_online=new_from_online )

  # contact is local-only (Prem matrices shipped in the repo) -> no online fetch
  data = load_data_contact( data=data, params=params, regenerate=regenerate )

  # ECDC demography lives in the ECDC internal SQL DB (only inside the ECDC network) -> dormant
  # by default; set params$use_ecdc_db=TRUE to query it live, otherwise the committed snapshot is used.
  use_ecdc_db = isTRUE(params$use_ecdc_db)
  data = load_data_demography_ECDC( data=data, params=params, regenerate=regenerate, new_from_online=use_ecdc_db )

  data = load_data_demography_respicast( data=data, params=params, regenerate=regenerate, new_from_online=new_from_online )

  return(data)
}
