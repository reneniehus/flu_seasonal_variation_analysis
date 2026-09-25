# joint_model.R -- data assembly, parameter packing, and the fitting strategy for the joint model.
# The model itself is MODEL.md; the log-posterior is joint_model.cpp; this file is the driver.
#
# Fitting strategy, and why. The likelihood is SEPARABLE: only R0_s, the season deviations and the one
# elderly susceptibility appear in more than one country's terms. Everything else is local to a single
# country. A flat finite-difference gradient over ~184 parameters costs 185 full-likelihood
# evaluations; a block sweep costs about 33, because a country's own likelihood is a twelfth of the
# joint one. So we alternate: every country's local block optimised in parallel with the shared block
# held fixed, then the shared block with the locals held fixed, repeated, then one joint polish.
#
# Two guards sit inside that loop, each because the failure it prevents actually happened. Every local
# block is MULTI-STARTED, because a block has a second optimum in which the dispersion collapses and
# the country flat-lines at its baseline instead of fitting a wave. And the FLAT-LINE PROTECTOR
# (jm_flat_check / jm_unflatten) runs before and after the polish, so the returned fit is always
# checked; see those functions for what it detects and what it refuses to do.

suppressMessages({library(Rcpp); library(dplyr)})

# ---- |-settings: the fixed constants and the priors, in one place ----
jm_settings = function(){
  list(
    # fixed, not fitted (MODEL.md)
    gamma_per_day = 0.2777778,     # 3.6-day mean infectious period, from the project notes
    ve_inf = 0.25, ve_ili = 0.20, ve_spread = 0.20,
    rate_per = 1e5, vax_day = 62L, # the 65+ vaccination pulse, 1 October
    # R0 IS FIXED, NOT FITTED (owner, 2026-09-25). Transmissibility and susceptibility enter the rise
    # rate as a product, so one wave identifies only their product; pinning R0 from the literature
    # makes S0 the single sensor of "how easily did this season spread here", and any real
    # season-to-season transmissibility variation is absorbed into S0's season effect by design.
    R0_fixed = 1.5,
    # THE SENSING SWITCH (2026-09-25 model comparison). Which quantity carries the season effect and
    # which the country effect, each "S0" or "R0". The working model is S0/S0. Whatever carries no
    # effect is pinned at its anchor: R0_fixed above, S0_fixed here. Same layout and parameter count
    # in every setting, so the four combinations compare by likelihood alone.
    season_on = "S0", country_on = "S0",
    S0_fixed = 0.75,
    # priors, all on the unconstrained scale the fit works in
    pr_x_sd = 0.5,                            # season effect on logit S0, centred on zero. At S0 ~ 0.8
                                              # one sd is ~ +/-0.09 on S0 itself; wide enough for the
                                              # data to dominate (contraction is reported), tight
                                              # enough to close the S0 -> 1 escape
    pr_S0_mean = qlogis(0.75), pr_S0_sd = 1,  # the COUNTRY level of logit S0 (the mean of the two-way
                                              # decomposition lives here, not in a separate slot)
    # the R0-side counterparts, used only when the switch puts an effect on R0. Both weak: the
    # comparison is by likelihood, and the priors must not be what decides it. r_s at 0.15 on the log
    # scale is the same width the previous fitted-R0_s model used; log R0_c at 0.3 spans 0.8-2.7.
    pr_R0c_mean = log(1.5), pr_R0c_sd = 0.3,
    pr_r_sd = 0.15,
    pr_sigma_mean = 1, pr_sigma_sd = 0.5,     # log2 sigma_eld ~ N(1, .5): centre 2x, 95% band 1.0-4.0x
    pr_delta_sd = 0.5,                        # season observation deviation, constrained to average 0
    pr_off_sd = 1,                            # log2 age reporting offsets
    pr_c_mean = log(0.05), pr_c_sd = 2,       # weak; closes the 'no epidemic' trap at c -> 0
    # The data's own week-to-week scatter around a 3-week mean implies phi ~ 3.7 (median over 257
    # country-season-age series), so that is the centre. The sd is deliberately WIDE: without the
    # filter the mean is rigid, so phi must be free to absorb model misfit and REPORT it. A tight
    # prior would force phi up and push the misfit into R0 and S0 instead, which is worse. The gap
    # between the fitted phi and 3.7 is therefore the model-adequacy diagnostic (jm_adequacy).
    pr_phi_mean = log(4), pr_phi_sd = 1.0,
    pr_b_mean = log(1), pr_b_sd = 2,          # NEW: the pilot left b unpenalised, which left a
                                              # single-season baseline effectively free
    pr_I0_mean = log(10^-6.5), pr_I0_sd = 3
  )
}

jm_load_cpp = function(dir = "output/joint_model/cpp_cache"){
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(here::here("code/07_joint_model/joint_model.cpp"), cacheDir = dir, verbose = FALSE)
  invisible(TRUE)
}

# ---- |-assemble every country into the flat structure the C++ expects ----
# Reuses the compartmental pilot's data layer (build_comp_data): the panel stitch, the 3-group
# collapse of the contact matrix, the populations and the 65+ coverage are already settled there.
# Counts are ROUNDED here: the panel holds rates, the count scale is a device (MODEL.md), and the
# negative binomial wants integers.
# min_seasons = 5L is the OWNER'S DECISION of 2026-09-12 ("5 keeps Spain in"), which is what gives the
# documented design of 12 countries / 86 country-seasons / 184 parameters. It was left at 6 here while
# only run_joint_model.R passed the override, so any other caller reproducing "the model" as the
# documents describe it silently got an 11-country / 81 / 173 design instead, announced by one buried
# line of verbose output. The default now IS the decision, and a test pins the resulting design.
jm_build_data = function(countries, models_in, demo, set = jm_settings(), min_seasons = 5L,
                         exclude_ambiguous_positivity = TRUE, ambiguous_min_unexplained = 1L,
                         verbose = TRUE){
  # The data layer is configured by comp_model_settings() (line below) while `set` supplies the same
  # two constants to the fitted object. They agree today only by coincidence -- nothing tied them --
  # so a change to either file alone would silently put d$y on one basis and d$rate_per on another.
  cms = comp_model_settings()
  if (!isTRUE(all.equal(set$rate_per, cms$rate_per)))
    stop(sprintf("rate_per disagrees between jm_settings() (%g) and comp_model_settings() (%g); d$y is built with the latter, so the two must match",
                 set$rate_per, cms$rate_per))
  if (anyDuplicated(countries)){                     # cds[[cc]] assigns by NAME, so a repeated code
    dup = unique(countries[duplicated(countries)])   # would overwrite and the design silently shrink
    if (verbose) cat("  dropping duplicate country code(s):", paste(dup, collapse = ", "), "\n")
    countries = unique(countries)
  }
  # PROVISIONAL EXCLUSION, PENDING CONFIRMATION BY SURVEILLANCE COLLEAGUES.
  # See documentation/to_confirm_with_surveillance.md (question 1) and erviss_encoding_ambiguous().
  #
  # ERVISS encodes a zero-detection week two ways -- "detections = 0" explicitly, or the detections
  # row absent with tests > 0 -- and our re-derived positivity turns the second into NA, which the
  # stitch deletes; a trailing run shortens the season instead of leaving holes.
  #
  # WHAT IS EXCLUDED IS DECIDED ON THE EVIDENCE, not on a count of affected weeks (owner,
  # 2026-09-16). For each affected week, erviss_encoding_ambiguous computes the probability that it
  # really was zero, from that week's test count and the local positivity of its published
  # neighbours. A country-season is excluded only if it contains weeks that CANNOT plausibly have
  # been zero -- implausible under the local positivity, or unknowable because nothing was published
  # nearby. That distinction is the whole point: it excludes CZ 2024/2025, where one week expected
  # ~4 detections (p = 0.014) and sixteen more have no published neighbour at all, while keeping
  # PL 2024/2025, whose single affected week sits among neighbours with 0 detections over 114 tests
  # (p = 1.000) and is almost certainly a genuine zero. A count threshold would have separated those
  # two only by luck.
  #
  # exclude_ambiguous_positivity = FALSE fits everything; ambiguous_min_unexplained raises how many
  # not-plausibly-zero weeks a country-season may carry before it is dropped.
  ambig = NULL
  if (exclude_ambiguous_positivity){
    ambig = tryCatch(erviss_encoding_ambiguous(models_in), error = function(e) NULL)
    if (is.null(ambig) && verbose)
      cat("  NOTE: could not evaluate the positivity-encoding flag; no season excluded for it\n")
    else ambig = ambig[ambig$n_not_plausibly_zero >= ambiguous_min_unexplained, , drop = FALSE]
  }
  excluded = data.frame(country = character(0), season = character(0), n_ambiguous = integer(0),
                        n_not_plausibly_zero = integer(0), n_no_neighbour = integer(0),
                        min_p_zero = numeric(0), stringsAsFactors = FALSE)

  cds = list()
  for (cc in countries){
    # Report WHY a country was dropped. "(no data)" was printed for every failure mode, including the
    # ones it cannot be -- a missing external CSV, a bad working directory, an unknown country code --
    # so a country with a full panel could vanish from the design while the log blamed its data.
    cd = tryCatch(build_comp_data(cc, models_in, demo, cms), error = function(e) conditionMessage(e))
    if (is.character(cd)){ if (verbose) cat("  skip", cc, "-- build_comp_data failed:", cd, "\n"); next }
    if (is.null(cd)){ if (verbose) cat("  skip", cc, "(no data)\n"); next }
    # PRUNE unusable seasons before anything counts them. A season with no finite observation
    # contributes nothing to the likelihood but would still claim a seed slot that no country-season
    # references (a dead parameter) and would break the data-driven starting values, since the
    # quantile of an empty set is NA. This has to happen before n_local and off_country are computed.
    keep = which(vapply(cd$y, function(m) any(is.finite(m)), logical(1)))
    # ... and, in the SAME prune, the seasons whose positivity encoding is ambiguous. Both have to
    # happen before n_local and off_country are computed, or a dropped season leaves an orphaned seed
    # slot that no country-season references.
    if (!is.null(ambig) && nrow(ambig)){
      sub = ambig[ambig$country_short == cc, , drop = FALSE]
      if (nrow(sub)){
        drop_i = which(cd$seasons %in% sub$season)
        if (length(drop_i)){
          j = match(cd$seasons[drop_i], sub$season)
          excluded = rbind(excluded, data.frame(
            country = cc, season = cd$seasons[drop_i],
            n_ambiguous = sub$n_ambiguous[j],
            n_not_plausibly_zero = sub$n_not_plausibly_zero[j],
            n_no_neighbour = sub$n_no_neighbour[j],
            min_p_zero = sub$min_p_zero[j],
            stringsAsFactors = FALSE))
          if (verbose) cat("  ", cc, ": excluding", length(drop_i),
                           "season(s) for ambiguous positivity encoding:",
                           paste(cd$seasons[drop_i], collapse = ", "), "\n")
          keep = setdiff(keep, drop_i)
        }
      }
    }
    if (length(keep) < length(cd$seasons)){
      if (verbose) cat("  ", cc, ": keeping", length(keep), "of", length(cd$seasons), "season(s)\n")
      cd$seasons = cd$seasons[keep]; cd$y = cd$y[keep]; cd$rates = cd$rates[keep]
      cd$vax = cd$vax[keep, , drop = FALSE]
      cd$source_by_season = cd$source_by_season[keep]
      if (!is.null(cd$n_weeks)) cd$n_weeks = cd$n_weeks[keep]
    }
    if (length(cd$seasons) < min_seasons){ if (verbose) cat("  skip", cc, sprintf("(%d usable seasons < %d)\n", length(cd$seasons), min_seasons)); next }
    cds[[cc]] = cd
  }
  stopifnot(length(cds) > 0)
  seasons_all = sort(unique(unlist(lapply(cds, `[[`, "seasons"))))
  C = length(cds); S = length(seasons_all)

  # per country: the data sources it actually has. The pilot created a slot for both sources
  # unconditionally, which left an inert parameter wherever a source was absent (ES, NO) and made the
  # Hessian exactly singular. Here the slots follow the data.
  srcs = lapply(cds, function(cd) sort(unique(unname(cd$source_by_season))))
  n_src = vapply(srcs, length, integer(1))
  n_cs_of_country = vapply(cds, function(cd) length(cd$seasons), integer(1))
  n_local = 5L + n_src + n_cs_of_country                 # S0, c, 2 offsets, phi, baselines, seeds
  n_shared = 2L * S - 1L                                 # S-1 free x, S-1 free delta, 1 log2 sigma
  off_country = as.integer(cumsum(c(n_shared, head(n_local, -1))))   # 0-based starts
  n_par = as.integer(n_shared + sum(n_local))

  y = list(); cs_country = integer(0); cs_season = integer(0); cs_src = integer(0)
  cs_pos = integer(0); vax_eld = numeric(0); lgamma_y1 = numeric(0)
  cs_label = character(0); n_weeks = integer(0)
  for (ic in seq_along(cds)){
    cd = cds[[ic]]
    for (k in seq_along(cd$seasons)){
      ym = round(cd$y[[k]]); ym[ym < 0] = 0
      y[[length(y) + 1]] = ym
      cs_country = c(cs_country, ic - 1L)
      cs_season  = c(cs_season, match(cd$seasons[k], seasons_all) - 1L)
      cs_src     = c(cs_src, match(unname(cd$source_by_season[k]), srcs[[ic]]) - 1L)
      cs_pos     = c(cs_pos, k - 1L)
      vax_eld    = c(vax_eld, cd$vax$coverage[k])
      lgamma_y1  = c(lgamma_y1, sum(lgamma(ym[is.finite(ym)] + 1)))
      cs_label   = c(cs_label, paste(names(cds)[ic], cd$seasons[k]))
      n_weeks    = c(n_weeks, nrow(ym))
    }
  }
  cs_of_country = lapply(seq_along(cds), function(ic) as.integer(which(cs_country == ic - 1L) - 1L))

  # FAIL LOUDLY IN R RATHER THAN ABORT IN C++. cs_season and cs_src come from match() - 1L, which
  # yields NA_integer_ for a label it cannot find -- and srcs is built with sort(unique(...)), which
  # DROPS NA, so a single unresolved source label produces exactly that. The C++ uses these as raw
  # indices (logb[cs_src[i]], sh.R0[cs_season[i]]), so NA_INTEGER = INT_MIN reads far out of bounds and
  # the R process dies with "an irrecoverable exception occurred", taking the whole fit with it and
  # reporting nothing about the cause. Unreachable from the committed panel (its source column has no
  # NA), so this guards a future panel build rather than today's.
  bad = c(cs_country = anyNA(cs_country), cs_season = anyNA(cs_season), cs_src = anyNA(cs_src),
          cs_pos = anyNA(cs_pos), vax_eld = anyNA(vax_eld), lgamma_y1 = anyNA(lgamma_y1))
  if (any(bad))
    stop("unresolved label(s) in the assembled design: ", paste(names(bad)[bad], collapse = ", "),
         ". A country-season carries a season or source label that is not in the design's own list; ",
         "fix the panel rather than fitting, because the C++ would read these as out-of-range indices.")

  d = c(list(
    n_country = C, n_season = S, n_cs = length(y), n_par = n_par,
    countries = names(cds), seasons = seasons_all, sources = srcs,
    Cn = lapply(cds, `[[`, "Cn"), N = lapply(cds, `[[`, "N"), groups = cds[[1]]$groups,
    # PROVENANCE OF THE MIXING PATTERN. build_comp_data records whether a country uses its own
    # Prem-derived contact matrix or the EU average; that was being discarded here, so nothing
    # downstream could disclose it. Norway has no matrix of its own and gets the EU average, which
    # matters because the age reporting offsets are the parameters most sensitive to the mixing
    # pattern (MODEL.md).
    contact_source = vapply(cds, function(cd)
      if (is.null(cd$contact_source)) NA_character_ else as.character(cd$contact_source), character(1)),
    n_src = as.integer(n_src), n_cs_of_country = as.integer(n_cs_of_country),
    n_local = as.integer(n_local), off_country = off_country,
    cs_of_country = cs_of_country, cs_country = cs_country, cs_season = cs_season,
    cs_src = cs_src, cs_pos = cs_pos, cs_label = cs_label, n_weeks = as.integer(n_weeks),
    vax_eld = vax_eld, lgamma_y1 = lgamma_y1, y = y,
    rates = unlist(lapply(cds, `[[`, "rates"), recursive = FALSE),
    gamma = set$gamma_per_day, ve_inf = set$ve_inf, ve_ili = set$ve_ili,
    ve_spread = set$ve_spread, rate_per = set$rate_per, vax_day = set$vax_day,
    R0_fixed = set$R0_fixed, S0_fixed = set$S0_fixed,
    season_on = set$season_on, country_on = set$country_on,
    # The DYNAMICS horizon, distinct from the observation windows. Each country-season is observed for
    # however long its surveillance series runs (33 to 53 weeks here), but the attack rate has to mean
    # the same thing in every cell to be comparable across them and against cohort evidence, so the
    # epidemic is integrated to a full season everywhere. Only the observed weeks enter the
    # likelihood, so this changes no fitted value -- see simulate_season in the C++.
    attack_weeks = max(53L, max(as.integer(n_weeks))),
    # what was dropped for the ambiguous positivity encoding, so the exclusion is visible in the
    # fitted object rather than only in a build log that scrolls away
    excluded_ambiguous = excluded
  ), set[grep("^pr_", names(set))])
  if (verbose){
    cat(sprintf("%d countries, %d seasons, %d country-seasons, %d parameters (%d shared, %d local)\n",
                C, S, length(y), n_par, n_shared, sum(n_local)))
    if (nrow(excluded)){
      cat(sprintf("EXCLUDED for ambiguous positivity encoding (%d country-season(s), PROVISIONAL --\n",
                  nrow(excluded)))
      cat("  awaiting confirmation from surveillance colleagues, see erviss_encoding_ambiguous()):\n")
      for (k in seq_len(nrow(excluded)))
        cat(sprintf("    %s %s: %d affected week(s), of which %d cannot plausibly have been zero (%d with no published neighbour at all)\n",
                    excluded$country[k], excluded$season[k], excluded$n_ambiguous[k],
                    excluded$n_not_plausibly_zero[k], excluded$n_no_neighbour[k]))
    }
  }
  d
}

# ---- |-the parameter vector: names, index blocks, starting values ----
jm_par_names = function(d){
  S = d$n_season
  sea = if (identical(d$season_on, "R0")) "r_" else "x_"          # r_s on log R0, x_s on logit S0
  cty = if (identical(d$country_on, "R0")) ":log_R0" else ":logit_S0"
  nm = c(paste0(sea, head(d$seasons, S - 1)), paste0("delta_", head(d$seasons, S - 1)), "log2_sigma_eld")
  for (ic in seq_len(d$n_country)){
    cc = d$countries[ic]
    nm = c(nm, paste0(cc, c(cty, ":log_c", ":off_young", ":off_eld", ":log_phi")),
           paste0(cc, ":log_b_", d$sources[[ic]]),
           paste0(cc, ":log_I0_", d$seasons[d$cs_season[d$cs_of_country[[ic]] + 1L] + 1L]))
  }
  nm
}
jm_blocks = function(d){
  S = d$n_season
  list(shared = seq_len(2L * S - 1L),
       local = lapply(seq_len(d$n_country), function(ic) d$off_country[ic] + seq_len(d$n_local[ic])))
}

jm_theta0 = function(d, set = jm_settings()){
  th = numeric(d$n_par)
  S = d$n_season
  th[seq_len(S - 1)] = 0                                # no season effect on susceptibility
  th[S - 1L + seq_len(S - 1)] = 0                       # no season effect on visibility
  th[2L * S - 1L] = set$pr_sigma_mean                   # elderly susceptibility at 2x
  r0 = set$gamma_per_day * (set$R0_fixed * 0.8 - 1)     # a plausible early growth rate, per day
  for (ic in seq_len(d$n_country)){
    base = d$off_country[ic]
    ics = d$cs_of_country[[ic]] + 1L
    rate_all = unlist(lapply(ics, function(i) as.numeric(d$rates[[i]])))
    peak = max(c(rate_all, 0), na.rm = TRUE)
    # na.rm on the OUTER max as well: the quantile of an empty positive set is NA, and max(NA, 1e-3)
    # is NA, which would put NA into log_b and leave that whole block silently un-optimised
    floor_r = max(c(as.numeric(quantile(rate_all[rate_all > 0], 0.10, na.rm = TRUE)), 1e-3), na.rm = TRUE)
    c0 = min(max(peak / d$rate_per / 0.02, 1e-4), 1)
    th[base + 1] = if (identical(d$country_on, "R0")) log(set$R0_fixed) else qlogis(0.80)
    th[base + 2] = log(c0)
    th[base + 3] = 0; th[base + 4] = 0
    th[base + 5] = set$pr_phi_mean
    th[base + 5 + seq_len(d$n_src[ic])] = log(floor_r)
    # The seed start comes from the OBSERVED onset: at growth rate r0 the default seed reaches onset
    # around week 12, so a season whose wave crosses 10% of its peak in week o is shifted by (o - 12)
    # weeks. Three things make it ROBUST, because the naive version was set by a single early spike
    # (DK 2025/2026 had a week-2 value of 19.5 against a peak of 128, which put onset at 2 and the
    # seed start 8.6 log units from its fitted value): the off-season floor is subtracted first, the
    # crossing is tested on a 3-week mean so one week cannot trigger it, and the result is clamped to
    # a plausible window. Measured: this start and a flat one at the prior mean converge to the same
    # optimum, so it buys speed and robustness rather than the answer.
    onset = vapply(ics, function(i){
      m = d$rates[[i]]; tot = rowSums(sweep(m, 2, d$N[[ic]] / sum(d$N[[ic]]), "*"), na.rm = TRUE)
      tot[rowSums(is.finite(m)) == 0] = NA
      if (!any(is.finite(tot))) return(12)
      z = tot - as.numeric(quantile(tot, 0.10, na.rm = TRUE))     # above the off-season floor
      n = length(z); sm = rep(NA_real_, n)                        # 3-week mean, NA-tolerant
      for (k in seq_len(n)){ w = z[max(1, k - 1):min(n, k + 1)]; if (all(is.finite(w))) sm[k] = mean(w) }
      pk = max(sm, na.rm = TRUE)
      if (!is.finite(pk) || pk <= 0) return(12)
      w = which(sm >= 0.1 * pk)[1]
      if (is.na(w)) 12 else min(max(w, 4), 34) }, numeric(1))
    th[base + 5 + d$n_src[ic] + seq_along(ics)] = set$pr_I0_mean - r0 * 7 * (onset - 12)
  }
  names(th) = jm_par_names(d)
  th
}

# ---- |-the dispersion the DATA's own scatter implies, per country ----
# The weekly series scatters around its own 3-week moving average by some amount that no mean could
# ever beat, so that scatter is a floor on measurement noise. Read as a negative-binomial dispersion
# it gives phi = 1/cv^2. Used in three places and therefore defined once: the model-adequacy
# diagnostic compares the fitted phi against it, the multi-start uses it as an alternative start, and
# the flat-line protector pins phi there while it escapes.
jm_phi_data = function(d){
  ma3 = function(v){ n = length(v); o = rep(NA_real_, n)
    if (n >= 3) for (i in 2:(n - 1)){ w = v[(i - 1):(i + 1)]; if (all(is.finite(w))) o[i] = mean(w) }; o }
  vapply(seq_len(d$n_country), function(ic){
    cvs = unlist(lapply(d$cs_of_country[[ic]] + 1L, function(i){
      y = d$y[[i]]; vapply(seq_len(ncol(y)), function(a){
        v = y[, a]; m = ma3(v); ok = is.finite(v) & is.finite(m) & m > 20
        if (sum(ok) < 5) NA_real_ else sd(v[ok] / m[ok]) / sqrt(2/3) }, numeric(1)) }))
    cv = median(cvs, na.rm = TRUE)
    if (!is.finite(cv) || cv <= 0) exp(d$pr_phi_mean) else 1 / cv^2
  }, numeric(1))
}

# ---- |-the FLAT-LINE detector ----
# The failure this guards against: a country's local block has a second optimum in which the
# dispersion collapses, the negative binomial becomes so diffuse that every curve fits about equally
# well, nothing pushes the model to place a wave, and the country sits at its off-season baseline for
# all eight seasons. It is a LOCAL optimum -- the first joint fit lost the Netherlands to it and a
# harder search later found a solution 406 nats better -- so it is an optimiser failure, not a
# statement about the data.
# Detection is MECHANISTIC rather than goodness-of-fit based, because a country can fit badly for
# honest reasons but cannot have an epidemic that never happened:
#   attack_max     the largest population-weighted attack rate over that country's seasons. The model
#                  cannot produce a few per cent at any plausible R0 and S0, so a value near zero
#                  means no epidemic was fitted at all.
#   peak_epi_frac  how much of the fitted peak is epidemic rather than baseline, max(mu) vs min(mu)
#                  over the season. A flat line is all baseline.
# phi_ratio (fitted dispersion over the data-implied one) is reported as CONTEXT, not as a trigger:
# a low value is the mechanism of the failure but is also seen in countries that fit fine.
# Aggregating over seasons HIDES a minority of flat ones: the seed is per country-season while the
# dispersion is per country, so a single wave can flat-line alone, and neither a max over attack rates
# nor a median over epidemic fractions will see it (measured: 3 of the Netherlands' 6 seasons flat left
# the country unflagged while costing 3084 nats). So the per-season minima are reported, and a PARTIAL
# flat line triggers too -- but only when the country demonstrably CAN have epidemics, i.e. at least
# two of its other seasons are healthy. Without that condition a design whose attack rates are
# genuinely small everywhere, as the driver-recovery truths can be, would false-positive.
jm_flat_check = function(th, d, attack_min = 0.03, epi_frac_min = 0.15){
  f = jm_fitted_cpp(th, d); p = jm_unpack(th, d); phid = jm_phi_data(d)
  rows = lapply(seq_len(d$n_country), function(ic){
    ics = d$cs_of_country[[ic]] + 1L; N = d$N[[ic]]; w = N / sum(N)
    attack = vapply(ics, function(i) sum(f$attack[i, ] * w), numeric(1))
    epi = vapply(ics, function(i){ m = f$mu[[i]]; tot = rowSums(m)
      if (!any(is.finite(tot)) || max(tot) <= 0) 0 else (max(tot) - min(tot)) / max(tot) }, numeric(1))
    cors = vapply(ics, function(i){ y = as.numeric(d$y[[i]]); m = as.numeric(f$mu[[i]])
      ok = is.finite(y) & is.finite(m); if (sum(ok) > 3) suppressWarnings(cor(y[ok], m[ok])) else NA_real_ }, numeric(1))
    data.frame(country = d$countries[ic], attack_max = max(attack), attack_med = median(attack),
               attack_min_season = min(attack), n_seasons = length(attack),
               n_seasons_flat = sum(attack < attack_min), n_seasons_healthy = sum(attack > 4 * attack_min),
               peak_epi_frac = median(epi), peak_epi_frac_min = min(epi),
               cor_med = median(cors, na.rm = TRUE),
               phi = p$country[[ic]]$phi, phi_data = phid[ic],
               phi_ratio = p$country[[ic]]$phi / phid[ic], stringsAsFactors = FALSE)
  })
  out = do.call(rbind, rows)
  out$flat_whole = out$attack_max < attack_min | out$peak_epi_frac < epi_frac_min
  out$flat_partial = !out$flat_whole & out$n_seasons_flat > 0 & out$n_seasons_healthy >= 2
  out$flat = out$flat_whole | out$flat_partial
  out
}

# ---- |-the FLAT-LINE PROTECTOR ----
# For every country the detector flags, re-optimise that country's local block from starts designed to
# make the flat mode expensive, and adopt a rescue only if it IMPROVES the objective. That keeps this
# an optimisation device rather than a way of overriding the likelihood.
# The ladder, cheapest first. Each rung first optimises the block with the dispersion HELD at the
# value the data's own scatter implies -- which makes ignoring a wave costly, because a tight
# dispersion cannot explain a missed peak -- and then releases it and re-optimises the whole block:
#   1. the fully data-driven start for that country (jm_theta0's own local block)
#   2. the same with the seed ten times larger (earlier arrival)
#   3. the same with the seed ten times smaller (later arrival)
#   4. the current values with only the dispersion reset
# If after all of them the best objective is STILL flat, the flat solution genuinely wins on the
# likelihood. That is then reported as unresolved rather than papered over: it means that country's
# series carries no identifiable wave, which is a finding about the data and must not be hidden by
# forcing a wave the likelihood does not want.
jm_unflatten = function(th, d, maxit = 400L, cores = max(1L, parallel::detectCores() - 1L),
                        attack_min = 0.03, epi_frac_min = 0.15, verbose = TRUE){
  chk = jm_flat_check(th, d, attack_min, epi_frac_min)
  bad = which(chk$flat)
  report = data.frame(country = character(0), action = character(0),
                      negll_before = numeric(0), negll_after = numeric(0),
                      attack_before = numeric(0), attack_after = numeric(0), stringsAsFactors = FALSE)
  if (!length(bad)){ if (verbose) cat("flat-line protector: no country flagged\n"); return(list(theta = th, report = report, check = chk, n_fixed = 0L, n_unresolved = 0L)) }
  if (verbose) cat(sprintf("flat-line protector: %d country(ies) flagged (%s)\n", length(bad),
                           paste(chk$country[bad], collapse = ", ")))
  bl = jm_blocks(d); th0_data = jm_theta0(d); phid = jm_phi_data(d)
  res = parallel::mclapply(bad, function(ic){
    idx = bl$local[[ic]]; i_phi = 5L
    i_seed = 5L + d$n_src[ic] + seq_len(d$n_cs_of_country[ic])
    f = function(x){ t2 = th; t2[idx] = x; jm_country_negll_cpp(t2, d, ic - 1L) }
    v0 = f(th[idx])
    starts = list(th0_data[idx], th0_data[idx], th0_data[idx], th[idx])
    starts[[2]][i_seed] = starts[[2]][i_seed] + log(10)
    starts[[3]][i_seed] = starts[[3]][i_seed] - log(10)
    best = list(par = th[idx], value = v0, conv = NA_integer_)
    for (st in starts){
      st[i_phi] = log(phid[ic])
      # (a) hold the dispersion, optimise everything else in the block
      free = setdiff(seq_along(idx), i_phi)
      g = function(z){ x = st; x[free] = z; f(x) }
      o1 = tryCatch(optim(st[free], g, method = "BFGS", control = list(maxit = maxit, reltol = 1e-10)),
                    error = function(e) NULL)
      st2 = st; if (!is.null(o1)) st2[free] = o1$par
      # (b) release it and optimise the whole block
      o2 = tryCatch(optim(st2, f, method = "BFGS", control = list(maxit = maxit, reltol = 1e-10)),
                    error = function(e) NULL)
      if (!is.null(o2) && is.finite(o2$value) && o2$value < best$value)
        best = list(par = o2$par, value = o2$value, conv = o2$convergence)
    }
    list(ic = ic, par = best$par, value = best$value, before = v0, conv = best$conv)
  }, mc.cores = min(cores, length(bad)))
  # mclapply returns NULL for a KILLED worker and a "try-error" object for one that raised; the
  # latter would abort the whole fit at r$value with "$ operator is invalid for atomic vectors".
  dead = vapply(res, function(r) is.null(r) || inherits(r, "try-error") || !is.list(r), logical(1))
  if (any(dead)){
    warning(sprintf("flat-line protector: the worker failed for %s, so those countries were NOT rescued",
                    paste(chk$country[bad[dead]], collapse = ", ")), call. = FALSE)
    res = res[!dead]
  }
  for (r in res) if (r$value < r$before) th[bl$local[[r$ic]]] = r$par
  chk2 = jm_flat_check(th, d, attack_min, epi_frac_min)
  for (r in res){
    ic = r$ic
    report = rbind(report, data.frame(country = chk$country[ic],
      action = if (!chk2$flat[ic]) "rescued" else if (r$value < r$before) "improved but still flat" else "unresolved",
      negll_before = r$before, negll_after = min(r$value, r$before),
      attack_before = chk$attack_max[ic], attack_after = chk2$attack_max[ic],
      conv = if (is.null(r$conv)) NA_integer_ else as.integer(r$conv), stringsAsFactors = FALSE))
  }
  n_fixed = sum(report$action == "rescued")
  n_un = sum(chk2$flat)     # from the RE-CHECK: a missing report row must not be able to hide one
  if (verbose){
    print(report, row.names = FALSE, digits = 4)
    cat(sprintf("flat-line protector: %d rescued, %d unresolved\n", n_fixed, n_un))
    if (n_un > 0) cat("  UNRESOLVED means the flat solution wins on the likelihood: that country's series\n",
                      " carries no identifiable wave. Report it; do not force one.\n")
  }
  list(theta = th, report = report, check = chk2, n_fixed = n_fixed, n_unresolved = n_un)
}

# ---- |-fit by block coordinate descent, then one joint polish ----
jm_fit = function(d, theta0 = jm_theta0(d), max_sweeps = 15L, tol = 0.05, maxit_local = 400L,
                  maxit_shared = 400L, polish_maxit = 300L, cores = max(1L, parallel::detectCores() - 1L),
                  attack_min = 0.03, epi_frac_min = 0.15, verbose = TRUE){
  bl = jm_blocks(d); th = theta0
  os = NULL; sw = 0L              # so the returned convergence fields exist even at max_sweeps = 0
  hit_tol = FALSE                 # did the sweep loop reach its tolerance, or run out of sweeps?
  last_gain = NA_real_; prev_gain = NA_real_
  phid = jm_phi_data(d)          # each country's own data-implied dispersion, for the extra starts
  t_start = Sys.time()
  obj = jm_negll_cpp(th, d)
  trace = data.frame(sweep = 0L, step = "start", negll = obj, seconds = 0)
  if (verbose) cat(sprintf("start: negll %.2f\n", obj))
  conv_local = integer(0)
  for (sw in seq_len(max_sweeps)){
    prev = obj
    # --- every country's local block, with the shared block fixed. Independent, so run in parallel.
    # MULTI-START, and why it is necessary. Each country's block has a second, WRONG optimum: let the
    # dispersion collapse, and the negative binomial becomes so diffuse that every curve fits equally
    # well, so there is no pressure to place a wave at all and the country flat-lines at its baseline.
    # Which country falls into it depends only on where BFGS starts -- tightening the dispersion prior
    # moved the failure from the Netherlands to Poland rather than removing it. So each block is
    # started from three points and the best is kept. The extra starts pin the dispersion at the value
    # the data's own scatter implies, which makes ignoring a wave expensive and pulls the block into
    # the wave-fitting mode; the dispersion is free again on the next sweep.
    fits = parallel::mclapply(seq_len(d$n_country), function(ic){
      idx = bl$local[[ic]]
      f = function(x){ t2 = th; t2[idx] = x; jm_country_negll_cpp(t2, d, ic - 1L) }
      i_phi = 5L                                       # position of log_phi within the local block
      i_seed = 5L + d$n_src[ic] + seq_len(d$n_cs_of_country[ic])
      starts = list(th[idx])                           # (1) warm: where the last sweep left it
      s2 = th[idx]; s2[i_phi] = log(phid[ic]); starts[[2]] = s2           # (2) dispersion at THIS country's scatter
      s3 = s2; s3[i_seed] = s3[i_seed] + log(10); starts[[3]] = s3       # (3) ... and an earlier arrival
      # seed with the value we ALREADY have, so a block can never be written back worse than it came
      # in: if the warm start errored and a cold start succeeded at a worse value, the old code
      # adopted the worse one and the objective stopped being monotone
      best = list(par = th[idx], value = f(th[idx]), conv = 0L)
      for (st in starts){
        o = tryCatch(optim(st, f, method = "BFGS", control = list(maxit = maxit_local, reltol = 1e-10)),
                     error = function(e) NULL)
        if (!is.null(o) && is.finite(o$value) && o$value < best$value)
          best = list(par = o$par, value = o$value, conv = o$convergence)
      }
      best
    }, mc.cores = min(cores, d$n_country))
    bad_w = vapply(fits, function(f) is.null(f) || inherits(f, "try-error") || !is.list(f), logical(1))
    if (any(bad_w)) warning(sprintf("sweep %d: the worker failed for %s, so those blocks were left as they were",
                                    sw, paste(d$countries[bad_w], collapse = ", ")), call. = FALSE)
    for (ic in seq_len(d$n_country)) if (!bad_w[ic]) th[bl$local[[ic]]] = fits[[ic]]$par
    conv_local = vapply(seq_len(d$n_country), function(ic) if (bad_w[ic]) 99L else as.integer(fits[[ic]]$conv), integer(1))
    obj_l = jm_negll_cpp(th, d)
    trace = rbind(trace, data.frame(sweep = sw, step = "local", negll = obj_l,
                                    seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
    # --- the shared block, with every local block fixed
    idx = bl$shared
    fs = function(x){ t2 = th; t2[idx] = x; jm_negll_cpp(t2, d) }
    os = tryCatch(optim(th[idx], fs, method = "BFGS", control = list(maxit = maxit_shared, reltol = 1e-10)),
                  error = function(e) NULL)
    if (!is.null(os)) th[idx] = os$par
    obj = jm_negll_cpp(th, d)
    trace = rbind(trace, data.frame(sweep = sw, step = "shared", negll = obj,
                                    seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
    if (verbose) cat(sprintf("sweep %2d: after locals %.2f, after shared %.2f  (gain %.2f)\n",
                             sw, obj_l, obj, prev - obj))
    prev_gain = last_gain; last_gain = prev - obj
    if (last_gain < tol){ hit_tol = TRUE; break }
  }
  # --- the flat-line protector, before the polish: a flat country would otherwise be polished into
  # a locally optimal flat solution and look converged
  prot1 = jm_unflatten(th, d, cores = cores, attack_min = attack_min, epi_frac_min = epi_frac_min,
                       verbose = verbose)
  th = prot1$theta; obj_p = jm_negll_cpp(th, d)
  if (obj_p < obj && verbose) cat(sprintf("  after rescue: negll %.2f (gained %.2f)\n", obj_p, obj - obj_p))
  obj = obj_p
  # --- one joint polish: all parameters together, to clear any residual cross-block curvature
  op = tryCatch(optim(th, function(x) jm_negll_cpp(x, d), method = "BFGS",
                      control = list(maxit = polish_maxit, reltol = 1e-11)), error = function(e) NULL)
  polish_gain = NA_real_
  if (!is.null(op) && is.finite(op$value) && op$value < obj){ polish_gain = obj - op$value; th = op$par; obj = op$value }
  names(th) = jm_par_names(d)
  # --- and once more after the polish, so the RETURNED fit is guaranteed checked
  # ADOPT UNCONDITIONALLY. jm_unflatten's theta is provably never worse -- it adopts a country's block
  # only when that block's objective improves -- so gating on n_fixed (which counts only rescues that
  # CLEARED the threshold) would discard a strictly better vector whenever a country improved without
  # clearing it. Measured: that path threw away a 5648-nat improvement, and because the reported check
  # was computed on the returned theta while the discarded one was kept, fit$flat then described a
  # different parameter vector from fit$theta.
  prot2 = jm_unflatten(th, d, cores = cores, attack_min = attack_min, epi_frac_min = epi_frac_min,
                       verbose = verbose)
  th = prot2$theta; obj = jm_negll_cpp(th, d); names(th) = jm_par_names(d)
  flat_final = jm_flat_check(th, d, attack_min, epi_frac_min)   # on the theta we actually return
  trace = rbind(trace, data.frame(sweep = 99L, step = "polish", negll = obj,
                                  seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  secs = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (verbose) cat(sprintf("done in %.1f s: negll %.2f (polish gained %.2f), loglik %.2f\n",
                           secs, obj, polish_gain, jm_loglik_cpp(th, d)))
  # --- what "converged" means here, because the obvious definition is misleading.
  # The fit has THREE stages: each country's local block, the shared block, then one joint polish. The
  # sweep loop over the first two is expected to run out of sweeps: its gain decays geometrically and
  # never reaches a 0.05-nat tolerance in 15 sweeps (measured: the 15th sweep still gains 2.3 nats,
  # extrapolating to a ~24-nat tail). The JOINT POLISH is what clears that tail -- it recovered 31 nats
  # on the real fit, i.e. slightly more than the extrapolated remainder, the difference being
  # cross-block curvature the sweeps cannot see by construction.
  # So reporting `converged = sweeps < max_sweeps` named stage one and read FALSE on a fit where all
  # three stages had succeeded, which trains a reader to ignore the flag. It also had an off-by-one:
  # breaking exactly at the last sweep counted as failure. `converged` now means what a reader
  # expects -- every optimiser stage returned success and no country is stuck flat -- and the sweep
  # loop's own status is reported separately, with the size of the tail it left behind.
  conv_shared_code = if (is.null(os)) 99L else os$convergence
  conv_polish_code = if (is.null(op)) 99L else op$convergence
  # geometric extrapolation of the block-descent tail from the last two sweep gains, so the reader can
  # judge what the polish had to absorb rather than take "ran out of sweeps" on trust
  r = if (is.finite(last_gain) && is.finite(prev_gain) && prev_gain > 0) last_gain / prev_gain else NA_real_
  sweep_tail = if (is.finite(r) && r > 0 && r < 1) last_gain * r / (1 - r) else NA_real_
  list(theta = th, negll = obj, loglik = jm_loglik_cpp(th, d), seconds = secs, trace = trace,
       sweeps = sw, conv_local = conv_local,
       conv_shared = conv_shared_code, conv_polish = conv_polish_code, polish_gain = polish_gain,
       flat = flat_final, n_flat_unresolved = sum(flat_final$flat),
       flat_report = rbind(if (nrow(prot1$report)) cbind(prot1$report, pass = "pre-polish"),
                           if (nrow(prot2$report)) cbind(prot2$report, pass = "post-polish")),
       flat_thresholds = c(attack_min = attack_min, epi_frac_min = epi_frac_min),
       # stage one on its own: did the sweep loop reach `tol`, and how big a tail did it leave?
       sweeps_hit_tol = hit_tol, sweep_last_gain = last_gain, sweep_tail_nats = sweep_tail,
       converged = isTRUE(all(conv_local == 0L) && conv_shared_code == 0L &&
                          conv_polish_code == 0L && sum(flat_final$flat) == 0L),
       d = d)
}

# ---- |-tidy the parameters ----
jm_unpack = function(th, d){
  S = d$n_season
  x_free = th[seq_len(S - 1)]; dev_free = th[S - 1L + seq_len(S - 1)]
  x = c(x_free, -sum(x_free))
  sea_S0 = !identical(d$season_on, "R0"); cty_S0 = !identical(d$country_on, "R0")
  out = list(season_on = if (sea_S0) "S0" else "R0", country_on = if (cty_S0) "S0" else "R0",
             R0 = unname(d$R0_fixed), S0_fixed = unname(d$S0_fixed),   # the anchors
             x = setNames(x, d$seasons),           # the season effect, on logit S0 or on log R0
             delta = setNames(c(dev_free, -sum(dev_free)), d$seasons),
             sigma_eld = unname(2^th[2L * S - 1L]))
  cty = lapply(seq_len(d$n_country), function(ic){
    base = d$off_country[ic]; ics = d$cs_of_country[[ic]] + 1L
    ss = d$cs_season[ics] + 1L; sn = d$seasons[ss]
    # unname the scalars: they would otherwise carry the parameter's own label ("DK:logit_S0") and
    # leak it into every data frame, summary column and plot label built from them.
    # Each quantity is its anchor unless it carries an effect: the country's LEVEL (at the average
    # season) and its per-season value, for both S0 and R0
    logit_S0_c = if (cty_S0) th[base + 1] else qlogis(d$S0_fixed)
    log_R0_c   = if (cty_S0) log(d$R0_fixed) else th[base + 1]
    list(S0 = unname(plogis(logit_S0_c)),
         # a zero PER SEASON when the effect lives elsewhere, so the vector keeps one value per season
         S0_season = setNames(plogis(logit_S0_c + if (sea_S0) x[ss] else 0 * ss), sn),
         R0 = unname(exp(log_R0_c)),
         R0_season = setNames(exp(log_R0_c + if (sea_S0) 0 * ss else x[ss]), sn),
         c = unname(exp(th[base + 2])),
         off_young = unname(th[base + 3]), off_eld = unname(th[base + 4]),
         phi = unname(exp(th[base + 5])),
         b = setNames(exp(th[base + 5 + seq_len(d$n_src[ic])]), d$sources[[ic]]),
         I0 = setNames(exp(th[base + 5 + d$n_src[ic] + seq_along(ics)]), d$seasons[d$cs_season[ics] + 1L]))
  })
  names(cty) = d$countries
  out$country = cty
  out
}

jm_summary_country = function(fit){
  p = jm_unpack(fit$theta, fit$d); d = fit$d
  data.frame(country = d$countries,
             S0 = vapply(p$country, `[[`, numeric(1), "S0"),
             c_adult = vapply(p$country, `[[`, numeric(1), "c"),
             rel_young = 2^vapply(p$country, `[[`, numeric(1), "off_young"),
             rel_elderly = 2^vapply(p$country, `[[`, numeric(1), "off_eld"),
             phi = vapply(p$country, `[[`, numeric(1), "phi"),
             n_seasons = d$n_cs_of_country,
             # whose mixing pattern this country's age offsets were fitted under: a country on the EU
             # average has no contact matrix of its own, and the offsets are the parameters most
             # sensitive to it, so the provenance travels with the number it qualifies
             contact = if (is.null(d$contact_source)) NA_character_ else unname(d$contact_source),
             row.names = NULL)
}
# A season-level number must never be read without its sample size. Season support is markedly
# uneven -- the last season rests on 7 of 12 countries and roughly half the weekly cells of the
# richest -- and it carries the highest fitted R0, so n_country and obs_cells belong in the same
# table as R0 rather than in a caveat somewhere else.
jm_summary_season = function(fit){
  d = fit$d; p = jm_unpack(fit$theta, d)
  per = function(s, f) { ii = which(d$cs_season == s - 1L); f(ii) }
  # S0_typical: the season's susceptibility for a typical country, i.e. at the median country level.
  # This is the number to read: x is a logit shift and hard to picture on its own.
  med_logit = median(vapply(p$country, function(q) qlogis(q$S0), numeric(1)))
  med_logR0 = median(vapply(p$country, function(q) log(q$R0), numeric(1)))
  sea_S0 = p$season_on == "S0"
  data.frame(season = d$seasons, x = unname(p$x),
             S0_typical = plogis(med_logit + if (sea_S0) unname(p$x) else 0),
             R0_typical = exp(med_logR0 + if (sea_S0) 0 else unname(p$x)),
             deviation = unname(p$delta), reporting_mult = exp(unname(p$delta)),
             n_country = vapply(seq_len(d$n_season), function(s) per(s, length), integer(1)),
             obs_cells = vapply(seq_len(d$n_season), function(s)
               per(s, function(ii) sum(vapply(ii, function(i) sum(is.finite(d$y[[i]])), numeric(1)))), numeric(1)),
             last_week_min = vapply(seq_len(d$n_season), function(s)
               per(s, function(ii) min(d$n_weeks[ii])), integer(1)),
             last_week_max = vapply(seq_len(d$n_season), function(s)
               per(s, function(ii) max(d$n_weeks[ii])), integer(1)),
             row.names = NULL)
}

# ---- |-base-R reference implementation of the SAME model, for the identity test ----
# Deliberately written straight from the maths rather than by translating the C++, and using R's own
# eigen() for the spectral radius rather than the C++ power iteration, so the two agree only if both
# are right. tests/testthat/test-joint-model.R requires 1e-10.
jm_negll_R = function(th, d){
  S = d$n_season; A = 3L
  x_free = th[seq_len(S - 1)]; xs = c(x_free, -sum(x_free))
  dev_free = th[S - 1L + seq_len(S - 1)]
  dev = exp(c(dev_free, -sum(dev_free)))
  sigma = c(1, 1, 2^th[2L * S - 1L])
  dn = function(x, m, s) -0.5 * ((x - m) / s)^2 - log(s) - 0.5 * log(2 * pi)
  sea_S0 = !identical(d$season_on, "R0"); cty_S0 = !identical(d$country_on, "R0")
  lp = sum(dn(xs, 0, if (sea_S0) d$pr_x_sd else d$pr_r_sd)) +
       sum(dn(c(dev_free, -sum(dev_free)), 0, d$pr_delta_sd)) +
       dn(th[2L * S - 1L], d$pr_sigma_mean, d$pr_sigma_sd)
  for (ic in seq_len(d$n_country)){
    base = d$off_country[ic]; nsrc = d$n_src[ic]
    slot0 = th[base + 1]; cc = exp(th[base + 2])
    logit_S0_c = if (cty_S0) slot0 else qlogis(d$S0_fixed)
    log_R0_c   = if (cty_S0) log(d$R0_fixed) else slot0
    c_age = cc * 2^c(th[base + 3], 0, th[base + 4])
    phi = exp(th[base + 5])
    b = exp(th[base + 5 + seq_len(nsrc)])
    I0v = exp(th[base + 5 + nsrc + seq_len(d$n_cs_of_country[ic])])
    Cs = sweep(d$Cn[[ic]], 1, sigma, "*")
    Cs = Cs / max(abs(eigen(Cs, only.values = TRUE)$values))
    N = d$N[[ic]]
    for (ics in d$cs_of_country[[ic]] + 1L){
      y = d$y[[ics]]; nw = nrow(y); s = d$cs_season[ics] + 1L
      beta = exp(log_R0_c + if (sea_S0) 0 else xs[s]) * d$gamma
      S0 = plogis(logit_S0_c + if (sea_S0) xs[s] else 0)   # this country, this season
      Su = rep(S0, A); Iu = rep(I0v[d$cs_pos[ics] + 1L], A); Sv = rep(0, A); Iv = rep(0, A)
      vax = c(0, 0, d$vax_eld[ics]); inc = matrix(0, nw, A); day = 0L
      for (t in seq_len(nw)){
        acc = rep(0, A)
        for (k in 1:7){
          day = day + 1L
          lam = beta * as.numeric(Cs %*% (Iu + (1 - d$ve_spread) * Iv))
          # cap the flow at what the pool holds, exactly as the .cpp does: without it the
          # accumulator banks infections implying a negative S_u that the clamp then undoes
          nu = pmin(lam * Su, Su); nv = pmin((1 - d$ve_inf) * lam * Sv, Sv)
          Su = Su - nu; Iu = Iu + nu - d$gamma * Iu
          Sv = Sv - nv; Iv = Iv + nv - d$gamma * Iv
          acc = acc + nu + (1 - d$ve_ili) * nv
          if (day == d$vax_day){ mv = vax * Su; Su = Su - mv; Sv = Sv + mv }
          Su = pmin(pmax(Su, 0), 1); Iu = pmin(pmax(Iu, 0), 1)
          Sv = pmin(pmax(Sv, 0), 1); Iv = pmin(pmax(Iv, 0), 1)
        }
        inc[t, ] = acc
      }
      mu = sweep(inc, 2, c_age * dev[s] * N, "*") +
           matrix(b[d$cs_src[ics] + 1L] * N / d$rate_per, nw, A, byrow = TRUE)
      mu[mu < 1e-10] = 1e-10
      ok = is.finite(y)
      lp = lp + sum(lgamma(y[ok] + phi) - lgamma(phi) + phi * (log(phi) - log(phi + mu[ok])) +
                    y[ok] * (log(mu[ok]) - log(phi + mu[ok]))) - d$lgamma_y1[ics]
      lp = lp + dn(log(I0v[d$cs_pos[ics] + 1L]), d$pr_I0_mean, d$pr_I0_sd)
    }
    lp = lp + (if (cty_S0) dn(th[base + 1], d$pr_S0_mean, d$pr_S0_sd)
               else dn(th[base + 1], d$pr_R0c_mean, d$pr_R0c_sd)) +
         dn(th[base + 2], d$pr_c_mean, d$pr_c_sd) +
         dn(th[base + 3], 0, d$pr_off_sd) + dn(th[base + 4], 0, d$pr_off_sd) +
         dn(th[base + 5], d$pr_phi_mean, d$pr_phi_sd) +
         sum(dn(th[base + 5 + seq_len(nsrc)], d$pr_b_mean, d$pr_b_sd))
  }
  unname(-lp)              # a log-posterior is a scalar: strip the name lp inherits from theta
}

# ---- |-convergence and identifiability indicators ----
# The likelihood-only Hessian is the curvature of the fit-quality landscape with the priors switched
# off: steep directions are combinations the data determine, near-flat ones are combinations the data
# cannot see. Comparing it with the penalised Hessian shows what the priors are holding up.
# `contraction` = 1 - sd_posterior/sd_prior: 0 means the posterior is just the prior, 1 means the data
# did the work. `sd_data` is the data-only precision with the other parameters held fixed.
jm_identifiability = function(fit, verbose = TRUE){
  d = fit$d; th = fit$theta; nm = names(th)
  if (verbose) cat("Hessians (", length(th), "parameters, ~", round(2 * length(th)^2 / 1000), "k evaluations )...\n")
  t0 = Sys.time()
  H_post = optimHess(th, function(x) jm_negll_cpp(x, d))
  H_lik  = optimHess(th, function(x) -jm_loglik_cpp(x, d))
  H_post = (H_post + t(H_post)) / 2; H_lik = (H_lik + t(H_lik)) / 2
  secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  S = d$n_season
  fam = ifelse(grepl("^x_", nm), "S0 season effect (shared)",
        ifelse(grepl("^r_", nm), "R0 season effect (shared)",
        ifelse(grepl(":log_R0$", nm), "R0 (country)",
        ifelse(grepl("^delta_", nm), "season deviation (shared)",
        ifelse(grepl("log2_sigma", nm), "elderly susceptibility (global)",
        ifelse(grepl(":logit_S0", nm), "S0 (country)",
        ifelse(grepl(":log_c$", nm), "reporting c (country)",
        ifelse(grepl(":off_", nm), "age reporting offset",
        ifelse(grepl(":log_phi", nm), "dispersion phi",
        ifelse(grepl(":log_b_", nm), "baseline b",
        ifelse(grepl(":log_I0_", nm), "seed I0 (country-season)", nm)))))))))))
  pr_sd = ifelse(grepl("^x_", nm), d$pr_x_sd,
          ifelse(grepl("^r_", nm), d$pr_r_sd,
          ifelse(grepl(":log_R0$", nm), d$pr_R0c_sd,
          ifelse(grepl("^delta_", nm), d$pr_delta_sd,
          ifelse(grepl("log2_sigma", nm), d$pr_sigma_sd,
          ifelse(grepl(":logit_S0", nm), d$pr_S0_sd,
          ifelse(grepl(":log_c$", nm), d$pr_c_sd,
          ifelse(grepl(":off_", nm), d$pr_off_sd,
          ifelse(grepl(":log_phi", nm), d$pr_phi_sd,
          ifelse(grepl(":log_b_", nm), d$pr_b_sd,
          ifelse(grepl(":log_I0_", nm), d$pr_I0_sd, NA_real_)))))))))))
  e_lik = eigen(H_lik, symmetric = TRUE, only.values = TRUE)$values
  e_post = eigen(H_post, symmetric = TRUE, only.values = TRUE)$values
  V_post = tryCatch(solve(H_post), error = function(e) NULL)
  sd_post = if (is.null(V_post)) rep(NA_real_, length(th)) else sqrt(pmax(diag(V_post), 0))
  sd_data = 1 / sqrt(pmax(diag(H_lik), 1e-300))
  # the shared block on its own, conditional on every local parameter: what pooling actually delivers
  sh = seq_len(2L * S)
  V_sh = tryCatch(solve(H_lik[sh, sh, drop = FALSE]), error = function(e) NULL)
  sd_shared_cond = if (is.null(V_sh)) rep(NA_real_, length(sh)) else sqrt(pmax(diag(V_sh), 0))
  tab = data.frame(parameter = nm, family = fam, prior_sd = pr_sd, sd_data = sd_data,
                   sd_post = sd_post, contraction = 1 - sd_post / pr_sd, row.names = NULL)
  famtab = tab %>% group_by(family) %>%
    summarise(n = n(), prior_sd = first(prior_sd), sd_data_med = median(sd_data),
              sd_post_med = median(sd_post), contraction_med = median(contraction),
              contraction_min = min(contraction), .groups = "drop") %>% arrange(desc(contraction_med))
  list(seconds = secs, eig_lik = e_lik, eig_post = e_post,
       flat_ratio = min(abs(e_lik)) / max(abs(e_lik)),
       n_near_flat = sum(abs(e_lik) < 1e-8 * max(abs(e_lik))),
       n_negative = sum(e_lik < -1e-6 * max(abs(e_lik))),
       post_pd = all(e_post > 0), table = tab, family = as.data.frame(famtab),
       sd_shared_cond = setNames(sd_shared_cond, nm[sh]),
       # the penalised Hessian itself: ~5 minutes to compute and the basis of every interval and
       # correlation downstream, so it travels with the result instead of being recomputed
       H_post = H_post)
}

# ---- |-model adequacy: is the fitted dispersion explainable as MEASUREMENT noise? ----
# The negative binomial has Var = mu + mu^2/phi, so at large counts the coefficient of variation is
# 1/sqrt(phi). The data's own scatter around a 3-week moving average is a lower bound on what any
# mean could achieve. If the fitted phi implies much more noise than that, the excess is the
# deterministic mean misfitting, not measurement -- and that is the case for adding the filter back.
jm_adequacy = function(fit){
  d = fit$d; p = jm_unpack(fit$theta, d)
  phid = jm_phi_data(d)                     # ONE definition, shared with the multi-start and the protector
  phi_fit = vapply(p$country, `[[`, numeric(1), "phi")
  data.frame(country = d$countries, phi_fitted = unname(phi_fit),
             cv_fitted = unname(1 / sqrt(phi_fit)), cv_data = 1 / sqrt(phid),
             phi_data = phid, excess = unname((1 / sqrt(phi_fit)) / (1 / sqrt(phid))),
             row.names = NULL)
}

# ---- |-observed vs fitted, as rates per 100 000, tidy ----
jm_tidy_fit = function(fit){
  d = fit$d; f = jm_fitted_cpp(fit$theta, d)
  do.call(rbind, lapply(seq_len(d$n_cs), function(i){
    ic = d$cs_country[i] + 1L; N = d$N[[ic]]; nw = d$n_weeks[i]
    per100k = d$rate_per / N
    obs = sweep(d$y[[i]], 2, per100k, "*"); fitv = sweep(f$mu[[i]], 2, per100k, "*")
    data.frame(country = d$countries[ic], season = d$seasons[d$cs_season[i] + 1L],
               week = rep(seq_len(nw), 2 * length(N)),
               group = rep(rep(d$groups, each = nw), 2),
               what = rep(c("observed", "model"), each = nw * length(N)),
               value = c(as.numeric(obs), as.numeric(fitv)), stringsAsFactors = FALSE)
  }))
}
