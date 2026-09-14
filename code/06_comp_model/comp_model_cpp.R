# comp_model_cpp.R -- load the C++ engine and expose it behind the same interface as the R reference
#
# cm_load_cpp() compiles code/06_comp_model/comp_model_core.cpp once (Rcpp::sourceCpp, cached under
# output/comp_model/cpp_cache) and defines cm_ekf_season_cpp(). cm_ekf_season_engine() dispatches on
# settings$engine ("R" = the reference in comp_model_core.R, "cpp" = the port); the fitter calls
# only this dispatcher. tests/testthat/test-comp-model-cpp.R pins the two engines to 1e-10 relative.

cm_load_cpp = function(cache = "output/comp_model/cpp_cache", quiet = TRUE){
  if (exists("cm_ekf_season_cpp", mode = "function")) return(invisible(TRUE))
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp("code/06_comp_model/comp_model_core.cpp", cacheDir = cache, verbose = !quiet, showOutput = !quiet)
  invisible(TRUE)
}

cm_ekf_season_engine = function(y, f, S0, R0, c, b, phi, q, vax_day = NA, vax_frac = NULL, engine = "R", I0 = f$I0){
  if (engine == "R") return(cm_ekf_season(y, f, S0, R0, c, b, phi, q, vax_day, vax_frac, I0 = I0))
  if (!exists("cm_ekf_season_cpp", mode = "function")) cm_load_cpp()
  cm_ekf_season_cpp(y, f$Cn, f$N, f$gamma, f$ve_inf, f$ve_spread, f$ve_ili, I0, f$p0, f$rate_per,
                    S0, R0, rep_len(as.numeric(c), f$n_age), b, phi, q, if (is.na(vax_day)) -1L else as.integer(vax_day),
                    if (is.null(vax_frac)) rep(0, f$n_age) else as.numeric(vax_frac))
}

# deterministic season through either engine (used by the stage-1 fit and the report)
cm_simulate_season_engine = function(f, S0, R0, n_weeks, vax_day = NA, vax_frac = NULL, I0 = f$I0, engine = "R"){
  if (engine == "R") return(cm_simulate_season(f, S0, R0, n_weeks, vax_day, vax_frac, I0))
  if (!exists("cm_simulate_season_cpp", mode = "function")) cm_load_cpp()
  cm_simulate_season_cpp(as.integer(n_weeks), f$Cn, f$N, f$gamma, f$ve_inf, f$ve_spread, f$ve_ili, I0, S0, R0,
                         if (is.na(vax_day)) -1L else as.integer(vax_day),
                         if (is.null(vax_frac)) rep(0, f$n_age) else as.numeric(vax_frac))
}
