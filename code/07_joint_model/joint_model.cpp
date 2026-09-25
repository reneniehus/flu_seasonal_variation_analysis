// joint_model.cpp -- the ENTIRE log-posterior of the joint EU/EEA influenza model, in one place.
//
// Why everything is here: profiling the compartmental pilot showed the C++ engine was only 58% of an
// objective evaluation, the surrounding R loop was the rest, and 28% of every call was spent building
// per-season trajectories that the optimiser threw away. So this file takes the parameter vector and
// the assembled data and returns ONE number. R does no per-evaluation work at all.
//
// Model (see MODEL.md): an age (3) x vaccination (2) SIR per country-season, daily forward Euler,
// weekly aggregation, negative-binomial observation. No Kalman filter, hence no Jacobian: the pilot's
// hand-derived Jacobian and its 245 lines are gone, which is most of why this file is short.
//
// The reference implementation of the SAME model in base R is jm_negll_R() in joint_model.R, and
// tests/testthat/test-joint-model.R requires the two to agree to 1e-10. Every model change must be
// made in both. That rule caught a real bug in the pilot (a clamp floor that re-seeded the vaccinated
// infectious every week) and is cheap to keep while the C++ has no derivatives to get wrong.
//
// Parameter vector layout (all on unconstrained scales), mirrored by jm_pack/jm_unpack in R:
//   [0 .. S-2]         x_s, s = 1..S-1               season effect on SUSCEPTIBILITY, logit scale,
//                                                    free values; x_S = -sum(free), so they average zero
//   [S-1 .. 2S-3]      delta_s, s = 1..S-1           season effect on VISIBILITY, log scale, free values;
//                                                    delta_S = -sum(free), i.e. they average zero
//   [2S-2]             log2 sigma_eld                one global elderly susceptibility
//   then per country c, starting at off_country[c]:
//   +0                 logit S0_c   the country's susceptibility at the average season;
//                                   logit S0_{c,s} = logit S0_c + x_s
//   +1                 log c_c      the country's reporting level at the average season;
//                                   log c_{c,s} = log c_c + delta_s
//
// WHICH QUANTITY SENSES WHAT is a switch, read from d["season_on"] and d["country_on"], each "S0" or
// "R0". The parameter LAYOUT is identical in every setting -- the same slots, the same count -- only
// the meaning of two slot sets changes:
//   season slots [0..S-2]  x_s, a logit shift on S0   when season_on  == "S0"
//                          r_s, a log   shift on R0   when season_on  == "R0"
//   country slot  +0       logit S0_c                 when country_on == "S0"
//                          log   R0_c                 when country_on == "R0"
// Whatever carries no effect is FIXED: d["R0_fixed"] (1.5) and/or d["S0_fixed"] (0.75). So the four
// combinations are four models with the same parameter count, comparable by likelihood:
//   S0/S0 (the working model), R0/R0, and the two hybrids. Transmissibility and susceptibility enter
// the rise rate as a product and are indistinguishable from one wave; they differ in whether the
// POOL left to infect changes too (S0) or only the speed (R0). That difference is what the
// comparison tests (MODEL.md).
//   +2                 off_young                     log2 reporting offset vs adults
//   +3                 off_eld
//   +4                 log phi_c
//   +5 .. 4+n_src[c]   log b_{c,src}                 one per data source the country actually has
//   then n_cs[c]       log I0_{c,s}                  one seed per season of that country

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static const int A = 3;                       // young, medium, elderly -- fixed, so all A-loops unroll
static const double LOG_2PI = 1.8378770664093453;

// ---- |-dominant eigenvalue of a 3x3 positive matrix by power iteration ----
// The contact matrix is rescaled to spectral radius 1 so that the realised R0 equals R0_s exactly
// (MODEL.md). Scaling its rows by the age susceptibilities changes that radius, and sigma_eld is
// fitted, so the radius has to be recomputed every evaluation. Power iteration on 3x3 costs nothing
// and avoids a dependency on an eigen solver.
static inline double spectral_radius3(const double M[A][A]){
  double v[A] = {1.0, 1.0, 1.0}, lam = 0.0;
  for (int it = 0; it < 400; ++it){
    double w[A];
    for (int a = 0; a < A; ++a){ w[a] = 0.0; for (int j = 0; j < A; ++j) w[a] += M[a][j] * v[j]; }
    double nrm = std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
    if (!(nrm > 0.0)) return 0.0;
    for (int a = 0; a < A; ++a) v[a] = w[a] / nrm;
    if (it > 3 && std::fabs(nrm - lam) < 1e-15 * (nrm > 1.0 ? nrm : 1.0)){ return nrm; }
    lam = nrm;
  }
  return lam;
}

// ---- |-one country-season: weekly observation-relevant incidence as a fraction of each age group ----
// inc is n_weeks x A in COLUMN-MAJOR order, matching an R matrix, so it can be handed straight back.
// attack is the per-age cumulative infected fraction (S only leaves S_u/S_v by infection).
//
// TWO HORIZONS, and why. n_weeks is however many weeks that country's surveillance series happens to
// run, 33 to 53 across the panel. Integrating only that far made the attack rate a quantity measured
// over the OBSERVED WINDOW while figure 11 called it "over the season" and compared it with cohort
// evidence: for 5 of 86 country-seasons the modelled epidemic was still running when the window ended,
// and the reported number moved by up to 15% (IT 2015/2016, window 38 weeks) purely because Italy's
// series stopped early. n_weeks_dyn (>= n_weeks) is the horizon the DYNAMICS run to, so attack is a
// season quantity for every country-season alike. Only the first n_weeks are written to inc, so the
// likelihood is bit-identical -- verified: padding every window to 53 weeks changed the negative
// log-likelihood by exactly 0.
static inline void simulate_season(int n_weeks, int n_weeks_dyn, const double Cs[A][A], double beta,
                                   double S0, double I0, double gamma,
                                   double ve_inf, double ve_ili, double ve_spread,
                                   int vax_day, double vax_eld,
                                   double* inc, double* attack){
  double Su[A], Iu[A], Sv[A], Iv[A], acc[A];
  for (int a = 0; a < A; ++a){ Su[a] = S0; Iu[a] = I0; Sv[a] = 0.0; Iv[a] = 0.0; }
  const double vax[A] = {0.0, 0.0, vax_eld};      // only the 65+ group is vaccinated
  const double s_spread = 1.0 - ve_spread, e_inf = 1.0 - ve_inf, w_ili = 1.0 - ve_ili;
  int day = 0;
  for (int t = 0; t < n_weeks_dyn; ++t){
    for (int a = 0; a < A; ++a) acc[a] = 0.0;    // the accumulator resets every observation week
    for (int k = 0; k < 7; ++k){
      ++day;
      double lam[A];                              // force of infection, from the state BEFORE the step
      for (int a = 0; a < A; ++a){
        double s = 0.0;
        for (int j = 0; j < A; ++j) s += Cs[a][j] * (Iu[j] + s_spread * Iv[j]);
        lam[a] = beta * s;
      }
      for (int a = 0; a < A; ++a){
        // cap the flow at what the pool holds: without this the accumulator banks infections that
        // imply a negative S_u, which the clamp below then silently undoes, so the weekly
        // observation could exceed the supply while attack[] stayed correctly bounded
        double nu = lam[a] * Su[a];               // new infections, unvaccinated
        double nv = e_inf * lam[a] * Sv[a];       // ... and vaccinated, who are partly protected
        if (nu > Su[a]) nu = Su[a];
        if (nv > Sv[a]) nv = Sv[a];
        Su[a] -= nu;  Iu[a] += nu - gamma * Iu[a];
        Sv[a] -= nv;  Iv[a] += nv - gamma * Iv[a];
        acc[a] += nu + w_ili * nv;                // vaccinated infections are less likely to be ILI
      }
      if (day == vax_day){
        for (int a = 0; a < A; ++a){ const double mv = vax[a] * Su[a]; Su[a] -= mv; Sv[a] += mv; }
      }
      for (int a = 0; a < A; ++a){                // clamp to [0,1]; the floor is exactly 0
        if (Su[a] < 0.0) Su[a] = 0.0; else if (Su[a] > 1.0) Su[a] = 1.0;
        if (Iu[a] < 0.0) Iu[a] = 0.0; else if (Iu[a] > 1.0) Iu[a] = 1.0;
        if (Sv[a] < 0.0) Sv[a] = 0.0; else if (Sv[a] > 1.0) Sv[a] = 1.0;
        if (Iv[a] < 0.0) Iv[a] = 0.0; else if (Iv[a] > 1.0) Iv[a] = 1.0;
      }
    }
    // only the observed weeks are recorded; the tail past n_weeks advances the state for attack[] only
    if (t < n_weeks) for (int a = 0; a < A; ++a) inc[t + a * n_weeks] = acc[a];
  }
  for (int a = 0; a < A; ++a) attack[a] = S0 - Su[a] - Sv[a];
}

// F1: every export must check the length before any read. country_lp indexes up to exactly n_par-1
// for the last country, so there is no slack and a short theta reads adjacent memory.
static inline void check_len(const NumericVector& theta, const List& d){
  if (theta.size() != as<int>(d["n_par"]))
    stop("theta has the wrong length: %d supplied, %d expected", theta.size(), as<int>(d["n_par"]));
}

static inline double dnorm_log(double x, double m, double s){
  const double z = (x - m) / s;
  return -0.5 * z * z - std::log(s) - 0.5 * LOG_2PI;
}

// ---- |-the shared block: the two season effects and the global elderly susceptibility ----
struct Shared {
  std::vector<double> xs, dev;                  // xs[s]: logit-scale season effect on S0;
  double sigma_eld;                             // dev[s] = exp(delta_s), the reporting multiplier
};
// finiteness of the TRANSFORMED shared values. theta = 800 is a finite number whose exp() is Inf,
// which then produced Inf*0 = NaN inside the dynamics, so checking theta alone is not enough.
static inline bool shared_ok(const Shared& sh){
  if (!R_finite(sh.sigma_eld) || sh.sigma_eld <= 0.0) return false;
  for (size_t s = 0; s < sh.xs.size(); ++s)
    if (!R_finite(sh.xs[s]) || !R_finite(sh.dev[s]) || sh.dev[s] < 0.0) return false;
  return true;
}

static inline Shared read_shared(const double* th, int S){
  Shared sh; sh.xs.resize(S); sh.dev.resize(S);
  double sum_x = 0.0, sum_d = 0.0;
  for (int s = 0; s < S - 1; ++s){ sh.xs[s] = th[s]; sum_x += th[s]; }
  sh.xs[S - 1] = -sum_x;                        // the sum-to-zero constraint, on the logit scale
  for (int s = 0; s < S - 1; ++s){ const double dd = th[S - 1 + s]; sh.dev[s] = std::exp(dd); sum_d += dd; }
  sh.dev[S - 1] = std::exp(-sum_d);             // the sum-to-zero constraint, on the log scale
  sh.sigma_eld = std::exp2(th[2 * S - 2]);
  return sh;
}

// ---- |-log-posterior contribution of ONE country: its likelihood plus its own priors ----
// Everything the block optimiser needs: with the shared block fixed this is a self-contained
// objective over that country's local parameters, and the countries do not interact.
static double country_lp(const double* th, const List& d, int ic, const Shared& sh,
                         bool want_fit, List* fit_out){
  const double gamma    = as<double>(d["gamma"]);
  const double ve_inf   = as<double>(d["ve_inf"]);
  const double ve_ili   = as<double>(d["ve_ili"]);
  const double ve_spread= as<double>(d["ve_spread"]);
  const double rate_per = as<double>(d["rate_per"]);
  const int    vax_day  = as<int>(d["vax_day"]);
  // the switch fields default to the working model when absent, so a data object built before the
  // switch existed (the cached fit) keeps evaluating as the S0/S0 model it is -- the same pattern as
  // attack_weeks
  const double R0_fixed = as<double>(d["R0_fixed"]);
  const double S0_fixed = d.containsElementNamed("S0_fixed") ? as<double>(d["S0_fixed"]) : 0.75;
  const bool season_on_S0  = !d.containsElementNamed("season_on")  || as<std::string>(d["season_on"])  == "S0";
  const bool country_on_S0 = !d.containsElementNamed("country_on") || as<std::string>(d["country_on"]) == "S0";
  const double logit_S0_fixed = std::log(S0_fixed / (1.0 - S0_fixed)), log_R0_fixed = std::log(R0_fixed);
  // the season horizon the dynamics run to, so the attack rate does not depend on where a country's
  // surveillance series happens to stop. Absent (an older cached d) falls back to 0 = the old
  // behaviour, which keeps a stale object readable rather than erroring on it.
  const int attack_weeks = d.containsElementNamed("attack_weeks") ? as<int>(d["attack_weeks"]) : 0;

  const List Cn_list = d["Cn"];  const List N_list = d["N"];
  const IntegerVector n_src = d["n_src"];  const IntegerVector off_country = d["off_country"];
  const List cs_of_country = d["cs_of_country"];
  const IntegerVector cs_season = d["cs_season"], cs_src = d["cs_src"], cs_pos = d["cs_pos"];
  const NumericVector vax_eld = d["vax_eld"], lgamma_y1 = d["lgamma_y1"];
  const List y_list = d["y"];

  const NumericMatrix Cn = Cn_list[ic];
  const NumericVector N  = N_list[ic];
  const IntegerVector mine = cs_of_country[ic];
  const int base = off_country[ic], nsrc = n_src[ic];

  // local parameters. Slot 0 is the country's level of whichever quantity carries the country
  // effect (logit S0_c or log R0_c); the season effect is added per country-season below.
  const double slot0 = th[base + 0];
  const double c_c   = std::exp(th[base + 1]);
  const double oy    = th[base + 2], oe = th[base + 3];
  const double phi   = std::exp(th[base + 4]);
  const double* logb = th + base + 5;
  const double* logI0= th + base + 5 + nsrc;

  // reporting proportion by age: adults are the reference, so their offset is fixed at zero
  double c_age[A];
  c_age[0] = c_c * std::exp2(oy);  c_age[1] = c_c;  c_age[2] = c_c * std::exp2(oe);

  // the contact matrix with the age susceptibilities, renormalised so R0 keeps its meaning
  double Cs[A][A];
  const double sig[A] = {1.0, 1.0, sh.sigma_eld};
  for (int a = 0; a < A; ++a) for (int j = 0; j < A; ++j) Cs[a][j] = sig[a] * Cn(a, j);
  const double rho = spectral_radius3(Cs);
  if (!(rho > 0.0) || !R_finite(rho)){
    // write an empty result first: the want_fit caller indexes fit_out unconditionally
    if (want_fit) *fit_out = List::create(_["mu"] = List(0), _["attack"] = NumericMatrix(0, A));
    return R_NegInf;
  }
  for (int a = 0; a < A; ++a) for (int j = 0; j < A; ++j) Cs[a][j] /= rho;

  const double lg_phi = std::lgamma(phi), log_phi = std::log(phi);
  double lp = 0.0;
  std::vector<double> inc;
  List mu_out; NumericMatrix attack_out;
  if (want_fit){ mu_out = List(mine.size()); attack_out = NumericMatrix(mine.size(), A); }

  for (int m = 0; m < mine.size(); ++m){
    const int ics = mine[m];
    const NumericMatrix y = y_list[ics];
    const int nw = y.nrow();
    const int s  = cs_season[ics];
    // this country, this season: each quantity is its fixed anchor unless it carries an effect
    const double logitS0 = (country_on_S0 ? slot0 : logit_S0_fixed) + (season_on_S0 ? sh.xs[s] : 0.0);
    const double logR0   = (country_on_S0 ? log_R0_fixed : slot0)   + (season_on_S0 ? 0.0 : sh.xs[s]);
    const double beta = std::exp(logR0) * gamma;
    const double S0 = 1.0 / (1.0 + std::exp(-logitS0));
    const double I0 = std::exp(logI0[cs_pos[ics]]);
    const double b  = std::exp(logb[cs_src[ics]]);
    const double dev = sh.dev[s];

    inc.assign((size_t)nw * A, 0.0);
    double attack[A];
    // the dynamics run to the season horizon so attack[] is comparable across country-seasons, while
    // only the nw observed weeks feed the likelihood (see simulate_season)
    const int nw_dyn = attack_weeks > nw ? attack_weeks : nw;
    simulate_season(nw, nw_dyn, Cs, beta, S0, I0, gamma, ve_inf, ve_ili, ve_spread, vax_day, vax_eld[ics],
                    inc.data(), attack);

    NumericMatrix mu_m;
    if (want_fit){ mu_m = NumericMatrix(nw, A); }
    for (int a = 0; a < A; ++a){
      const double scale = c_age[a] * dev * N[a];        // infections -> expected positive consultations
      const double floor_a = b * N[a] / rate_per;        // the off-season floor of the observed series
      for (int t = 0; t < nw; ++t){
        double mu = scale * inc[t + a * nw] + floor_a;
        if (ISNAN(mu)) return R_NegInf;           // reject a NaN trajectory; do not floor it
        if (mu < 1e-10) mu = 1e-10;               // and floor a genuinely tiny mean
        if (want_fit) mu_m(t, a) = mu;
        const double yy = y(t, a);
        if (ISNAN(yy)) continue;
        // negative binomial with mean mu and dispersion phi: Var = mu + mu^2/phi. The lgamma(y+1)
        // term is constant in the parameters and is added back once, outside, from lgamma_y1.
        lp += std::lgamma(yy + phi) - lg_phi + phi * (log_phi - std::log(phi + mu))
              + yy * (std::log(mu) - std::log(phi + mu));
      }
    }
    lp -= lgamma_y1[ics];
    if (want_fit){ mu_out[m] = mu_m; for (int a = 0; a < A; ++a) attack_out(m, a) = attack[a]; }
    // per-season prior: the seed
    lp += dnorm_log(logI0[cs_pos[ics]], as<double>(d["pr_I0_mean"]), as<double>(d["pr_I0_sd"]));
  }

  // the country's own priors
  if (country_on_S0) lp += dnorm_log(th[base + 0], as<double>(d["pr_S0_mean"]),  as<double>(d["pr_S0_sd"]));
  else               lp += dnorm_log(th[base + 0], as<double>(d["pr_R0c_mean"]), as<double>(d["pr_R0c_sd"]));
  lp += dnorm_log(th[base + 1], as<double>(d["pr_c_mean"]),  as<double>(d["pr_c_sd"]));
  lp += dnorm_log(oy, 0.0, as<double>(d["pr_off_sd"]));
  lp += dnorm_log(oe, 0.0, as<double>(d["pr_off_sd"]));
  lp += dnorm_log(th[base + 4], as<double>(d["pr_phi_mean"]), as<double>(d["pr_phi_sd"]));
  for (int k = 0; k < nsrc; ++k) lp += dnorm_log(logb[k], as<double>(d["pr_b_mean"]), as<double>(d["pr_b_sd"]));

  if (want_fit){ *fit_out = List::create(_["mu"] = mu_out, _["attack"] = attack_out); }
  return lp;
}

// ---- |-priors on the shared block ----
// Both season-effect sets are centred on zero, and the CONSTRAINED last member of each is penalised
// too, so all S values of each set sit under the same prior.
static double shared_lp(const double* th, const List& d, int S){
  double lp = 0.0;
  // the season effect's prior sd depends on which scale it lives on
  const bool season_on_S0 = !d.containsElementNamed("season_on") || as<std::string>(d["season_on"]) == "S0";
  const double s_x = season_on_S0 ? as<double>(d["pr_x_sd"]) : as<double>(d["pr_r_sd"]);
  const double s_dev = as<double>(d["pr_delta_sd"]);
  double sum_x = 0.0, sum_d = 0.0;
  for (int s = 0; s < S - 1; ++s){ lp += dnorm_log(th[s], 0.0, s_x); sum_x += th[s]; }
  lp += dnorm_log(-sum_x, 0.0, s_x);
  for (int s = 0; s < S - 1; ++s){ lp += dnorm_log(th[S - 1 + s], 0.0, s_dev); sum_d += th[S - 1 + s]; }
  lp += dnorm_log(-sum_d, 0.0, s_dev);
  lp += dnorm_log(th[2 * S - 2], as<double>(d["pr_sigma_mean"]), as<double>(d["pr_sigma_sd"]));
  return lp;
}

// ============================ exported ============================

//' Negative log-posterior of the whole joint model. One call, one number.
// [[Rcpp::export]]
double jm_negll_cpp(NumericVector theta, List d){
  const int S = as<int>(d["n_season"]), C = as<int>(d["n_country"]);
  check_len(theta, d);
  const double* th = theta.begin();
  for (int i = 0; i < theta.size(); ++i) if (!R_finite(th[i])) return 1e10;
  const Shared sh = read_shared(th, S);
  if (!shared_ok(sh)) return 1e10;
  double lp = shared_lp(th, d, S);
  for (int ic = 0; ic < C; ++ic){
    const double v = country_lp(th, d, ic, sh, false, nullptr);
    if (!R_finite(v)) return 1e10;
    lp += v;
  }
  if (!R_finite(lp)) return 1e10;
  return -lp;
}

//' Negative log-posterior contribution of ONE country (0-based), its own priors included and the
//' shared-block priors excluded. This is the objective of the per-country step of the block optimiser.
// [[Rcpp::export]]
double jm_country_negll_cpp(NumericVector theta, List d, int ic){
  const int S = as<int>(d["n_season"]);
  check_len(theta, d);
  if (ic < 0 || ic >= as<int>(d["n_country"])) stop("country index out of range");
  const double* th = theta.begin();
  for (int i = 0; i < theta.size(); ++i) if (!R_finite(th[i])) return 1e10;
  const Shared sh = read_shared(th, S);
  if (!shared_ok(sh)) return 1e10;
  const double v = country_lp(th, d, ic, sh, false, nullptr);
  if (!R_finite(v)) return 1e10;
  return -v;
}

//' Pure log-likelihood, priors excluded. For model comparison across variants, where the prior
//' constants differ with the parameter count and would otherwise contaminate the comparison.
// [[Rcpp::export]]
double jm_loglik_cpp(NumericVector theta, List d){
  const int S = as<int>(d["n_season"]), C = as<int>(d["n_country"]);
  check_len(theta, d);
  const double* th = theta.begin();
  const Shared sh = read_shared(th, S);
  if (!shared_ok(sh)) return R_NegInf;
  const double nl = jm_negll_cpp(theta, d);
  // PROPAGATE the rejection instead of stripping priors off the sentinel: subtracting a log-prior
  // from 1e10 yields a finite number that looks exactly like a log-likelihood
  if (nl >= 1e10) return R_NegInf;
  double ll = -nl;                                     // posterior
  ll -= shared_lp(th, d, S);                           // strip the shared priors
  // strip every local prior by recomputing them
  const IntegerVector n_src = d["n_src"], off_country = d["off_country"];
  const List cs_of_country = d["cs_of_country"];
  const IntegerVector cs_pos = d["cs_pos"];
  for (int ic = 0; ic < C; ++ic){
    const int base = off_country[ic], nsrc = n_src[ic];
    const IntegerVector mine = cs_of_country[ic];
    if (!d.containsElementNamed("country_on") || as<std::string>(d["country_on"]) == "S0")
         ll -= dnorm_log(th[base + 0], as<double>(d["pr_S0_mean"]),  as<double>(d["pr_S0_sd"]));
    else ll -= dnorm_log(th[base + 0], as<double>(d["pr_R0c_mean"]), as<double>(d["pr_R0c_sd"]));
    ll -= dnorm_log(th[base + 1], as<double>(d["pr_c_mean"]),  as<double>(d["pr_c_sd"]));
    ll -= dnorm_log(th[base + 2], 0.0, as<double>(d["pr_off_sd"]));
    ll -= dnorm_log(th[base + 3], 0.0, as<double>(d["pr_off_sd"]));
    ll -= dnorm_log(th[base + 4], as<double>(d["pr_phi_mean"]), as<double>(d["pr_phi_sd"]));
    for (int k = 0; k < nsrc; ++k) ll -= dnorm_log(th[base + 5 + k], as<double>(d["pr_b_mean"]), as<double>(d["pr_b_sd"]));
    for (int m = 0; m < mine.size(); ++m)
      ll -= dnorm_log(th[base + 5 + nsrc + cs_pos[mine[m]]], as<double>(d["pr_I0_mean"]), as<double>(d["pr_I0_sd"]));
  }
  return ll;
}

//' Fitted quantities at a parameter vector: expected counts per country-season and per-age attack
//' rates. Returned only when asked for, never during optimisation.
// [[Rcpp::export]]
List jm_fitted_cpp(NumericVector theta, List d){
  const int S = as<int>(d["n_season"]), C = as<int>(d["n_country"]);
  check_len(theta, d);
  const double* th = theta.begin();
  const Shared sh = read_shared(th, S);
  if (!shared_ok(sh))
    stop("cannot compute fitted values: the shared block is not usable (a season effect or the "
         "elderly susceptibility is non-finite or negative)");
  const List cs_of_country = d["cs_of_country"];
  const int n_cs = as<int>(d["n_cs"]);
  List mu_all(n_cs); NumericMatrix attack_all(n_cs, A);
  for (int ic = 0; ic < C; ++ic){
    List fit;
    country_lp(th, d, ic, sh, true, &fit);
    const List mu = fit["mu"]; const NumericMatrix at = fit["attack"];
    const IntegerVector mine = cs_of_country[ic];
    for (int m = 0; m < mine.size(); ++m){
      mu_all[mine[m]] = mu[m];
      for (int a = 0; a < A; ++a) attack_all(mine[m], a) = at(m, a);
    }
  }
  NumericVector xs(S), dev(S);
  for (int s = 0; s < S; ++s){ xs[s] = sh.xs[s]; dev[s] = sh.dev[s]; }
  return List::create(_["mu"] = mu_all, _["attack"] = attack_all,
                      _["x"] = xs, _["dev"] = dev, _["sigma_eld"] = sh.sigma_eld);
}
