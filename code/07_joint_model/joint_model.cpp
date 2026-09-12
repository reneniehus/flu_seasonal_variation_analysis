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
//   [0 .. S-1]         log R0_s                      one per season, shared across countries
//   [S .. 2S-2]        delta_s, s = 1..S-1           season observation deviation, free values
//                                                    delta_S = -sum(free), i.e. they average zero
//   [2S-1]             log2 sigma_eld                one global elderly susceptibility
//   then per country c, starting at off_country[c]:
//   +0                 logit S0_c
//   +1                 log c_c
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
static inline void simulate_season(int n_weeks, const double Cs[A][A], double beta,
                                   double S0, double I0, double gamma,
                                   double ve_inf, double ve_ili, double ve_spread,
                                   int vax_day, double vax_eld,
                                   double* inc, double* attack){
  double Su[A], Iu[A], Sv[A], Iv[A], acc[A];
  for (int a = 0; a < A; ++a){ Su[a] = S0; Iu[a] = I0; Sv[a] = 0.0; Iv[a] = 0.0; }
  const double vax[A] = {0.0, 0.0, vax_eld};      // only the 65+ group is vaccinated
  const double s_spread = 1.0 - ve_spread, e_inf = 1.0 - ve_inf, w_ili = 1.0 - ve_ili;
  int day = 0;
  for (int t = 0; t < n_weeks; ++t){
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
        const double nu = lam[a] * Su[a];         // new infections, unvaccinated
        const double nv = e_inf * lam[a] * Sv[a]; // ... and vaccinated, who are partly protected
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
    for (int a = 0; a < A; ++a) inc[t + a * n_weeks] = acc[a];
  }
  for (int a = 0; a < A; ++a) attack[a] = S0 - Su[a] - Sv[a];
}

static inline double dnorm_log(double x, double m, double s){
  const double z = (x - m) / s;
  return -0.5 * z * z - std::log(s) - 0.5 * LOG_2PI;
}

// ---- |-the shared block: R0_s, the season deviations, the global elderly susceptibility ----
struct Shared {
  std::vector<double> R0, dev;                  // dev[s] = exp(delta_s), the reporting multiplier
  double sigma_eld;
};
static inline Shared read_shared(const double* th, int S){
  Shared sh; sh.R0.resize(S); sh.dev.resize(S);
  double sum_free = 0.0;
  for (int s = 0; s < S; ++s) sh.R0[s] = std::exp(th[s]);
  for (int s = 0; s < S - 1; ++s){ const double d = th[S + s]; sh.dev[s] = std::exp(d); sum_free += d; }
  sh.dev[S - 1] = std::exp(-sum_free);          // the sum-to-zero constraint, in the log scale
  sh.sigma_eld = std::exp2(th[2 * S - 1]);
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

  // local parameters
  const double S0    = 1.0 / (1.0 + std::exp(-th[base + 0]));
  const double c_c   = std::exp(th[base + 1]);
  const double oy    = th[base + 2], oe = th[base + 3];
  const double phi   = std::exp(th[base + 4]);
  const double* logb = th + base + 5;
  const double* logI0= th + base + 5 + nsrc;

  // reporting proportion by age: adults are the reference, so their offset is fixed at zero
  double c_age[A];
  c_age[0] = c_c * std::exp2(oy);  c_age[1] = c_c;  c_age[2] = c_c * std::exp2(oe);

  // the contact matrix with the age susceptibilities, renormalised so R0_s keeps its meaning
  double Cs[A][A];
  const double sig[A] = {1.0, 1.0, sh.sigma_eld};
  for (int a = 0; a < A; ++a) for (int j = 0; j < A; ++j) Cs[a][j] = sig[a] * Cn(a, j);
  const double rho = spectral_radius3(Cs);
  if (!(rho > 0.0) || !R_finite(rho)) return R_NegInf;
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
    const double beta = sh.R0[s] * gamma;
    const double I0 = std::exp(logI0[cs_pos[ics]]);
    const double b  = std::exp(logb[cs_src[ics]]);
    const double dev = sh.dev[s];

    inc.assign((size_t)nw * A, 0.0);
    double attack[A];
    simulate_season(nw, Cs, beta, S0, I0, gamma, ve_inf, ve_ili, ve_spread, vax_day, vax_eld[ics],
                    inc.data(), attack);

    NumericMatrix mu_m;
    if (want_fit){ mu_m = NumericMatrix(nw, A); }
    for (int a = 0; a < A; ++a){
      const double scale = c_age[a] * dev * N[a];        // infections -> expected positive consultations
      const double floor_a = b * N[a] / rate_per;        // the off-season floor of the observed series
      for (int t = 0; t < nw; ++t){
        double mu = scale * inc[t + a * nw] + floor_a;
        if (!(mu > 1e-10)) mu = 1e-10;
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
  lp += dnorm_log(th[base + 0], as<double>(d["pr_S0_mean"]), as<double>(d["pr_S0_sd"]));
  lp += dnorm_log(th[base + 1], as<double>(d["pr_c_mean"]),  as<double>(d["pr_c_sd"]));
  lp += dnorm_log(oy, 0.0, as<double>(d["pr_off_sd"]));
  lp += dnorm_log(oe, 0.0, as<double>(d["pr_off_sd"]));
  lp += dnorm_log(th[base + 4], as<double>(d["pr_phi_mean"]), as<double>(d["pr_phi_sd"]));
  for (int k = 0; k < nsrc; ++k) lp += dnorm_log(logb[k], as<double>(d["pr_b_mean"]), as<double>(d["pr_b_sd"]));

  if (want_fit){ *fit_out = List::create(_["mu"] = mu_out, _["attack"] = attack_out); }
  return lp;
}

// ---- |-priors on the shared block ----
static double shared_lp(const double* th, const List& d, int S){
  double lp = 0.0;
  const double m_R0 = as<double>(d["pr_R0_mean"]), s_R0 = as<double>(d["pr_R0_sd"]);
  const double s_dev = as<double>(d["pr_delta_sd"]);
  for (int s = 0; s < S; ++s) lp += dnorm_log(th[s], m_R0, s_R0);
  double sum_free = 0.0;
  for (int s = 0; s < S - 1; ++s){ lp += dnorm_log(th[S + s], 0.0, s_dev); sum_free += th[S + s]; }
  lp += dnorm_log(-sum_free, 0.0, s_dev);                 // the constrained last deviation
  lp += dnorm_log(th[2 * S - 1], as<double>(d["pr_sigma_mean"]), as<double>(d["pr_sigma_sd"]));
  return lp;
}

// ============================ exported ============================

//' Negative log-posterior of the whole joint model. One call, one number.
// [[Rcpp::export]]
double jm_negll_cpp(NumericVector theta, List d){
  const int S = as<int>(d["n_season"]), C = as<int>(d["n_country"]);
  const double* th = theta.begin();
  if (theta.size() != as<int>(d["n_par"])) stop("theta has the wrong length");
  for (int i = 0; i < theta.size(); ++i) if (!R_finite(th[i])) return 1e10;
  const Shared sh = read_shared(th, S);
  if (!R_finite(sh.sigma_eld) || sh.sigma_eld <= 0.0) return 1e10;
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
  const double* th = theta.begin();
  for (int i = 0; i < theta.size(); ++i) if (!R_finite(th[i])) return 1e10;
  const Shared sh = read_shared(th, S);
  if (!R_finite(sh.sigma_eld) || sh.sigma_eld <= 0.0) return 1e10;
  const double v = country_lp(th, d, ic, sh, false, nullptr);
  if (!R_finite(v)) return 1e10;
  return -v;
}

//' Pure log-likelihood, priors excluded. For model comparison across variants, where the prior
//' constants differ with the parameter count and would otherwise contaminate the comparison.
// [[Rcpp::export]]
double jm_loglik_cpp(NumericVector theta, List d){
  const int S = as<int>(d["n_season"]), C = as<int>(d["n_country"]);
  const double* th = theta.begin();
  const Shared sh = read_shared(th, S);
  double ll = -jm_negll_cpp(theta, d);                 // posterior
  ll -= shared_lp(th, d, S);                           // strip the shared priors
  // strip every local prior by recomputing them
  const IntegerVector n_src = d["n_src"], off_country = d["off_country"];
  const List cs_of_country = d["cs_of_country"];
  const IntegerVector cs_pos = d["cs_pos"];
  for (int ic = 0; ic < C; ++ic){
    const int base = off_country[ic], nsrc = n_src[ic];
    const IntegerVector mine = cs_of_country[ic];
    ll -= dnorm_log(th[base + 0], as<double>(d["pr_S0_mean"]), as<double>(d["pr_S0_sd"]));
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
  const double* th = theta.begin();
  const Shared sh = read_shared(th, S);
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
  NumericVector R0(S), dev(S);
  for (int s = 0; s < S; ++s){ R0[s] = sh.R0[s]; dev[s] = sh.dev[s]; }
  return List::create(_["mu"] = mu_all, _["attack"] = attack_all,
                      _["R0"] = R0, _["dev"] = dev, _["sigma_eld"] = sh.sigma_eld);
}
