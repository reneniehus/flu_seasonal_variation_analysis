// comp_model_core.cpp -- C++ (Rcpp) port of the compartmental model's Kalman loop
//
// THE R VERSION IS THE REFERENCE. code/06_comp_model/comp_model_core.R defines the model; this file
// re-implements cm_day / cm_week / cm_ekf_season for speed (the joint all-countries fit needs it).
// Every change to the model must be made in BOTH files, and tests/testthat/test-comp-model-cpp.R
// requires the two likelihoods and filtered trajectories to agree to 1e-10 relative on real data.
// Operation order mirrors the R code line by line (same clamps, same symmetrisation, same
// Cholesky-based innovation likelihood) so that the agreement is numerical, not approximate.
//
// State layout (n = 5A): S_u[A], I_u[A], S_v[A], I_v[A], C[A]  (fractions of each age group).

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

struct Fixed {
  int A, n;
  std::vector<double> Cn;   // row-major A x A
  std::vector<double> N;
  double gamma, ve_inf, ve_spread, ve_ili, I0, p0, rate_per;
};

static inline double clamp01(double v){ return v < 0 ? 0 : (v > 1 ? 1 : v); }

// one day of forward Euler; J (n x n, row-major) receives the analytical Jacobian if jac
static void cm_day(std::vector<double>& x, const Fixed& f, double beta,
                   const double* vf, bool jac, std::vector<double>& J){
  const int A = f.A, n = f.n;
  const int oSu = 0, oIu = A, oSv = 2*A, oIv = 3*A, oC = 4*A;
  const double s = 1.0 - f.ve_spread, e = 1.0 - f.ve_inf, w = 1.0 - f.ve_ili, g = f.gamma;
  std::vector<double> lam(A), new_u(A), new_v(A), Su(A), Sv(A);
  for (int a = 0; a < A; ++a){
    double acc = 0.0;
    for (int j = 0; j < A; ++j) acc += f.Cn[a*A + j] * (x[oIu + j] + s * x[oIv + j]);
    lam[a] = beta * acc;
    Su[a] = x[oSu + a]; Sv[a] = x[oSv + a];
    new_u[a] = lam[a] * Su[a];
    new_v[a] = e * lam[a] * Sv[a];
  }
  std::vector<double> xn(x);
  for (int a = 0; a < A; ++a){
    xn[oSu + a] = Su[a] - new_u[a];
    xn[oIu + a] = x[oIu + a] + new_u[a] - g * x[oIu + a];
    xn[oSv + a] = Sv[a] - new_v[a];
    xn[oIv + a] = x[oIv + a] + new_v[a] - g * x[oIv + a];
    xn[oC  + a] = x[oC + a] + new_u[a] + w * new_v[a];
  }
  if (vf){
    for (int a = 0; a < A; ++a){
      double moved = vf[a] * xn[oSu + a];
      xn[oSu + a] -= moved; xn[oSv + a] += moved;
    }
  }
  for (int k = 0; k < 4*A; ++k) xn[k] = clamp01(xn[k]);
  if (jac){
    // J = identity, then the blocks of the (pre-clamp, pre-pulse) update
    std::fill(J.begin(), J.end(), 0.0);
    for (int k = 0; k < n; ++k) J[k*n + k] = 1.0;
    for (int a = 0; a < A; ++a){
      // diagonal (lam) terms
      J[(oSu+a)*n + (oSu+a)] = 1.0 - lam[a];
      J[(oIu+a)*n + (oSu+a)] = lam[a];
      J[(oIu+a)*n + (oIu+a)] = 1.0 - g;             // + DSu * dlam_dIu added below
      J[(oSv+a)*n + (oSv+a)] = 1.0 - e * lam[a];
      J[(oIv+a)*n + (oSv+a)] = e * lam[a];
      J[(oIv+a)*n + (oIv+a)] = 1.0 - g;
      J[(oC +a)*n + (oSu+a)] = lam[a];
      J[(oC +a)*n + (oSv+a)] = w * e * lam[a];
      for (int j = 0; j < A; ++j){
        double dIu = beta * f.Cn[a*A + j];            // d lam_a / d I_u[j]
        double dIv = beta * s * f.Cn[a*A + j];        // d lam_a / d I_v[j]
        J[(oSu+a)*n + (oIu+j)] += -Su[a] * dIu;      J[(oSu+a)*n + (oIv+j)] += -Su[a] * dIv;
        J[(oIu+a)*n + (oIu+j)] +=  Su[a] * dIu;      J[(oIu+a)*n + (oIv+j)] +=  Su[a] * dIv;
        J[(oSv+a)*n + (oIu+j)] += -e * Sv[a] * dIu;  J[(oSv+a)*n + (oIv+j)] += -e * Sv[a] * dIv;
        J[(oIv+a)*n + (oIu+j)] +=  e * Sv[a] * dIu;  J[(oIv+a)*n + (oIv+j)] +=  e * Sv[a] * dIv;
        J[(oC +a)*n + (oIu+j)] +=  Su[a] * dIu + w * e * Sv[a] * dIu;
        J[(oC +a)*n + (oIv+j)] +=  Su[a] * dIv + w * e * Sv[a] * dIv;
      }
    }
    if (vf){   // J = Jv %*% J with Jv: S_u' = (1-v) S_u ; S_v' = S_v + v S_u  (row ops on J)
      for (int a = 0; a < A; ++a){
        double v = vf[a];
        for (int c = 0; c < n; ++c){
          double rSu = J[(oSu+a)*n + c];
          J[(oSv+a)*n + c] += v * rSu;
          J[(oSu+a)*n + c]  = (1.0 - v) * rSu;
        }
      }
    }
  }
  x.swap(xn);
}

// one week: reset C, 7 days, pulse on its day; Jw = product of daily Jacobians (with C rows zeroed at reset)
static void cm_week(std::vector<double>& x, const Fixed& f, double beta, int day0, int vax_day,
                    const double* vf, bool jac, std::vector<double>& Jw){
  const int n = f.n, A = f.A;
  for (int a = 0; a < A; ++a) x[4*A + a] = 0.0;
  std::vector<double> Jd, tmp;
  if (jac){
    std::fill(Jw.begin(), Jw.end(), 0.0);
    for (int k = 0; k < n; ++k) Jw[k*n + k] = 1.0;
    for (int a = 0; a < A; ++a) for (int c = 0; c < n; ++c) Jw[(4*A + a)*n + c] = 0.0;
    Jd.assign(n*n, 0.0); tmp.assign(n*n, 0.0);
  }
  for (int d = 1; d <= 7; ++d){
    int day = day0 + d;
    const double* v = (vax_day > 0 && day == vax_day) ? vf : nullptr;
    cm_day(x, f, beta, v, jac, Jd);
    if (jac){   // Jw = Jd %*% Jw
      for (int r = 0; r < n; ++r) for (int c = 0; c < n; ++c){
        double acc = 0.0;
        for (int k = 0; k < n; ++k) acc += Jd[r*n + k] * Jw[k*n + c];
        tmp[r*n + c] = acc;
      }
      Jw.swap(tmp);
    }
  }
}

// [[Rcpp::export]]
List cm_ekf_season_cpp(NumericMatrix y, NumericMatrix Cn, NumericVector N,
                       double gamma, double ve_inf, double ve_spread, double ve_ili, double I0, double p0, double rate_per,
                       double S0, double R0, NumericVector c, double b, double phi, double q,
                       int vax_day, NumericVector vax_frac){
  Fixed f; f.A = Cn.nrow(); f.n = 5 * f.A;
  f.Cn.assign(f.A * f.A, 0.0);
  for (int a = 0; a < f.A; ++a) for (int j = 0; j < f.A; ++j) f.Cn[a*f.A + j] = Cn(a, j);
  f.N.assign(N.begin(), N.end());
  f.gamma = gamma; f.ve_inf = ve_inf; f.ve_spread = ve_spread; f.ve_ili = ve_ili; f.I0 = I0; f.p0 = p0; f.rate_per = rate_per;
  const int A = f.A, n = f.n, T = y.nrow();
  const double beta = R0 * gamma;
  std::vector<double> vf(vax_frac.begin(), vax_frac.end());

  std::vector<double> x(n, 0.0), P(n*n, 0.0), Jw(n*n, 0.0), Ppr(n*n, 0.0), tmp(n*n, 0.0);
  for (int a = 0; a < A; ++a){ x[a] = S0; x[A + a] = I0; }
  for (int a = 0; a < A; ++a){ P[a*n + a] = (p0*S0)*(p0*S0); P[(A+a)*n + (A+a)] = (p0*I0)*(p0*I0); }
  std::vector<double> cN(A), b_cnt(A);
  for (int a = 0; a < A; ++a){ cN[a] = c[a] * f.N[a]; b_cnt[a] = b * f.N[a] / rate_per; }

  NumericMatrix mu_pred(T, A), I_filt(T, A), S_filt(T, A);
  std::fill(mu_pred.begin(), mu_pred.end(), NA_REAL); std::fill(I_filt.begin(), I_filt.end(), NA_REAL); std::fill(S_filt.begin(), S_filt.end(), NA_REAL);
  double ll = 0.0; const double LOG2PI = std::log(2.0 * M_PI);

  for (int t = 0; t < T; ++t){
    cm_week(x, f, beta, 7*t, vax_day, vf.data(), true, Jw);          // x is now xpr
    // Ppr = Jw P Jw^T
    for (int r = 0; r < n; ++r) for (int k = 0; k < n; ++k){ double acc = 0.0; for (int m = 0; m < n; ++m) acc += Jw[r*n + m] * P[m*n + k]; tmp[r*n + k] = acc; }
    for (int r = 0; r < n; ++r) for (int k = 0; k < n; ++k){ double acc = 0.0; for (int m = 0; m < n; ++m) acc += tmp[r*n + m] * Jw[k*n + m]; Ppr[r*n + k] = acc; }
    for (int a = 0; a < A; ++a){                                       // multiplicative process noise on I
      double iu = x[A + a], iv = x[3*A + a];
      Ppr[(A+a)*n + (A+a)]     += (q*iu)*(q*iu);
      Ppr[(3*A+a)*n + (3*A+a)] += (q*iv)*(q*iv);
    }
    std::vector<double> mu(A);
    for (int a = 0; a < A; ++a){ mu[a] = cN[a] * x[4*A + a] + b_cnt[a]; mu_pred(t, a) = mu[a]; }
    std::vector<int> ok;
    for (int a = 0; a < A; ++a) if (R_finite(y(t, a))) ok.push_back(a);
    const int m = ok.size();
    if (m > 0){
      // S = H Ppr H' + diag(Rt), with H[i, C_a] = cN[a]
      std::vector<double> S(m*m, 0.0), innov(m), PHt(n*m, 0.0);
      for (int i = 0; i < m; ++i){
        int ci = 4*A + ok[i];
        for (int r = 0; r < n; ++r) PHt[r*m + i] = Ppr[r*n + ci] * cN[ok[i]];
      }
      for (int i = 0; i < m; ++i) for (int j = 0; j < m; ++j){
        int cj = 4*A + ok[j];
        S[i*m + j] = PHt[cj*m + i] * cN[ok[j]];
      }
      for (int i = 0; i < m; ++i) S[i*m + i] += mu[ok[i]] + mu[ok[i]]*mu[ok[i]] / phi;
      for (int i = 0; i < m; ++i) for (int j = i+1; j < m; ++j){ double v = 0.5*(S[i*m+j] + S[j*m+i]); S[i*m+j] = v; S[j*m+i] = v; }
      // Cholesky S = L L^T
      std::vector<double> Lc(m*m, 0.0); bool okchol = true;
      for (int i = 0; i < m && okchol; ++i){
        for (int j = 0; j <= i; ++j){
          double acc = S[i*m + j];
          for (int k = 0; k < j; ++k) acc -= Lc[i*m + k] * Lc[j*m + k];
          if (i == j){ if (acc <= 0){ okchol = false; break; } Lc[i*m + i] = std::sqrt(acc); }
          else Lc[i*m + j] = acc / Lc[j*m + j];
        }
      }
      if (!okchol) return List::create(_["loglik"] = -1e10, _["mu_pred"] = mu_pred, _["I"] = I_filt, _["S"] = S_filt);
      for (int i = 0; i < m; ++i) innov[i] = y(t, ok[i]) - mu[ok[i]];
      // z = S^{-1} innov via forward/back substitution
      std::vector<double> fw(m), z(m);
      for (int i = 0; i < m; ++i){ double acc = innov[i]; for (int k = 0; k < i; ++k) acc -= Lc[i*m + k] * fw[k]; fw[i] = acc / Lc[i*m + i]; }
      for (int i = m-1; i >= 0; --i){ double acc = fw[i]; for (int k = i+1; k < m; ++k) acc -= Lc[k*m + i] * z[k]; z[i] = acc / Lc[i*m + i]; }
      // Sinv (m x m) for K = PHt Sinv
      std::vector<double> Sinv(m*m, 0.0);
      for (int col = 0; col < m; ++col){
        std::vector<double> e1(m, 0.0), f1(m), g1(m); e1[col] = 1.0;
        for (int i = 0; i < m; ++i){ double acc = e1[i]; for (int k = 0; k < i; ++k) acc -= Lc[i*m + k] * f1[k]; f1[i] = acc / Lc[i*m + i]; }
        for (int i = m-1; i >= 0; --i){ double acc = f1[i]; for (int k = i+1; k < m; ++k) acc -= Lc[k*m + i] * g1[k]; g1[i] = acc / Lc[i*m + i]; }
        for (int i = 0; i < m; ++i) Sinv[i*m + col] = g1[i];
      }
      std::vector<double> K(n*m, 0.0);
      for (int r = 0; r < n; ++r) for (int j = 0; j < m; ++j){ double acc = 0.0; for (int i = 0; i < m; ++i) acc += PHt[r*m + i] * Sinv[i*m + j]; K[r*m + j] = acc; }
      // x = xpr + K innov
      for (int r = 0; r < n; ++r){ double acc = 0.0; for (int i = 0; i < m; ++i) acc += K[r*m + i] * innov[i]; x[r] += acc; }
      // P = (I - K H) Ppr ;  (K H)[r, k] = sum_i K[r,i] * cN[ok_i] * [k == C_{ok_i}]
      for (int r = 0; r < n; ++r) for (int k = 0; k < n; ++k){
        double kh = 0.0;
        for (int i = 0; i < m; ++i) if (k == 4*A + ok[i]) kh += K[r*m + i] * cN[ok[i]];
        tmp[r*n + k] = (r == k ? 1.0 : 0.0) - kh;
      }
      for (int r = 0; r < n; ++r) for (int k = 0; k < n; ++k){ double acc = 0.0; for (int q2 = 0; q2 < n; ++q2) acc += tmp[r*n + q2] * Ppr[q2*n + k]; P[r*n + k] = acc; }
      double logdet = 0.0, quad = 0.0;
      for (int i = 0; i < m; ++i){ logdet += 2.0 * std::log(Lc[i*m + i]); quad += innov[i] * z[i]; }
      ll -= 0.5 * (logdet + quad + m * LOG2PI);
    } else {
      P.swap(Ppr);
    }
    for (int k = 0; k < 4*A; ++k){ if (x[k] < 1e-12) x[k] = 1e-12; if (x[k] > 1) x[k] = 1; }
    for (int a = 0; a < A; ++a){ I_filt(t, a) = x[A + a] + x[3*A + a]; S_filt(t, a) = x[a] + x[2*A + a]; }
  }
  return List::create(_["loglik"] = ll, _["mu_pred"] = mu_pred, _["I"] = I_filt, _["S"] = S_filt);
}

// [[Rcpp::export]]
List cm_simulate_season_cpp(int n_weeks, NumericMatrix Cn, NumericVector N,
                            double gamma, double ve_inf, double ve_spread, double ve_ili, double I0,
                            double S0, double R0, int vax_day, NumericVector vax_frac){
  // deterministic season (no filter): the weekly observation-relevant incidence C per age group,
  // the end state and the per-age attack rate -- mirrors cm_simulate_season() in the R reference
  Fixed f; f.A = Cn.nrow(); f.n = 5 * f.A;
  f.Cn.assign(f.A * f.A, 0.0);
  for (int a = 0; a < f.A; ++a) for (int j = 0; j < f.A; ++j) f.Cn[a*f.A + j] = Cn(a, j);
  f.N.assign(N.begin(), N.end());
  f.gamma = gamma; f.ve_inf = ve_inf; f.ve_spread = ve_spread; f.ve_ili = ve_ili; f.I0 = I0; f.p0 = 0.0; f.rate_per = 1.0;
  const int A = f.A, n = f.n; const double beta = R0 * gamma;
  std::vector<double> vf(vax_frac.begin(), vax_frac.end());
  std::vector<double> x(n, 0.0), Jdummy;
  for (int a = 0; a < A; ++a){ x[a] = S0; x[A + a] = I0; }
  NumericMatrix inc(n_weeks, A);
  for (int t = 0; t < n_weeks; ++t){
    cm_week(x, f, beta, 7*t, vax_day, vf.data(), false, Jdummy);
    for (int a = 0; a < A; ++a) inc(t, a) = x[4*A + a];
  }
  NumericVector x_end(n), attack(A);
  for (int k = 0; k < n; ++k) x_end[k] = x[k];
  for (int a = 0; a < A; ++a) attack[a] = S0 - x[a] - x[2*A + a];
  return List::create(_["inc"] = inc, _["x_end"] = x_end, _["attack"] = attack);
}
