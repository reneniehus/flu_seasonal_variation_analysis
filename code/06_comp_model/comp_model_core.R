# comp_model_core.R -- the compartmental model's engine (base R, the REFERENCE implementation)
#
# Implements code/06_comp_model/ASSUMPTIONS.md sections 1, 3, 4, 5, 6, 7 exactly:
#   state per age group a (fractions of the group): S_u, I_u, S_v, I_v and the weekly accumulator C
#   of observation-relevant new infections (unvaccinated + (1 - ve_ili_cond_inf) x vaccinated);
#   R_u, R_v are implied (they feed neither transmission nor observation) and are not carried.
#   Daily forward Euler; force of infection lambda_a = beta * sum_j Cn[a,j] * (I_u[j] + (1-ve_spread) I_v[j])
#   with Cn scaled to spectral radius 1 (contact_matrix.R) so beta = R0 * gamma realises R0 exactly;
#   vaccination as a single pulse (fraction v of S_u moved to S_v on the pulse day); seasons reset.
#   Observation: mu[t,a] = c * N_a * C_a(t) + b * N_a / 1e5 on a COUNT scale, Var = mu + mu^2/phi.
#   EKF with multiplicative process noise on the infected fractions (sd = q * I).
# The Jacobians are ANALYTICAL (finite differences are ~12x slower and this is the inner loop).
# The C++ port (stage iii) must reproduce every function here to 1e-10; see ASSUMPTIONS.md I2.

# ---- |-state layout: 5 blocks of n_age entries ----
cm_layout = function(n_age){
  blk = function(k) (k - 1) * n_age + seq_len(n_age)
  list(n_age = n_age, n = 5 * n_age, S_u = blk(1), I_u = blk(2), S_v = blk(3), I_v = blk(4), C = blk(5))
}

# ---- |-fixed inputs of one country (everything the dynamics need besides the fitted parameters) ----
# Cn: contact matrix with spectral radius 1; N: group populations; the rest from comp_model_settings().
cm_fixed = function(Cn, N, settings){
  n_age = nrow(Cn)
  stopifnot(length(N) == n_age, abs(spectral_radius(Cn) - 1) < 1e-8)
  list(n_age = n_age, L = cm_layout(n_age), Cn = as.matrix(Cn), N = as.numeric(N),
       gamma = settings$gamma_per_day, ve_inf = settings$ve_inf, ve_spread = settings$ve_spread,
       ve_ili = settings$ve_ili_cond_inf, I0 = settings$I0_fraction, p0 = settings$p0_frac,
       rate_per = settings$rate_per)
}

# ---- |-one day of forward Euler (dt = 1 day), optionally with the analytical Jacobian ----
# vax_frac: per-age fraction of S_u moved to S_v AFTER the infection/recovery update (0 = no pulse).
cm_day = function(x, f, beta, vax_frac = NULL, jac = FALSE){
  L = f$L; A = f$n_age; Cn = f$Cn
  S_u = x[L$S_u]; I_u = x[L$I_u]; S_v = x[L$S_v]; I_v = x[L$I_v]; C = x[L$C]
  s = 1 - f$ve_spread; e = 1 - f$ve_inf; w = 1 - f$ve_ili; g = f$gamma
  lam   = as.numeric(beta * (Cn %*% (I_u + s * I_v)))     # force of infection per age group
  new_u = lam * S_u
  new_v = e * lam * S_v
  xn = x
  xn[L$S_u] = S_u - new_u
  xn[L$I_u] = I_u + new_u - g * I_u
  xn[L$S_v] = S_v - new_v
  xn[L$I_v] = I_v + new_v - g * I_v
  xn[L$C]   = C + new_u + w * new_v
  if (!is.null(vax_frac)){                                # vaccination pulse (infectious are not vaccinated)
    moved = vax_frac * xn[L$S_u]
    xn[L$S_u] = xn[L$S_u] - moved
    xn[L$S_v] = xn[L$S_v] + moved
  }
  xn[c(L$S_u, L$I_u, L$S_v, L$I_v)] = pmin(pmax(xn[c(L$S_u, L$I_u, L$S_v, L$I_v)], 0), 1)
  if (!jac) return(xn)

  # analytical Jacobian d(xn)/d(x) of the update above (before clamping)
  J = diag(L$n)
  dlam_dIu = beta * Cn; dlam_dIv = beta * s * Cn          # A x A
  Dl = diag(lam, A); DSu = diag(S_u, A); DSv = diag(S_v, A); I_A = diag(A)
  # S_u' = S_u - lam*S_u
  J[L$S_u, L$S_u] = I_A - Dl;            J[L$S_u, L$I_u] = -DSu %*% dlam_dIu;      J[L$S_u, L$I_v] = -DSu %*% dlam_dIv
  # I_u' = I_u + lam*S_u - g*I_u
  J[L$I_u, L$S_u] = Dl;                  J[L$I_u, L$I_u] = (1 - g) * I_A + DSu %*% dlam_dIu;  J[L$I_u, L$I_v] = DSu %*% dlam_dIv
  # S_v' = S_v - e*lam*S_v
  J[L$S_v, L$S_v] = I_A - e * Dl;        J[L$S_v, L$I_u] = -e * DSv %*% dlam_dIu;  J[L$S_v, L$I_v] = -e * DSv %*% dlam_dIv
  # I_v' = I_v + e*lam*S_v - g*I_v
  J[L$I_v, L$S_v] = e * Dl;              J[L$I_v, L$I_u] = e * DSv %*% dlam_dIu;   J[L$I_v, L$I_v] = (1 - g) * I_A + e * DSv %*% dlam_dIv
  # C' = C + lam*S_u + w*e*lam*S_v
  J[L$C, L$S_u] = Dl;                    J[L$C, L$S_v] = w * e * Dl
  J[L$C, L$I_u] = DSu %*% dlam_dIu + w * e * DSv %*% dlam_dIu
  J[L$C, L$I_v] = DSu %*% dlam_dIv + w * e * DSv %*% dlam_dIv
  if (!is.null(vax_frac)){                                # linear pulse: S_u' = (1-v) S_u' ; S_v' = S_v' + v S_u'
    V = diag(vax_frac, A)
    Jv = diag(L$n)
    Jv[L$S_u, L$S_u] = I_A - V
    Jv[L$S_v, L$S_u] = V
    J = Jv %*% J
  }
  list(x = xn, J = J)
}

# ---- |-one observation week: reset C, integrate 7 days; pulse applied on its day; Jacobian optional ----
# day0: the season day index of the day BEFORE this week (week t covers days day0+1 .. day0+7).
cm_week = function(x, f, beta, day0, vax_day = NA, vax_frac = NULL, jac = FALSE){
  L = f$L
  x[L$C] = 0
  J = if (jac) diag(L$n) else NULL
  if (jac) J[L$C, ] = 0                                  # C is reset: it carries no memory of x
  for (d in 1:7){
    day = day0 + d
    vf = if (!is.na(vax_day) && day == vax_day) vax_frac else NULL
    if (jac){ st = cm_day(x, f, beta, vf, jac = TRUE); x = st$x; J = st$J %*% J }
    else      x = cm_day(x, f, beta, vf)
  }
  if (jac) list(x = x, J = J) else x
}

# ---- |-initial state of a season ----
# I0: the seed (infected fraction of every age group on season day 1). Fixed at settings$I0_fraction
# by default; when settings$I0_by_season is TRUE the fit supplies a per-season value, which is the
# smooth way to set each season's ARRIVAL TIME (a seed k times smaller arrives log(k)/r days later
# at growth rate r) without touching the growth rate itself (R0_s * S0_c). See ASSUMPTIONS.md E3.
cm_init = function(f, S0, I0 = f$I0){
  L = f$L; x = numeric(L$n)
  x[L$S_u] = S0; x[L$I_u] = I0                          # implied R_u = 1 - S0 - I0
  x
}

# ---- |-deterministic season: weekly expected observation-relevant incidence (fraction of group) ----
# Returns list(inc = n_weeks x A matrix of C per week, S_end, attack = per-age cumulative infections).
cm_simulate_season = function(f, S0, R0, n_weeks, vax_day = NA, vax_frac = NULL, I0 = f$I0){
  L = f$L; beta = R0 * f$gamma; x = cm_init(f, S0, I0)
  inc = matrix(NA_real_, n_weeks, f$n_age)
  for (t in seq_len(n_weeks)){
    x = cm_week(x, f, beta, day0 = 7 * (t - 1), vax_day, vax_frac)
    inc[t, ] = x[L$C]
  }
  list(inc = inc, x_end = x, attack = S0 - x[L$S_u] - x[L$S_v])   # S only leaves S_u/S_v by infection
}

# ---- |-expected observed counts from incidence fractions ----
# c: the reporting proportion, a scalar (age-invariant, F2) or a vector of length n_age (c_by_age)
cm_mu = function(inc, f, c, b) sweep(inc, 2, rep_len(c, f$n_age) * f$N, "*") + matrix(b * f$N / f$rate_per, nrow(inc), f$n_age, byrow = TRUE)

# ---- |-EKF over one season: returns the innovation log-likelihood + filtered/predicted quantities ----
# y: n_weeks x A matrix of observed COUNTS (NA = not observed). q: multiplicative process-noise sd on I.
cm_ekf_season = function(y, f, S0, R0, c, b, phi, q, vax_day = NA, vax_frac = NULL, I0 = f$I0){
  L = f$L; A = f$n_age; n_weeks = nrow(y); beta = R0 * f$gamma
  x = cm_init(f, S0, I0)
  P = matrix(0, L$n, L$n)
  P[cbind(L$S_u, L$S_u)] = (f$p0 * S0)^2                # tight initial covariance: trust the fitted S0, I0
  P[cbind(L$I_u, L$I_u)] = (f$p0 * I0)^2
  cA = rep_len(c, A)                                    # scalar or per-age reporting proportion
  Hbase = matrix(0, A, L$n); Hbase[cbind(seq_len(A), L$C)] = cA * f$N   # mu = c_a N_a C_a + b N_a/rate_per
  b_cnt = b * f$N / f$rate_per
  ll = 0; mu_pred = matrix(NA_real_, n_weeks, A); I_filt = matrix(NA_real_, n_weeks, A); S_filt = I_filt
  for (t in seq_len(n_weeks)){
    st = cm_week(x, f, beta, day0 = 7 * (t - 1), vax_day, vax_frac, jac = TRUE)
    xpr = st$x; Ppr = st$J %*% P %*% t(st$J)
    Ipr = xpr[c(L$I_u, L$I_v)]
    Ppr[cbind(c(L$I_u, L$I_v), c(L$I_u, L$I_v))] = Ppr[cbind(c(L$I_u, L$I_v), c(L$I_u, L$I_v))] + (q * Ipr)^2   # multiplicative noise on I
    mu = as.numeric(Hbase %*% xpr) + b_cnt
    mu_pred[t, ] = mu
    ok = which(is.finite(y[t, ]))
    if (length(ok)){
      H = Hbase[ok, , drop = FALSE]; Rt = mu[ok] + mu[ok]^2 / phi
      S = H %*% Ppr %*% t(H) + diag(Rt, length(ok))
      S = (S + t(S)) / 2
      cS = tryCatch(chol(S), error = function(e) NULL)
      if (is.null(cS)) return(list(loglik = -1e10, mu_pred = mu_pred, I = I_filt, S = S_filt))
      innov = y[t, ok] - mu[ok]
      Sinv_innov = backsolve(cS, forwardsolve(t(cS), innov))
      K = Ppr %*% t(H) %*% chol2inv(cS)
      x = as.numeric(xpr + K %*% innov)
      P = (diag(L$n) - K %*% H) %*% Ppr
      ll = ll - 0.5 * (2 * sum(log(diag(cS))) + sum(innov * Sinv_innov) + length(ok) * log(2 * pi))
    } else { x = xpr; P = Ppr }
    # clamp the UPDATED state to [0, 1] exactly as the deterministic step does. A positive floor (an
    # earlier 1e-12) would re-seed every compartment every week -- I_v before the vaccination pulse
    # in particular -- and make the EKF drift off the deterministic model by ~1e-12/I0 per week
    x[c(L$S_u, L$I_u, L$S_v, L$I_v)] = pmin(pmax(x[c(L$S_u, L$I_u, L$S_v, L$I_v)], 0), 1)
    I_filt[t, ] = x[L$I_u] + x[L$I_v]; S_filt[t, ] = x[L$S_u] + x[L$S_v]
  }
  list(loglik = ll, mu_pred = mu_pred, I = I_filt, S = S_filt)
}
