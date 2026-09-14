# contact_matrix.R -- contact structure for the compartmental model (ASSUMPTIONS.md B3, C3)
#
# Two jobs, both base R:
#   1. collapse the 4-group Prem-derived matrices in models_in$contacts (0-4 | 5-14 | 15-64 | 65+,
#      mean contacts per person per day, rows = the contacting group) into the model's three groups
#      young 0-14 | medium 15-64 | elderly 65+, under the directed-contact-ends convention;
#   2. scale the matrix so that the force of infection realises EXACTLY the intended R0 at full
#      susceptibility, given the age mixing.
#
# THE R0 CALIBRATION RULE (owner decision 2026-09; the clear note the owner asked for):
#   With the force of infection  lambda_a = beta * sum_j Cn[a,j] * I_j  (I_j = prevalence FRACTION of
#   group j) the next-generation matrix at S = 1 is K = (beta/gamma) * Cn = R0 * Cn, so the realised
#   basic reproduction number is  R0 * rho(Cn)  with rho the dominant eigenvalue. Setting
#   Cn = C / rho(C) makes it exactly R0. This is coordinate-free: rho(C) == rho(D C D^-1) for the
#   count-based matrix, so it does not matter whether compartments are fractions or counts.
#   The Stan model instead used beta * a_factor[a] * rowNormalised(C), which equals (beta / cbar) * C
#   with cbar the population-weighted mean contacts per person; since rho(C) >= cbar for any
#   assortative matrix, its realised R0 was 1.5 * rho(C)/cbar = 1.58-1.73 across the 28 country
#   matrices (median 1.65), never 1.5, and DIFFERENT BY COUNTRY -- a contact-data artefact that would
#   read as between-country transmissibility differences. method = "stan_cbar" is kept only so tests
#   can reproduce that number; the model uses "spectral_radius".

# ---- |-collapse a contact matrix to coarser groups via directed contact-ends ----
# C: k x k mean contacts per person (row = contacting group); N: the k group populations;
# groups: list of integer index vectors (one per new group). Ends E[a,j] = C[a,j] * N[a] add over
# blocks; the merged per-person rate is (ends in block) / (population of the new row group).
collapse_contact_groups = function(C, N, groups){
  C = as.matrix(C); stopifnot(nrow(C) == ncol(C), length(N) == nrow(C))
  E  = C * N                                   # row a scaled by N[a]: directed ends a -> j
  Ng = vapply(groups, function(i) sum(N[i]), numeric(1))
  Cg = matrix(0, length(groups), length(groups), dimnames = list(names(groups), names(groups)))
  for (a in seq_along(groups)) for (j in seq_along(groups))
    Cg[a, j] = sum(E[groups[[a]], groups[[j]], drop = FALSE]) / Ng[a]
  list(C = Cg, N = Ng)
}

# ---- |-dominant eigenvalue of a non-negative matrix ----
spectral_radius = function(C) max(Re(eigen(as.matrix(C), only.values = TRUE)$values))

# ---- |-scale the matrix so the force of infection realises R0 exactly (or reproduce the old rule) ----
# Returns the scaled matrix with attribute "scaling" documenting what was applied.
contact_matrix_normalised = function(C, N = NULL, method = c("spectral_radius", "stan_cbar")){
  method = match.arg(method); C = as.matrix(C)
  if (method == "spectral_radius"){
    Cn = C / spectral_radius(C)
    stopifnot(abs(spectral_radius(Cn) - 1) < 1e-10)      # the assertion the settings file promises
  } else {                                                # the Stan model's convention, for tests only
    stopifnot(!is.null(N))
    cbar = sum(N * rowSums(C)) / sum(N)
    Cn = C / cbar
  }
  attr(Cn, "scaling") = method
  Cn
}

# ---- |-realised R0 at full susceptibility for a scaled matrix (the NGM spectral radius) ----
# lambda_a = R0 * gamma * sum_j Cn[a,j] I_j and recovery gamma give K = R0 * Cn at S = 1, so the
# realised R0 is R0 * rho(Cn); gamma cancels. Exposed so tests and reports can state it.
realised_R0 = function(Cn, R0) R0 * spectral_radius(Cn)

# ---- |-country contact matrix for the model: 3 groups, scaled ----
# models_in: output/models_in.rds (contacts = 4-group matrices by country name); pop4: the 4-group
# population vector in the order 0-4, 5-14, 15-64, 65+; settings: comp_model_settings().
model_contact_matrix = function(C4, pop4, settings){
  groups = list(young = 1:2, medium = 3, elderly = 4)     # 0-4 + 5-14 | 15-64 | 65+ (the 4-group cut nests the 17 bands)
  cg = collapse_contact_groups(C4, pop4, groups)
  list(C = contact_matrix_normalised(cg$C, cg$N, settings$contact_scaling), N = cg$N)
}
