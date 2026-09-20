*=======================================================================
* r6_df.do -- what a Wald/t built on an e(V) with N_V draws, but referred
* to e(df_r) = G-1, actually delivers.  Pure Monte Carlo, no wooldid, no
* data: the only wooldid inputs are the MEASURED pairs (N_V, df_r) from
* my own runs -- (6, 35) at reps(30), (4, 35) at reps(60) with six
* subgroups, plus (14, 29) from the original finding and (60, 35).
* SELF-CONTAINED, seconds.
*=======================================================================
clear all
set more off
set processors 2
set seed 20260917
mata:
real scalar wd_size(real scalar q, real scalar NV, real scalar dfr, real scalar nsim, real scalar lvl)
{
  real scalar i, rej, W, Fst, B
  real matrix Z, D, V, S
  real colvector b
  B = NV
  rej = 0
  S = I(q)                       // null: true covariance of the estimator
  for (i=1; i<=nsim; i++) {
     Z = rnormal(B, q, 0, 1)     // the B complete-case replication draws
     D = Z :- mean(Z)
     V = (D'*D)/(B-1)            // exactly wooldid_bootfin's listwise pass
     b = rnormal(q, 1, 0, 1)     // the observed estimate under a true null
     if (rank(V) < q) continue
     W = b' * invsym(V) * b      // what `test' forms
     Fst = W/q
     if (Ftail(q, dfr, Fst) < lvl) rej = rej + 1
  }
  return(rej/nsim)
}
real scalar wd_cover(real scalar NV, real scalar dfr, real scalar nsim, real scalar lvl)
{
  real scalar i, cov, se, b, tc
  real colvector z
  tc = invttail(dfr, lvl/2)
  cov = 0
  for (i=1; i<=nsim; i++) {
     z = rnormal(NV, 1, 0, 1)
     se = sqrt(variance(z))
     b  = rnormal(1, 1, 0, 1)
     if (abs(b) <= tc*se) cov = cov + 1
  }
  return(cov/nsim)
}
printf("--- joint Wald size at nominal 5%%, F(q, df_r) reference ---\n")
printf("  NV  dfr   q   actual_size\n")
NVl = (4, 6, 14, 30, 60)
dfl = (35, 35, 29, 35, 35)
for (a=1; a<=5; a++) {
  for (q=1; q<=6; q++) {
     if (q > NVl[a]-2) continue
     printf(" %3.0f %4.0f %3.0f   %6.4f\n", NVl[a], dfl[a], q, wd_size(q, NVl[a], dfl[a], 20000, .05))
  }
}
printf("--- single-coefficient CI coverage, nominal 95%%, t(df_r) x sqrt(V_kk) ---\n")
printf("  NV  dfr   coverage\n")
for (a=1; a<=5; a++) {
  printf(" %3.0f %4.0f   %6.4f\n", NVl[a], dfl[a], wd_cover(NVl[a], dfl[a], 40000, .05))
}
end
dis "R6_DONE"
