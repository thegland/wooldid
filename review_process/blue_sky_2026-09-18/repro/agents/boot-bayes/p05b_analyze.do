* p05b: GOF analysis of the multipliers already dumped by p05 (equal arm, R=1000)
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

mata:
void gofrun(string scalar Wname, string scalar Nname, string scalar am)
{
  real matrix W, P, CR
  real colvector nn, al, mu, vr, mhat, vhat, z, rs, v, u, i1, i0, off, sc, x, yy, cc
  real scalar R, C, N, a0, c, r, d, ksmax, kswhich, aa, bb
  W  = st_matrix(Wname)
  nn = st_matrix(Nname)
  R = rows(W); C = cols(W); N = sum(nn)
  if (am == "size") P = W :* nn'
  else              P = W
  rs = rowsum(P)
  printf("\n[%s] R=%f C=%f N=%f\n", am, R, C, N)
  printf("[%s] rowsum(p): maxdev from 1 = %g\n", am, max(abs(rs:-1)))
  if (am == "size") al = nn
  else              al = J(C,1,1)
  a0 = sum(al)
  mu = al/a0
  vr = (al:*(a0:-al)) :/ (a0^2*(a0+1))
  mhat = mean(P)'
  vhat = diagonal(variance(P))
  z = (mhat - mu) :/ sqrt(vr/R)
  printf("[%s] cluster means vs Dirichlet theory: max|z| = %6.3f  (|z|>3.5 suspicious)\n", am, max(abs(z)))
  printf("[%s] var ratio hat/theory: min=%6.4f max=%6.4f\n", am, min(vhat:/vr), max(vhat:/vr))
  ksmax = 0; kswhich = 0
  for (c=1; c<=C; c++) {
    v = sort(P[.,c], 1)
    u = J(R,1,.)
    for (r=1; r<=R; r++) u[r] = ibeta(al[c], a0-al[c], v[r])
    i1 = (1::R):/R
    i0 = (0::(R-1)):/R
    d = max((max(abs(u-i1)), max(abs(u-i0))))
    if (d > ksmax) { ksmax = d; kswhich = c }
  }
  printf("[%s] KS vs Beta(alpha_c, alpha0-alpha_c): max D = %6.4f at cluster %f  (5%% crit %6.4f, Bonf 5%%/C %6.4f)\n", am, ksmax, kswhich, 1.358/sqrt(R), 1.628/sqrt(R))
  CR = correlation(P)
  off = J(0,1,.)
  for (aa=1; aa<=C-1; aa++) for (bb=aa+1; bb<=C; bb++) off = off \ CR[aa,bb]
  printf("[%s] offdiag corr: mean=%8.5f   (Dirichlet(1..1) theory %8.5f)\n", am, mean(off), -1/(C-1))
  sc = J(C,1,.)
  for (c=1; c<=C; c++) {
    x  = P[1..(R-1), c]
    yy = P[2..R, c]
    cc = correlation((x, yy))
    sc[c] = cc[1,2]
  }
  printf("[%s] lag-1 serial corr across reps: max|rho| = %6.4f  (2/sqrt(R)=%6.4f)\n", am, max(abs(sc)), 2/sqrt(R))
  cc = correlation((nn, mhat))
  printf("[%s] corr(n_c, mean p_c) = %8.5f\n", am, cc[1,2])
}
end

use p05_panel, clear
qui bysort st (year): gen long nc = _N
qui bysort st (year): keep if _n==1
qui sort st
mkmat nc, matrix(NCraw)

use p05_equal, clear
qui keep if _rep >= 1
qui sort _rep
qui ds w_*
mkmat `r(varlist)', matrix(WF)
mata: gofrun("WF", "NCraw", "equal")

*--- now the size arm, 1000 reps -------------------------------------------
use p05_panel, clear
cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(size) reps(1000) rseed(20260917) bootsave(p05_size) replace dots(0)
dis "RUN_size_RC=" _rc
use p05_size, clear
qui keep if _rep >= 1
qui sort _rep
qui ds w_*
mkmat `r(varlist)', matrix(WS)
mata: gofrun("WS", "NCraw", "size")
dis "PROBE_DONE"
