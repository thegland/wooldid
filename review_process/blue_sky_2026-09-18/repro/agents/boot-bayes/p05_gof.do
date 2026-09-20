* p05: distributional goodness-of-fit on the ACTUAL multipliers wooldid applied.
*   equal: p_c ~ Beta(1, C-1),  sum_c p_c = 1
*   size : p_c = n_c*w_c ~ Beta(n_c, N-n_c), sum_c n_c*w_c = 1
* Tests: row sums, per-cluster mean/var vs theory, KS vs the Beta marginal,
*        cross-cluster correlation vs -1/(C-1), lag-1 serial correlation,
*        exchangeability (cluster means equal), size-independence under equal.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-bayes"

*--- small unbalanced 12-cluster panel --------------------------------------
clear
set obs 12
gen int st = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 2005 if st <= 3
replace onset = 2008 if st > 3 & st <= 6
replace onset = 2011 if st > 6 & st <= 8
expand 14
bysort st: gen int year = 2000 + _n
* deliberately unequal cluster row counts: 6..14 rows
gen byte keepn = _n
bysort st: drop if _n > (6 + st/2)
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.03*(year-2000) + 1.0*post + rnormal()*0.5
drop keepn
save p05_panel, replace
qui bysort st: gen long nc = _N if _n==1
mkmat nc if !missing(nc), matrix(NCraw)
qui su nc
dis "C=" r(N) "  Nrows=" _N "  minn=" r(min) " maxn=" r(max)
drop nc

local R 1000
foreach am in equal size {
  use p05_panel, clear
  cap noi wooldid y st year onset, att cluster(st) bootstrap(bayes) bootalpha(`am') reps(`R') rseed(20260917) bootsave(p05_`am') replace dots(0)
  dis "RUN_`am'_RC=" _rc

  use p05_`am', clear
  qui keep if _rep >= 1
  qui sort _rep
  qui ds w_*
  local wv `r(varlist)'
  mkmat `wv', matrix(WF)
  mata {
    W  = st_matrix("WF")            // R x C
    nn = st_matrix("NCraw")         // C x 1
    R  = rows(W); C = cols(W); N = sum(nn)
    am = "`am'"
    // p_c: the Dirichlet mass itself
    P = (am=="size" ? W :* nn' : W)
    rs = rowsum(P)
    printf("\n[`am'] R=%f C=%f N=%f\n", R, C, N)
    printf("[`am'] rowsum(p): min=%21x max=%21x  maxdev=%g\n", min(rs), max(rs), max(abs(rs:-1)))
    // theory
    al = (am=="size" ? nn : J(C,1,1))
    a0 = sum(al)
    mu = al/a0
    vr = (al:*(a0:-al)) :/ (a0^2*(a0+1))
    mhat = mean(P)'
    vhat = diagonal(variance(P))
    z = (mhat - mu) :/ sqrt(vr/R)
    printf("[`am'] cluster means vs theory: max|z| = %g  (|z|>3.5 suspicious)\n", max(abs(z)))
    printf("[`am'] var ratio hat/theory: min=%g max=%g\n", min(vhat:/vr), max(vhat:/vr))
    // KS per cluster against Beta(al_c, a0-al_c)
    ksmax = 0; kswhich = 0
    for (c=1; c<=C; c++) {
      v = sort(P[.,c], 1)
      u = J(R,1,.)
      for (r=1; r<=R; r++) u[r] = ibeta(al[c], a0-al[c], v[r])
      i1 = (1::R)/R
      i0 = (0::(R-1))/R
      d = max((max(abs(u-i1)), max(abs(u-i0))))
      if (d > ksmax) { ksmax = d; kswhich = c }
    }
    printf("[`am'] KS vs Beta marginal: max D = %g at cluster %f   (5%% crit ~ %g, Bonferroni 5%%/C ~ %g)\n", ksmax, kswhich, 1.358/sqrt(R), 1.628/sqrt(R))
    // cross-cluster correlation
    CR = correlation(P)
    off = J(0,1,.)
    for (a=1; a<=C-1; a++) for (b=a+1; b<=C; b++) off = off \ CR[a,b]
    printf("[`am'] offdiag corr: mean=%g  (Dirichlet theory for equal = %g)\n", mean(off), -1/(C-1))
    // lag-1 serial correlation across replications
    sc = J(C,1,.)
    for (c=1; c<=C; c++) {
      x = P[1..(R-1), c]; yy = P[2..R, c]
      sc[c] = correlation((x, yy))[1,2]
    }
    printf("[`am'] lag-1 serial corr across reps: max|rho| = %g  (2/sqrt(R) = %g)\n", max(abs(sc)), 2/sqrt(R))
    // does the EQUAL marginal depend on cluster size? (it must not)
    printf("[`am'] corr(n_c, mean p_c) = %g\n", correlation((nn, mhat))[1,2])
  }
}
dis "PROBE_DONE"
