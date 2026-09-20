* ============================================================================
* SELF-CONTAINED REPRO -- gteffects cell standard errors that are numerically
* zero are reported as REAL standard errors, giving p-values of order 1e-45 and
* confidence intervals of width 1e-14.
*
* The gteffects help (wooldid.sthlp, the gteffects option paragraph) states:
*   "... degenerate zero variances are always reported as missing, never zero"
* The guard that implements this (wooldid.ado:13496, inside mata wooldid_gtx) is
*   se = (v < . & v > 0) ? sqrt(v) : .
* i.e. EXACT zero only. A rank-deficient cluster-robust VCE leaves floating-point
* residue of order 1e-30 on the affected diagonal entries, which passes the guard.
*
* Trigger: cluster() at the level of a COARSE i (each i is exactly one cluster).
* Run:  stata-mp -e do repro_gtse_degenerate.do   then read repro_gtse_degenerate.log
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

* 40 states, 20 years, three adoption waves + never-treated. i = onset cohort,
* cluster at that same level. Clean, well-behaved, fully balanced data.
clear
set obs 40
gen long state = _n
gen byte wave  = mod(_n-1, 4)
gen double ttre = .
replace ttre = 2006 if wave==1
replace ttre = 2010 if wave==2
replace ttre = 2014 if wave==3
gen long cohid = wave + 1
gen double sfe = rnormal()*2
expand 20
bysort state: gen int year = 1999 + _n
gen double tfe = 0.05*(year-1999)
gen byte trt = (ttre<.) & (year >= ttre)
gen double rel = year - ttre
gen double y = 10 + sfe + tfe + cond(trt==1, 1.5 + 0.10*rel, 0) + rnormal()

dis "########## ARM A: cluster AT the level of i (cohid) ##########"
cap noi wooldid y cohid year ttre, att esprelength(3) espostlength(3) gteffects cluster(cohid)
dis "ARM_A_RC=" _rc
matrix GA = e(gteffects)
mata:
  G  = st_matrix("GA")
  cn = st_matrixcolstripe("GA")[.,2]
  js = selectindex(cn :== "se"); jp = selectindex(cn :== "p")
  jb = selectindex(cn :== "estimate"); ji = selectindex(cn :== "i")
  jr = selectindex(cn :== "relyr")
  jl = selectindex(cn :== "lb_95ci"); ju = selectindex(cn :== "ub_95ci")
  se = G[.,js]
  printf("A|rows=%g  se_missing=%g  se_exactly_zero=%g  se_below_1e-8=%g\n", ///
         rows(se), sum(se :>= .), sum(se :== 0), sum((se :< 1e-8) :& (se :> 0)))
  bad = select(G, (G[.,js] :< 1e-8) :& (G[.,js] :> 0))
  for (r=1; r<=rows(bad); r++) {
    printf("A_BAD|i=%g relyr=%g  b=%12.8f  se=%10.3e  p=%10.3e  ci=[%12.8f , %12.8f]  width=%10.3e\n", ///
      bad[r,ji], bad[r,jr], bad[r,jb], bad[r,js], bad[r,jp], bad[r,jl], bad[r,ju], bad[r,ju]-bad[r,jl])
  }
end

dis "########## ARM B: same data, cluster BELOW i (state) -- control ##########"
cap noi wooldid y cohid year ttre, att esprelength(3) espostlength(3) gteffects cluster(state)
dis "ARM_B_RC=" _rc
matrix GB = e(gteffects)
mata:
  G = st_matrix("GB"); cn = st_matrixcolstripe("GB")[.,2]
  js = selectindex(cn :== "se"); se = G[.,js]
  printf("B|rows=%g  se_missing=%g  se_exactly_zero=%g  se_below_1e-8=%g  min_se=%10.3e\n", ///
         rows(se), sum(se :>= .), sum(se :== 0), sum((se :< 1e-8) :& (se :> 0)), min(select(se, se :< .)))
end

dis "########## ARM C: same data, robust -- control ##########"
cap noi wooldid y cohid year ttre, att esprelength(3) espostlength(3) gteffects robust
dis "ARM_C_RC=" _rc
matrix GC = e(gteffects)
mata:
  G = st_matrix("GC"); cn = st_matrixcolstripe("GC")[.,2]
  js = selectindex(cn :== "se"); se = G[.,js]
  printf("C|rows=%g  se_missing=%g  se_exactly_zero=%g  se_below_1e-8=%g  min_se=%10.3e\n", ///
         rows(se), sum(se :>= .), sum(se :== 0), sum((se :< 1e-8) :& (se :> 0)), min(select(se, se :< .)))
end

dis "########## VERDICT ##########"
dis "Arm A prints ordinary-looking rows whose se is ~1e-15, whose p is ~1e-45, and"
dis "whose 95% CI has width ~1e-14.  The help promises those rows carry a MISSING se."
dis "===== REPRO DONE ====="
