* probe4_permgate.do -- GATE ARCHAEOLOGY on T1.
* The permute() nesting refusal (old ado:5990, now ado:6071-6081) is NOT removed: it
* is conditioned on `if "`clus'" != "`ivar'" & "`onc'" == ""'.  T1's own motivation
* was: "a user who follows both pieces of advice - coarse i, finer cluster - cannot
* use permute() at all."  This file asks whether the relief actually covers that user.
*
* permonsetcohorts requires a BIJECTION between i and the onset partition (ado:6002-6027):
* onset constant within i AND i constant within onset.  The second direction is the one
* that bites: as soon as TWO cohorts share an onset year -- the ordinary case in a
* staggered-adoption panel with 50 states and a handful of adoption years -- i is finer
* than the onset partition and the option is refused too.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 200
which wooldid

* 12 states, 4 counties each.  States 1-4 adopt in t=4, states 5-8 in t=6, 9-12 never.
* i = state (the help's "coarse i"), cluster = county (the help's "finer cluster").
prog def mkdat
  clear
  qui set obs 12
  gen long state = _n
  qui expand 4
  bysort state: gen byte cw = _n
  gen long county = _n
  gen int ttre = .
  qui replace ttre = 4 if state <= 4
  qui replace ttre = 6 if state > 4 & state <= 8
  gen double tau = 0.4 if !missing(ttre)
  qui replace tau = 0 if missing(tau)
  qui expand 8
  bysort county: gen int t = _n
  set seed 77
  qui bysort county (t): gen double mu = rnormal(0,0.4) if _n==1
  qui by county: replace mu = mu[1]
  gen double eps = rnormal(0,0.6)
  sort county t
  gen double y = tau*(t>=ttre & ttre<.) + 0.03*t + mu + eps
end
qui mkdat

dis "@@@ P1: i=state, cluster(county) finer, plain permute()"
cap noi wooldid y state t ttre, att att_i cluster(county) permute(shuffle) permreps(50) rseed(3) dots(0)
dis "@@@ P1_RC=" _rc

dis "@@@ P2: the same + permonsetcohorts (T1's remedy)"
cap noi wooldid y state t ttre, att att_i cluster(county) permute(shuffle) permreps(50) rseed(3) permonsetcohorts dots(0)
dis "@@@ P2_RC=" _rc

dis "@@@ P3: the SAME DESIGN with i = the onset grouping (one level per onset)"
qui mkdat
qui gen long cohid = cond(missing(ttre), 1, ttre)
cap noi wooldid y cohid t ttre, att att_i cluster(county) permute(shuffle) permreps(50) rseed(3) permonsetcohorts dots(0)
dis "@@@ P3_RC=" _rc

dis "@@@ P4: bootstrap(standard) on the P1 configuration that permute() refuses"
qui mkdat
cap noi wooldid y state t ttre, att att_i cluster(county) bootstrap(standard) reps(50) rseed(3) dots(0)
dis "@@@ P4_RC=" _rc
