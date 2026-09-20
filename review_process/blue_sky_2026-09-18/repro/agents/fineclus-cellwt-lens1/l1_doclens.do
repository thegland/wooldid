*! l1 — DOCUMENTED-BEHAVIOUR lens on FCW-02.
*! Not a re-measurement of the magnitude (prober + skeptic already paired that).
*! Three lens-specific questions:
*!   Q1  does the flip reproduce independently on one A/D pair?                (yes/no)
*!   Q2  is it SILENT — does the printed output of arm D mention the regime,
*!       and does e() record which regime ran?                                (silence?)
*!   Q3  the same phantom row on the PERMUTE side: is it loud there?          (r459?)
*! Q3 is the documentary contrast: sthlp:431 documents design-vs-markout for permute()
*! and the consequence there is a refusal; sthlp:241 states the bootstrap regime rule
*! with no sample qualifier and the consequence there is a silently smaller SE.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkpanel
prog def mkpanel
  syntax , dseed(int)
  local G 15
  local U 4
  local T 7
  clear
  set seed `dseed'
  set sortseed `dseed'
  set obs `G'
  gen long cty = _n
  gen double a_g = rnormal(0,.4)
  gen onset = .
  replace onset = 4 if cty <= 5
  replace onset = 6 if cty > 5 & cty <= 10
  expand `U'
  sort cty
  gen long unit = _n
  gen double b_i = rnormal(0,.4)
  gen double uf = rnormal(0,1)
  expand `T'
  bysort unit (cty): gen int t = _n
  gen double tau = 1 + a_g + b_i
  gen byte D = (onset < . & t >= onset)
  gen double y = 2 + uf + .1*t + tau*D + rnormal(0,1)
  * THE PHANTOM ROW: unit 1 at t = 3 has a missing outcome in BOTH arms, so the
  * estimator provably discards it either way.
  replace y = . if unit == 1 & t == 3
  keep y unit cty t onset
end

************ ARM A : phantom row keeps its own county. unit 1 nests in cty 1.
mkpanel , dseed(7201)
qui count if unit == 1 & t == 3
dis "A_PHANTOMROWS=" r(N)
log using "l1_armA.log", replace text name(armA)
cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(199) rseed(41) dots(0)
log close armA
dis "A_RC=" _rc " A_N=" e(N)
mat bA = e(b)
mat VA = e(V)
dis "A_b_att   = " %21x bA[1,1]
dis "A_b_atti  = " %21x bA[1,2]
dis "A_se_att  = " %21x sqrt(VA[1,1])
dis "A_se_atti = " %21x sqrt(VA[2,2])
qui est store eA

************ ARM D : the SAME unestimable row is recorded under county 2.
mkpanel , dseed(7201)
qui replace cty = 2 if unit == 1 & t == 3
qui count if unit == 1 & t == 3 & cty == 2 & missing(y)
dis "D_PHANTOMROWS=" r(N)
log using "l1_armD.log", replace text name(armD)
cap noi wooldid y unit t onset, att att_i cluster(cty) bootstrap(standard) reps(199) rseed(41) dots(0)
log close armD
dis "D_RC=" _rc " D_N=" e(N)
mat bD = e(b)
mat VD = e(V)
dis "D_b_att   = " %21x bD[1,1]
dis "D_b_atti  = " %21x bD[1,2]
dis "D_se_att  = " %21x sqrt(VD[1,1])
dis "D_se_atti = " %21x sqrt(VD[2,2])
dis "RATIO_se_atti = " %12.8f sqrt(VD[2,2])/sqrt(VA[2,2])
dis "RATIO_se_att  = " %12.8f sqrt(VD[1,1])/sqrt(VA[1,1])
dis "DELTA_b_atti  = " %12.0g (bD[1,2]-bA[1,2])

************ Q2 : does e() record the regime anywhere?
dis "--- EBEGIN ---"
qui ereturn list
ereturn list
dis "--- EEND ---"

************ Q3 : the same phantom row under permute(). Loud, or silent?
dis "--- PERMUTE ON THE D DATA ---"
cap noi wooldid y unit t onset, att att_i cluster(cty) permute(shuffle) permreps(19) rseed(41) dots(0)
dis "PERM_D_RC=" _rc

************ Q3b : control — permute() on the A data (no phantom straddler)
mkpanel , dseed(7201)
dis "--- PERMUTE ON THE A DATA ---"
cap noi wooldid y unit t onset, att att_i cluster(cty) permute(shuffle) permreps(19) rseed(41) dots(0)
dis "PERM_A_RC=" _rc
dis "DONE_L1"
