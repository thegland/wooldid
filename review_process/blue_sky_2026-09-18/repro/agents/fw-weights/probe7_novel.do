* probe7: probes nobody suggested, all of them metamorphic identities on weights.
*  Z  ZERO WEIGHT == DELETED ROW.  An fweight of 0 means "zero copies".  A user who
*     zeroes weights instead of dropping rows must get EXACTLY the dropped-row answer,
*     including N, df and SE, because wooldid marks out w<=0 itself.
*  C  customweightmultiplier() SCALE.  cwm enters the margins weight only; a CONSTANT
*     multiplier cancels in every weighted average, so cwm(1), cwm(2) and no cwm must
*     give identical estimates.
*  P  POISSON + pweight vs physical duplication (a completely different estimator,
*     ppmlhdfe; aweights are refused there so pw is the only weighted poisson path).
*  V  does pweight actually force a robust VCE, as sthlp:268 states?
*  S  more surface for the [aw=k] == duplication identity: unconditionalse, impvalues,
*     anticipation, oostreatedcontrols, esfixedbaseperiod, esrelativeto, timetrends.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

mata:
void wd_rel(string scalar a, string scalar b, string scalar tag) {
  real matrix A, B, D
  A = st_matrix(a); B = st_matrix(b)
  if (rows(A)==0 | rows(B)==0) {
    printf("REL_%s=EMPTY\n", tag)
  }
  else if (rows(A)!=rows(B) | cols(A)!=cols(B)) {
    printf("REL_%s=SHAPEMISMATCH %g x %g vs %g x %g\n", tag, rows(A),cols(A),rows(B),cols(B))
  }
  else {
    D = abs(A:-B) :/ (abs(A):*(abs(A):>1e-12) :+ (abs(A):<=1e-12))
    printf("REL_%s=%10.3e  maxabs=%10.3e\n", tag, colmax(rowmax(D)), colmax(rowmax(abs(A:-B))))
  }
}
end

program define pairrun
  syntax , fa(string) ca(string asis) fb(string) cb(string asis) tag(string)
  tempname b1 b2 m1 m2 v1 v2
  qui use "`fa'", clear
  cap noi `ca'
  local r1 = _rc
  if `r1' == 0 {
    matrix `b1' = e(b)
    cap matrix `m1' = e(wdidmainresults)
    cap matrix `v1' = e(V)
    local n1 = e(N)
    local d1 = e(df_r)
  }
  qui use "`fb'", clear
  cap noi `cb'
  local r2 = _rc
  if `r2' == 0 {
    matrix `b2' = e(b)
    cap matrix `m2' = e(wdidmainresults)
    cap matrix `v2' = e(V)
    local n2 = e(N)
    local d2 = e(df_r)
  }
  dis "PR_`tag'_RC=" `r1' "," `r2'
  if `r1' == 0 & `r2' == 0 {
    dis "PR_`tag'_N=" `n1' "," `n2' " DF=" `d1' "," `d2'
    mata: wd_rel("`b1'","`b2'","`tag'_b")
    mata: wd_rel("`m1'","`m2'","`tag'_main")
    mata: wd_rel("`v1'","`v2'","`tag'_V")
  }
end

*--- build the zero-weight / deleted-row pair ---------------------------------
use "p2_A.dta", clear
gen byte keeprow = !(mod(_n,7) == 0)
gen double wz = keeprow            /* 0 for the rows the other frame deletes */
gen double wzk = k*keeprow
save "p7_zero.dta", replace
keep if keeprow == 1
save "p7_del.dta", replace

*--- Z: zero weight must equal deleting the row -------------------------------
pairrun , fa(p7_zero.dta) ca(wooldid y i t ttre [aw=wz], att att_i att_it att_itime) fb(p7_del.dta) cb(wooldid y i t ttre, att att_i att_it att_itime) tag(ZERO)
pairrun , fa(p7_zero.dta) ca(wooldid y i t ttre [aw=wzk], att att_i att_it att_itime) fb(p7_del.dta) cb(wooldid y i t ttre [aw=k], att att_i att_it att_itime) tag(ZEROK)
pairrun , fa(p7_zero.dta) ca(wooldid y i t ttre [aw=wz], att att_it esprelength(2) espostlength(3) cluster(i)) fb(p7_del.dta) cb(wooldid y i t ttre, att att_it esprelength(2) espostlength(3) cluster(i)) tag(ZEROES)
* the same thing expressed as an `if' - the documented route
pairrun , fa(p7_zero.dta) ca(wooldid y i t ttre if keeprow == 1, att att_i att_it att_itime) fb(p7_del.dta) cb(wooldid y i t ttre, att att_i att_it att_itime) tag(ZEROIF)

*--- C: a constant customweightmultiplier() must be a no-op --------------------
use "p2_A.dta", clear
gen double two = 2
gen double c1  = 1
save "p7_cwm.dta", replace
pairrun , fa(p7_cwm.dta) ca(wooldid y i t ttre [aw=k], att att_i att_it att_itime customweightmultiplier(c1)) fb(p7_cwm.dta) cb(wooldid y i t ttre [aw=k], att att_i att_it att_itime customweightmultiplier(two)) tag(CWMSCALE)
pairrun , fa(p7_cwm.dta) ca(wooldid y i t ttre [aw=k], att att_i att_it att_itime) fb(p7_cwm.dta) cb(wooldid y i t ttre [aw=k], att att_i att_it att_itime customweightmultiplier(c1)) tag(CWMNONE)
* cwm == the weight ratio: [aw=k] with cwm(k) must hit the same margins weights as [aw=k*k]
*   for the att estimand only (att weight is w*cwm); the regression weight differs, so only
*   the ESTIMAND is comparable - recorded, not asserted.
pairrun , fa(p7_cwm.dta) ca(wooldid y i t ttre [aw=k], att customweightmultiplier(two)) fb(p7_cwm.dta) cb(wooldid y i t ttre [aw=k], att) tag(CWMTWO)

*--- P: poisson + pweight vs physical duplication ------------------------------
use "p2_A.dta", clear
gen double ypos = exp(0.3 + 0.2*post + 0.1*t + 0.3*ufe)
save "p7_poisA.dta", replace
expand k
drop k
save "p7_poisB.dta", replace
pairrun , fa(p7_poisA.dta) ca(wooldid ypos i t ttre [pw=k], att poisson) fb(p7_poisB.dta) cb(wooldid ypos i t ttre, att poisson) tag(POIS)
pairrun , fa(p7_poisA.dta) ca(wooldid ypos i t ttre [pw=k], att att_i att_it att_itime poisson cluster(i)) fb(p7_poisB.dta) cb(wooldid ypos i t ttre, att att_i att_it att_itime poisson cluster(i)) tag(POISALL)
* aweights must be refused under poisson (sthlp:710)
use "p7_poisA.dta", clear
cap noi wooldid ypos i t ttre [aw=k], att poisson
dis "POIS_AW_RC=" _rc

*--- V: does pweight force a robust VCE? --------------------------------------
use "p2_A.dta", clear
qui wooldid y i t ttre [pw=k], att
matrix VPW = e(V)
dis "VCE_PW_type=[" e(wdid_vcetype) "] [" e(vcetype) "] [" e(vce) "]"
qui wooldid y i t ttre [aw=k], att robust
matrix VAWR = e(V)
dis "VCE_AWR_type=[" e(wdid_vcetype) "] [" e(vcetype) "] [" e(vce) "]"
qui wooldid y i t ttre [aw=k], att
matrix VAW = e(V)
mata: wd_rel("VPW","VAWR","PW_VS_AWROBUST")
mata: wd_rel("VPW","VAW","PW_VS_AWIID")

*--- S: more surface for the aw==duplication identity -------------------------
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att impvalues cluster(i)) fb(p2_B.dta) cb(wooldid y i t ttre, att impvalues cluster(i)) tag(IMPV)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att anticipation(1)) fb(p2_B.dta) cb(wooldid y i t ttre, att anticipation(1)) tag(ANTIC)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att oostreatedcontrols(1)) fb(p2_B.dta) cb(wooldid y i t ttre, att oostreatedcontrols(1)) tag(OOS)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att att_it esprelength(2) espostlength(3) esfixedbaseperiod) fb(p2_B.dta) cb(wooldid y i t ttre, att att_it esprelength(2) espostlength(3) esfixedbaseperiod) tag(ESFIX)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att att_it esprelength(2) espostlength(3) esrelativeto(-2)) fb(p2_B.dta) cb(wooldid y i t ttre, att att_it esprelength(2) espostlength(3) esrelativeto(-2)) tag(ESREL)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att timetrends) fb(p2_B.dta) cb(wooldid y i t ttre, att timetrends) tag(TT)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att unconditionalse robust) fb(p2_B.dta) cb(wooldid y i t ttre, att unconditionalse robust) tag(UNCOND)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att att_it coarsecohortcontrols(sg)) fb(p2_B.dta) cb(wooldid y i t ttre, att att_it coarsecohortcontrols(sg)) tag(CCC)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att interactivecontrols(x)) fb(p2_B.dta) cb(wooldid y i t ttre, att interactivecontrols(x)) tag(IXC)
pairrun , fa(p2_A.dta) ca(wooldid y i t ttre [aw=k], att att_it pretrendpower(0.5) esprelength(2) espostlength(2)) fb(p2_B.dta) cb(wooldid y i t ttre, att att_it pretrendpower(0.5) esprelength(2) espostlength(2)) tag(PTP)

dis "PROBE7_DONE"
