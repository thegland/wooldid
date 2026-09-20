* RF01_fastnose_permanalyticse.do
* SELF-CONTAINED repro for PX-01 (skeptic's independent version).
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do RF01_fastnose_permanalyticse.do
*
* CLAIM: under permute(), a user-typed `fastnose' silently overrides `permanalyticse'.
* The run still posts e(perm_analyticse)=1 and still prints a header and a footer that
* describe analytic standard errors and an analytic e(V)/e(V_es_*) that are not there.
* The twin request permtstat is REFUSED with r(198) on reasoning that applies verbatim.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 771103
set sortseed 771103

* --- a plain staggered panel: 60 firms x 12 years, onsets 2002/2004/2006, 15 never-treated
set obs 60
gen int firm = _n
gen double fe_i = rnormal(0, 1.4)
gen byte grp = mod(_n, 4)
gen int onset = .
replace onset = 2002 if grp == 1
replace onset = 2004 if grp == 2
replace onset = 2006 if grp == 3
gen double size = 3 + rnormal(0, .8)
expand 12
bysort firm: gen int year = 1997 + _n
gen double fe_t = .12*(year - 1997)
gen double eps = rnormal(0, .9)
bysort firm (year): replace eps = .5*eps[_n-1] + eps if _n > 1
gen int rel = year - onset
gen double effect = 0
replace effect = 0.35 + 0.18*rel + 0.25*(mod(firm,3)) if onset < . & rel >= 0
gen double y = 2.0 + fe_i + fe_t + 0.4*size + effect + eps
drop rel effect

cap prog drop rpt
prog def rpt
  args tag
  dis "`tag'_permase=" cond("`e(perm_analyticse)'"=="", "ABSENT", "`e(perm_analyticse)'")
  dis "`tag'_vcvstate=" "`e(wdid_vcvstate)'"
  foreach m in V V_es_att b_es_att {
    cap confirm matrix e(`m')
    dis "`tag'_has_`m'=" (_rc==0)
  }
  cap confirm matrix e(wdidmainresults)
  if _rc == 0 {
    matrix WR = e(wdidmainresults)
    dis "`tag'_se_analytic_att=" WR[1,2]
    dis "`tag'_att=" WR[1,1]
    dis "`tag'_p_perm_att=" WR[1,4]
  }
  dis "`tag'_header=" "`e(wdid_vcespec)'"
end

dis "===== A: permanalyticse alone -- the documented behaviour ====="
ereturn clear
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) permanalyticse rseed(4242)
dis "A_RC=" _rc
rpt A

dis "===== B: permanalyticse + fastnose -- ACCEPTED, but fastnose wins silently ====="
ereturn clear
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) permanalyticse fastnose rseed(4242)
dis "B_RC=" _rc
rpt B

dis "===== C: the twin request permtstat + fastnose is REFUSED ====="
ereturn clear
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) permtstat fastnose rseed(4242)
dis "C_RC=" _rc

dis "===== D: the program ALREADY decides this pair is not the fast-mode default ====="
dis "----- D1: fastnose + jointtests + permute -> expected rc 0 (fnperm stand-down)"
ereturn clear
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) fastnose jointtests rseed(4242)
dis "D1_RC=" _rc
dis "----- D2: fastnose + permanalyticse + jointtests + permute -> expected r(198)"
ereturn clear
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) fastnose permanalyticse jointtests rseed(4242)
dis "D2_RC=" _rc

dis "===== E: the instruction the footer gives cannot be followed under B ====="
ereturn clear
qui wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) permanalyticse fastnose rseed(4242)
cap noi lincom _b[main_att]
dis "E_lincom_RC=" _rc
dis "---- RF01 done ----"
