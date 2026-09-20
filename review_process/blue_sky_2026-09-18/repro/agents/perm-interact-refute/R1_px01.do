* R1_px01.do -- SELF-CONTAINED. Skeptic's independent test of PX-01 and PX-02.
* PX-01: does fastnose silently override permanalyticse under permute()?
* PX-02: is permute()+customvce() reachable?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 771103
set sortseed 771103

*------------------------------------------------------------------ my own DGP
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

*--------------------------------------------------- reporting helper
cap prog drop rpt
prog def rpt
  args tag
  dis "== `tag' : rc=" _rc
  dis "`tag'_ecmd=" "`e(cmd)'"
  dis "`tag'_permase=" cond("`e(perm_analyticse)'"=="", "ABSENT", "`e(perm_analyticse)'")
  dis "`tag'_vcvstate=" "`e(wdid_vcvstate)'"
  cap confirm matrix e(V)
  dis "`tag'_hasV=" (_rc==0)
  cap confirm matrix e(V_es_att)
  dis "`tag'_hasVesatt=" (_rc==0)
  cap confirm matrix e(b_es_att)
  dis "`tag'_hasbesatt=" (_rc==0)
  cap confirm matrix e(b)
  dis "`tag'_hasb=" (_rc==0)
  cap confirm matrix e(wdidmainresults)
  if _rc == 0 {
    matrix WR = e(wdidmainresults)
    dis "`tag'_colnames=" "`: colnames WR'"
    matrix list WR
  }
  dis "`tag'_vcespec=" "`e(wdid_vcespec)'"
end

*=========================== A: permanalyticse ALONE (the promised behaviour)
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) permanalyticse rseed(4242)
dis "A_RC=" _rc
rpt A

*=========================== B: fastnose + permanalyticse (the accused pair)
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) permanalyticse fastnose rseed(4242)
dis "B_RC=" _rc
rpt B

*=========================== C: fast-mode DEFAULT (no permanalyticse, no fastnose)
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) rseed(4242)
dis "C_RC=" _rc
rpt C

*=========================== D: fastnose + permtstat  (documented refusal)
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) permtstat fastnose rseed(4242)
dis "D_RC=" _rc

*=========================== E: fastnose + jointtests + permute  (claimed ACCEPTED)
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) fastnose jointtests rseed(4242)
dis "E_RC=" _rc

*=========================== F: fastnose + permanalyticse + jointtests (claimed REFUSED)
cap noi wooldid y firm year onset, att cluster(firm) esprelength(2) espostlength(2) ///
   permute(shuffle) permreps(30) fastnose permanalyticse jointtests rseed(4242)
dis "F_RC=" _rc

*=========================== G: does the B run's p-value equal A's or C's?
* (checks whether the silent override moves any NUMBER)
dis "---- done PX01 block ----"

*=========================== PX-02 reachability attempts
cap noi wooldid y firm year onset, att cluster(firm) permute(shuffle) permreps(30) customvce(cluster firm)
dis "H_RC=" _rc
cap noi wooldid y firm year onset, att permute(shuffle) permreps(30) customvce(cluster firm)
dis "I_RC=" _rc
cap noi wooldid y firm year onset, att permute(shuffle) permreps(30) customvce(robust) cluster(firm) safety(off)
dis "J_RC=" _rc
* empty customvce()
cap noi wooldid y firm year onset, att cluster(firm) permute(shuffle) permreps(30) customvce()
dis "K_RC=" _rc
* customvce with a single space
cap noi wooldid y firm year onset, att cluster(firm) permute(shuffle) permreps(30) customvce( )
dis "L_RC=" _rc
dis "---- done PX02 block ----"
