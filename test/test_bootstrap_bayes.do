* test_bootstrap_bayes.do - commit-C engine validation: bayes bootstrap end-to-end.
* Run from test/: stata-mp -e do test_bootstrap_bayes.do
* O1/O2: BITWISE oracle vs a manual rwgen-bayes-driven rep loop (both alpha
* conventions - proves draws, weight injection and accumulation exactly reproduce
* the rwgen/bayesboot weight algorithm). O3: point estimates never change. O4 aw;
* O5 poisson pw; O6 subgroups (SG tables rewired, pwsgs blanked); O7 esf + the A8
* raw-esprelength emission; O8 replay; O9 state hygiene; O10 d3-class contreat
* parity with the analytic scrub. Prints BOOTBAYES_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0

* panel builders
cap prog drop mkpan
prog def mkpan
  clear
  set obs 12
  gen unit = _n
  gen sg = mod(unit, 2)
  gen ttre = cond(unit<=4, 4, cond(unit<=8, 6, .))
  expand 8
  bysort unit: gen t = _n
  set seed 11
  gen y = 2*(t>=ttre & ttre<.) + 0.3*(t>=ttre & ttre<.)*sg + 0.1*unit + 0.05*t + rnormal()
  gen w0 = 0.5 + mod(unit, 3)
  gen ycnt = max(0, round(exp(0.15*(t>=ttre & ttre<.) + 0.02*t + 0.05*mod(unit,4))))
  replace ycnt = ycnt + 1 + mod(unit+t, 3)
end

* ---- O1: ORACLE - engine vs manual rwgen loop, same seed, unweighted, size alpha
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(12) rseed(424)
matrix ENG = e(bootdist)
local engK = colsof(ENG)
mkpan
rwgen bayes RW, reps(12) id(unit) rseed(424)
matrix MAN = J(12, `engK', .)
forval r = 1/12 {
  qui wooldid y unit t ttre [aw=RW`r'], att cluster(unit) esprelength(1) espostlength(1) fastnose
  matrix MAN[`r', 1] = e(b)
}
local w1 0
forval r = 1/12 {
  forval c = 1/`engK' {
    local w1 = max(`w1', reldif(ENG[`r',`c'], MAN[`r',`c']))
  }
}
di "O1_ORACLE_SIZE_MAXRELDIF=" `w1'
if `w1' > 1e-12 {
  di as error "O1 FAIL"
  local ++fails
}

* ---- O2: equal alpha vs rwgen priorpowers(-n_c)
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(8) rseed(555) bootalpha(equal)
matrix ENG2 = e(bootdist)
mkpan
bysort unit: gen double nc = _N
gen double pv = -nc
rwgen bayes RE, reps(8) id(unit) priorpowers(pv) rseed(555)
matrix MAN2 = J(8, `engK', .)
forval r = 1/8 {
  qui wooldid y unit t ttre [aw=RE`r'], att cluster(unit) esprelength(1) espostlength(1) fastnose
  matrix MAN2[`r', 1] = e(b)
}
local w2 0
forval r = 1/8 {
  forval c = 1/`engK' {
    local w2 = max(`w2', reldif(ENG2[`r',`c'], MAN2[`r',`c']))
  }
}
di "O2_ORACLE_EQUAL_MAXRELDIF=" `w2'
if `w2' > 1e-12 {
  di as error "O2 FAIL"
  local ++fails
}

* ---- O3: point estimates NEVER change under bootstrap (replace-mode invariant)
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2)
matrix OB = e(b)
matrix OM = e(wdidmainresults)
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2) bootstrap(bayes) reps(10) rseed(2)
matrix BB = e(b)
matrix BM = e(wdidmainresults)
local o3b = mreldif(OB, BB)
local o3m 0
forval rr = 1/`=rowsof(OM)' {
  local o3m = max(`o3m', reldif(OM[`rr',1], BM[`rr',1]))
}
di "O3_POINTS_B=" `o3b' " M=" `o3m'
if `o3b' > 0 | `o3m' > 0 {
  di as error "O3 FAIL: bootstrap moved a point estimate"
  local ++fails
}

* ---- O4: weighted aw run
mkpan
cap qui wooldid y unit t ttre [aw=w0], att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(15) rseed(9)
if _rc != 0 | e(boot_completed) != 15 {
  di as error "O4 FAIL rc " _rc
  local ++fails
}
else di "O4_AW_OK se=" el(e(wdidmainresults),1,2)

* ---- O5: poisson pw run
mkpan
cap qui wooldid ycnt unit t ttre [pw=w0], att cluster(unit) esprelength(1) espostlength(1) poisson bootstrap(bayes) reps(10) rseed(31)
di "O5_POIS_RC=" _rc " completed=" e(boot_completed)
if _rc != 0 {
  di as error "O5 FAIL"
  local ++fails
}

* ---- O6: subgroups - SG tables rewired, pwsgs blanked
mkpan
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) subgroups(sg) bootstrap(bayes) reps(15) rseed(77)
if _rc != 0 {
  di as error "O6a FAIL rc " _rc
  local ++fails
}
else {
  matrix S0 = e(wdidsgresults_0)
  local sgse = S0[1, 2]
  if `sgse' >= . | `sgse' <= 0 {
    di as error "O6b FAIL: sg se not rewired"
    local ++fails
  }
  matrix PW = e(pwsgs)
  if PW[1, 2] < . {
    di as error "O6c FAIL: pwsgs not blanked"
    local ++fails
  }
  di "O6_SG_OK sgse=" `sgse'
}

* ---- O7: esf + A8 raw-esprelength emission (deep pre-history)
clear
set obs 10
gen unit = _n
gen ttre = cond(unit<=5, 8, .)
expand 12
bysort unit: gen t = _n
set seed 3
gen y = 1.5*(t>=ttre & ttre<.) + 0.1*unit + rnormal()
qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(2) esfixedbaseperiod esrelativeto(-2)
local nb0 = rowsof(e(wdidmainresults))
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(2) esfixedbaseperiod esrelativeto(-2) bootstrap(bayes) reps(10) rseed(5)
di "O7_ESF_RC=" _rc " rows " rowsof(e(wdidmainresults)) " vs " `nb0'
if _rc != 0 | rowsof(e(wdidmainresults)) != `nb0' {
  di as error "O7 FAIL"
  local ++fails
}

* ---- O8: replay after bootstrap
cap qui wooldid
di "O8_REPLAY_RC=" _rc
if _rc != 0 {
  di as error "O8 FAIL"
  local ++fails
}

* ---- O9: state hygiene
if "${WOOLDID_LASTRC}" != "" {
  di as error "O9 FAIL: LASTRC set"
  local ++fails
}

* ---- O10: d3-class contreat - bootstrap must match the analytic presentation
* (points at the structural zero, inference missing), never fabricate confidence
clear
set obs 20
gen unit = _n
gen ttre = cond(unit<=10, 5, .)
gen double dose = 0.5 + 0.1*unit
expand 10
bysort unit: gen t = _n
set seed 6
gen y = 0.4*dose*(t>=ttre & ttre<.) + 0.05*unit + 0.02*t + rnormal()
tempfile d3
save `d3'
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) contreat(dose) contreatcontrols(basic) safety(off)
matrix MC = e(wdidmainresults)
use `d3', clear
cap noi wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) contreat(dose) contreatcontrols(basic) safety(off) bootstrap(bayes) reps(12) rseed(44)
local o10rc = _rc
di "O10_AMESCRUB_RC=" `o10rc'
if `o10rc' != 0 {
  di as error "O10 FAIL rc"
  local ++fails
}
else {
  matrix MB = e(wdidmainresults)
  local m10rows : rowfullnames MC
  local m10k 0
  local o10bad 0
  foreach rn of local m10rows {
    local ++m10k
    if reldif(MB[`m10k',1], MC[`m10k',1]) > 1e-12 local o10bad 1
    if strpos("`rn'", "ame") {
      * analytic inference is scrubbed on these rows; bootstrap must not fabricate it
      if MC[`m10k',2] >= . & MB[`m10k',2] < . local o10bad 1
    }
  }
  di "O10_PARITY_BAD=" `o10bad'
  if `o10bad' {
    di as error "O10 FAIL: bootstrap diverged from the analytic scrub presentation"
    local ++fails
  }
}

if `fails' == 0 di as result "BOOTBAYES_ALL_PASS"
else di as error "BOOTBAYES_FAILS=`fails'"
