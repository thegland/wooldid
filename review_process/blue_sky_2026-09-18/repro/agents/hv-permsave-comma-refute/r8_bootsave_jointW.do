* r8_bootsave_jointW.do - NEW finding (mine, not seeded), SELF-CONTAINED.
*
* CLAIM UNDER TEST: sthlp:557 (the bootsave() option entry) says
*   "Under boottstat the replication joint-test statistics are added as jointW_name and
*    the adversarial statistic as advW, with the observed statistics on the _rep = 0 row."
* and the writer's own comment (ado:14905) says "with the observed statistic in row 0".
* But wooldid_bootsavemk writes  JO[1, k]  into row 0 of jointW_<k>, where JO = `jobs'
* = wdbjobs is the OBSERVED EVENT-STUDY COEFFICIENT VECTOR that wooldid_bootjwobs()
* consumes to COMPUTE W_obs - not a vector of observed W statistics. k indexes the
* columns of e(boot_jointW) (one per test family), so element k of a coefficient vector
* is written where the k-th family's observed Wald statistic should be.
*
* Prediction if that reading is right:
*   (a) row 0 of jointW_att can be NEGATIVE (a Wald statistic cannot),
*   (b) it equals one of the observed event-study coefficients in the same file, and
*   (c) the add-one rank p recomputed from the file does NOT equal the reported
*       e(joint_horiz).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 13579
set sortseed 13579
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma-refute/w8"
which wooldid
cap shell rm -f *.dta

clear
set obs 60
gen long unit = _n
gen byte grp = mod(_n, 5)
gen ttre = .
replace ttre = 2004 if grp == 1
replace ttre = 2006 if grp == 2
replace ttre = 2008 if grp == 3
gen double ue = rnormal()*0.8
expand 11
bys unit: gen int year = 1999 + _n
gen double y = 1.2 + 0.05*(year-1999) + ue + rnormal()*0.5
replace y = y + 0.4*(year - ttre + 1)/3 if year >= ttre & ttre < .
save "panel.dta", replace

wooldid y unit year ttre, att cluster(unit) jointtests esprelength(2) espostlength(2) ///
    bootstrap(standard) reps(25) boottstat bootsave(bj) replace dots(0)

local jh = e(joint_horiz)
local nrep = e(boot_completed)
dis "REPORTED_joint_horiz=" %21x `jh' "  = " `jh'
dis "BOOT_COMPLETED=" `nrep'
dis "BOOT_JT_FAMILY=`e(boot_jt_family)'"
tempname JW
matrix `JW' = e(boot_jointW)
dis "e(boot_jointW) dim: " rowsof(`JW') " x " colsof(`JW')
local jwc : colnames `JW'
dis "e(boot_jointW) colnames=`jwc'"
* the observed ES coefficients, for comparison
tempname BES
capture matrix `BES' = e(b_es_att)
if _rc == 0 {
  local bc : colnames `BES'
  dis "e(b_es_att) colnames=`bc'"
  matrix list `BES', format(%21x)
}

use "bj.dta", clear
qui ds
dis "FILE_COLS=" r(varlist)
dis "FILE_N=" _N

* ---- (a) row 0 of jointW_att
qui su jointW_att if _rep == 0, meanonly
local w0 = r(mean)
dis "JOINTW_ATT_ROW0=" %21x `w0' " = " `w0'
dis "IS_NEGATIVE=" (`w0' < 0)

* ---- replication values, for sanity (these should all be >= 0)
qui su jointW_att if _rep > 0
dis "JOINTW_ATT_REPS n=" r(N) " min=" r(min) " max=" r(max) " mean=" r(mean)

* ---- (b) does row 0 equal one of the observed coefficient columns?
qui ds _rep jointW_* se_* advW* w_*, not
local coefs `r(varlist)'
foreach v of local coefs {
  qui su `v' if _rep == 0, meanonly
  local c0 = r(mean)
  if abs(`c0' - `w0') <= 1e-14*max(abs(`w0'),1e-12) dis "ROW0_MATCHES_COEF `v' = " %21x `c0'
  else dis "  (no match) `v' row0 = " %21x `c0'
}

* ---- (c) is the reported p recoverable from the pair, as the help says?
qui count if _rep > 0 & jointW_att < . & jointW_att >= `w0'
local ge = r(N)
qui count if _rep > 0 & jointW_att < .
local n = r(N)
local precon = (1 + `ge')/(`n' + 1)
dis "RECON_p_from_file=" `precon' "   (ge=`ge' of n=`n')"
dis "REPORTED_p=" `jh'
dis "MATCH=" (abs(`precon' - `jh') < 1e-12)

* what W_obs WOULD reproduce the reported p?
local kneed = round(`jh'*(`n' + 1) - 1)
dis "TO_GET_REPORTED_p_NEED ge=" `kneed' " of n=" `n'
sort jointW_att
list _rep jointW_att if _rep > 0, clean noobs

dis "ALLDONE"
