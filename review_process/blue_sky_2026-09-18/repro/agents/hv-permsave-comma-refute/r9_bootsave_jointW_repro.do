* r9_bootsave_jointW_repro.do  -- SELF-CONTAINED repro, NEW finding (not seeded).
*
* sthlp:557 (the bootsave() entry):  "Under boottstat the replication joint-test
*   statistics are added as jointW_name and the adversarial statistic as advW, with the
*   observed statistics on the _rep = 0 row."
* ado:14905 (wooldid_bootsavemk):    "...the per-replication joint and adversarial
*   statistics under boottstat, with the observed statistic in row 0"
*
* ACTUAL: row _rep = 0 of jointW_<family> holds element k of the OBSERVED EVENT-STUDY
* COEFFICIENT VECTOR (wdbjobs), not the observed joint Wald statistic. wooldid_bootsavemk
* (ado:14919-14925) does
*        JO = st_matrix(jobs)          // jobs = wdbjobs = the observed b, one column
*                                      //        per pre-treatment horizon per family
*        st_store(., idx, (cols(JO) >= k ? JO[1, k] : .) \ JW[., k])
* pairing the k-th COLUMN OF e(boot_jointW) (one per test FAMILY) with the k-th column of
* a COEFFICIENT vector. With one family requested, k = 1, so row 0 gets the first
* pre-treatment horizon's observed coefficient.
*
* Run:  stata-mp -e do r9_bootsave_jointW_repro.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 13579
set sortseed 13579
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma-refute/w9"
which wooldid
cap shell rm -f *.dta

* ---------------- an ordinary staggered panel: 60 units, 11 years, 3 onset cohorts
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

local jh   = e(joint_horiz)
local nrep = e(boot_completed)
dis "REPORTED_joint_horiz = " `jh'   "   (" %21x `jh' ")"
dis "BOOT_JT_FAMILY = `e(boot_jt_family)'    BOOT_COMPLETED = " `nrep'
tempname JW BES
matrix `JW' = e(boot_jointW)
local jwc : colnames `JW'
dis "e(boot_jointW): " rowsof(`JW') " x " colsof(`JW') "  colnames = `jwc'"
matrix `BES' = e(b_es_att)
local bc : colnames `BES'
dis "e(b_es_att) colnames = `bc'"
matrix list `BES', format(%21x)
local b_pre2 = `BES'[1, 1]

use "bj.dta", clear
qui ds
dis "FILE_COLS = " r(varlist)
dis "FILE_N = " _N

* ---- (a) the value the help calls "the observed statistic"
qui su jointW_att if _rep == 0, meanonly
local w0 = r(mean)
dis "jointW_att[_rep==0] = " `w0' "   (" %21x `w0' ")"
dis "IS_NEGATIVE = " (`w0' < 0) "   <-- a Wald statistic cannot be negative"

qui su jointW_att if _rep > 0
dis "jointW_att replications: n = " r(N) "  min = " r(min) "  max = " r(max)
dis "   (all replication values are >= 0, as Wald statistics must be)"

* ---- (b) it is the first pre-treatment observed COEFFICIENT
dis "e(b_es_att)[pre2]   = " `b_pre2' "   (" %21x `b_pre2' ")"
dis "ROW0_EQUALS_b_pre2 = " (abs(`w0' - `b_pre2') <= 1e-13*max(abs(`b_pre2'), 1e-12))
qui su ES_att_pre2 if _rep == 0, meanonly
dis "file ES_att_pre2[_rep==0] = " r(mean) "   (" %21x r(mean) ")"

* ---- (c) the p the help says is reconstructable from the file
qui count if _rep > 0 & jointW_att < . & jointW_att >= `w0'
local ge = r(N)
qui count if _rep > 0 & jointW_att < .
local n = r(N)
local precon = (1 + `ge')/(`n' + 1)
dis "RECONSTRUCTED p from the file = " `precon' "   (" `ge' " of " `n' " replications at or above row 0)"
dis "REPORTED p                    = " `jh'
dis "AGREE = " (abs(`precon' - `jh') < 1e-12)

* ---- the adversarial twin, same mechanism (advW <- wdbaobs)
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) esprelength(2) espostlength(2) ///
    esfixedbaseperiod adversarial bootstrap(standard) reps(25) boottstat bootsave(ba) replace dots(0)
dis "ADV_RC = " _rc
dis "ADV_pval_att = " e(adversarial_pval_att)
capture confirm file "ba.dta"
if _rc == 0 {
  use "ba.dta", clear
  qui ds
  dis "ADVFILE_COLS = " r(varlist)
  capture confirm variable advW
  dis "advW_present_rc = " _rc
  if _rc == 0 {
    qui su advW if _rep == 0, meanonly
    dis "advW[_rep==0] = " r(mean) "   (" %21x r(mean) ")"
    qui su advW if _rep > 0
    dis "advW replications: n = " r(N) " min = " r(min) " max = " r(max)
  }
}
dis "ALLDONE"
