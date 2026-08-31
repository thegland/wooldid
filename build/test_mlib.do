* build/test_mlib.do - the mata library (UNWIRED: not under test/, so the cert component
* does not glob it, because no library is shipped yet and the design is undecided).
* Run from test/ once a library exists: stata-mp -e do ../build/test_mlib.do
* the shipped mata library: discovery, clear-mata survival, and the
* cleanmatrices cure (SI-18). Run from test/: stata-mp -e do test_mlib.do
*
* M1 the library is discoverable and reports the expected function set
* M2 functions resolve after `clear mata`, repeatably, with no inline definitions
* M3-M6 the four previously-failing cleanmatrices cases now complete
* M7 a normal run with no mata-dependent options is unaffected
* M8 BISECTABILITY: with the library hidden, M3 fails again - proving these checks
*    test the library rather than passing vacuously
* Prints MLIB_ALL_PASS iff all pass.
clear all
local REPO ..
local fails 0

* ---- M1: discovery + contents
adopath ++ "`REPO'/wooldidinstall"
mata: mata mlib index
if strpos("`c(matalibs)'", "lwooldid") == 0 {
  di as error "M1 FAIL: lwooldid not in c(matalibs)"
  local ++fails
}
qui log using mlib_contents.txt, replace text
cap mata: mata describe using lwooldid
qui log close
local mlibn = 0
foreach f in wooldid_gtx wooldid_advx wooldid_permx wooldid_permp wooldid_permxp ///
             wooldid_bootdw wooldid_bootgrab wooldid_bootfin wooldid_bootves ///
             wooldid_bootveskey wooldid_bootdist wooldid_bootjw wooldid_bootpw wooldid_bootci {
  if strpos(fileread("mlib_contents.txt"), "`f'(") > 0 local ++mlibn
}
di "M1_functions_found=" `mlibn' " of 14 spot-checked"
if `mlibn' < 14 {
  di as error "M1 FAIL: library is missing expected functions"
  local ++fails
}
erase mlib_contents.txt

* ---- M2: survives clear mata with NO inline definitions in the session
clear mata
cap mata: st_numscalar("m2a", wooldid_bootveskey("att", "pre1") == "ES_att_pre1")
if _rc != 0 | m2a != 1 {
  di as error "M2a FAIL: library function did not resolve after clear mata (rc " _rc ")"
  local ++fails
}
clear mata
cap mata: st_numscalar("m2b", wooldid_bootveskey("ame_it_s3", "post0") == "ES_ame_it_SG3_post0")
if _rc != 0 | m2b != 1 {
  di as error "M2b FAIL: second clear mata (rc " _rc ")"
  local ++fails
}
di "M2_survives_clear_mata=1"

* ---- the cure cases. cleanmatrices wipes mata mid-run; the library reloads on demand.
do "`REPO'/wooldidinstall/wooldid.ado"
use "`REPO'/test/data/dta_binary.dta", clear
tempfile base
save `base'

* M3: gteffects calls mata AFTER the clear, in the same run
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) gteffects cleanmatrices
di "M3_gteffects_cleanmatrices_rc=" _rc
if _rc != 0 {
  di as error "M3 FAIL"
  local ++fails
}

* M4: the adversarial line test likewise
use `base', clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod esrelativeto(-2) adversarial cleanmatrices
di "M4_adversarial_cleanmatrices_rc=" _rc
if _rc != 0 {
  di as error "M4 FAIL"
  local ++fails
}

* M5: session poisoning - a bootstrap run AFTER an earlier cleanmatrices run
use `base', clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) cleanmatrices
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) bootstrap(bayes) reps(15) rseed(3)
di "M5_bootstrap_after_cleanmatrices_rc=" _rc
if _rc != 0 {
  di as error "M5 FAIL"
  local ++fails
}

* M6: and permute after cleanmatrices
use `base', clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) cleanmatrices
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20)
di "M6_permute_after_cleanmatrices_rc=" _rc
if _rc != 0 {
  di as error "M6 FAIL"
  local ++fails
}

* ---- M7: a normal run is untouched
use `base', clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
di "M7_plain_rc=" _rc
if _rc != 0 {
  di as error "M7 FAIL"
  local ++fails
}

* ---- M8: bisectability. Hide the library from the adopath and confirm M3 fails again,
* so these checks cannot pass vacuously if the library ever stops being shipped.
adopath - "`REPO'/wooldidinstall"
mata: mata mlib index
clear mata
do "`REPO'/wooldidinstall/wooldid.ado"
use `base', clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) gteffects cleanmatrices
di "M8_without_library_rc=" _rc " (expect 3499)"
if _rc == 0 {
  di as error "M8 FAIL: the cure cases pass even without the library - they prove nothing"
  local ++fails
}

if `fails' == 0 di as result "MLIB_ALL_PASS"
else di as error "MLIB_FAILS=`fails'"
