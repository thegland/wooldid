* r7_W_lean.do - lean, self-contained repro for hv-permsave-comma-02 plus two new probes.
* Which of the joint/equality replication-statistic matrices are ereturned under
* boottstat, under the fixed metric, and under permute()?  Small panel, no adversarial,
* no esfixedbaseperiod, so it finishes in a couple of minutes.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 24680
set sortseed 24680
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/hv-permsave-comma-refute/w7"
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
gen byte sg = mod(_n, 2)
gen double ue = rnormal()*0.8
expand 11
bys unit: gen int year = 1999 + _n
gen double y = 1.2 + 0.05*(year-1999) + ue + rnormal()*0.5
replace y = y + (0.4 + 0.2*sg)*(year - ttre + 1)/3 if year >= ttre & ttre < .
save "panel.dta", replace

prog def whichW
  syntax , tag(string)
  foreach m in boot_jointW boot_jointWobs boot_sgW boot_sgWobs boot_advW boot_advWobs ///
               perm_jointW perm_jointWobs perm_sgW perm_sgWobs perm_advW perm_advWobs {
    capture confirm matrix e(`m')
    dis "`tag'|`m'|" cond(_rc==0,"PRESENT","absent")
  }
  dis "`tag'|jt_family=`e(boot_jt_family)'"
  dis "`tag'|joint_horiz=" e(joint_horiz)
  dis "`tag'|joint_att_allsgssame=" e(joint_att_allsgssame)
end

* ---- A: boottstat + jointtests, no subgroups
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) jointtests esprelength(2) espostlength(2) ///
    bootstrap(standard) reps(25) boottstat bootsave(bA) replace dots(0)
dis "A_RC=" _rc
if _rc == 0 whichW, tag(A_BTS)

* ---- B: same, fixed metric
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) jointtests esprelength(2) espostlength(2) ///
    bootstrap(standard) reps(25) dots(0)
dis "B_RC=" _rc
if _rc == 0 whichW, tag(B_FIX)

* ---- C: boottstat + jointtests + subgroups
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) subgroups(sg) jointtests esprelength(2) espostlength(2) ///
    bootstrap(standard) reps(25) boottstat dots(0)
dis "C_RC=" _rc
if _rc == 0 whichW, tag(C_BTS_SG)

* ---- D: same, fixed metric
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) subgroups(sg) jointtests esprelength(2) espostlength(2) ///
    bootstrap(standard) reps(25) dots(0)
dis "D_RC=" _rc
if _rc == 0 whichW, tag(D_FIX_SG)

* ---- E: permute() twin
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) subgroups(sg) jointtests esprelength(2) espostlength(2) ///
    permute(shuffle) permreps(25) dots(0)
dis "E_RC=" _rc
if _rc == 0 whichW, tag(E_PERM_SG)

* ---- NEW: is the observed joint W recoverable from the bootsave() file under boottstat?
capture noisily use "bA.dta", clear
dis "F_useRC=" _rc
if _rc == 0 {
  qui ds
  dis "F_COLS=" r(varlist)
  capture confirm variable jointW_att
  dis "F_jointW_att_RC=" _rc
  if _rc == 0 {
    qui su jointW_att if _rep == 0, meanonly
    dis "F_jointW_att_row0=" %21x r(mean) "  n=" r(N)
  }
}

* ---- NEW S3: permute() and bootstrap() given together
use "panel.dta", clear
capture noisily wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) ///
    bootstrap(standard) reps(20) dots(0)
dis "S3_RC=" _rc
dis "S3_perm_reps=" e(perm_reps) " S3_boot_reps=" e(boot_reps)
dis "S3_perm_method=`e(perm_method)' S3_boot_method=`e(boot_method)'"

dis "ALLDONE"
