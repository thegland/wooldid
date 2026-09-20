* r5_floor.do -- SELF-CONTAINED: what is the smallest p-value a permexact run can print?
* The help's p-value MAP (wooldid.sthlp line 607, the row "permute() + permexact") says
*   "Plain shares over the M enumerated assignments ... Floor 1/M."
* The help's own floors paragraph (sthlp:512) says instead
*   "under permexact with M assignments it is 2/M under the default and 1/M under the
*    absolute-value conventions".
* This probe measures it. Design: 6 clusters, onsets 5,5,7,.,.,. -> M = 6!/(2!1!3!) = 60.
* The data give a large, clean treatment effect, so the observed statistic is the extreme
* of its own enumeration and the reported p sits AT the floor.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 6
gen long unit = _n
gen double onset = .
replace onset = 5 in 1
replace onset = 5 in 2
replace onset = 7 in 3
gen double ue = rnormal()*0.2
expand 9
bysort unit: gen int year = _n
gen double te = rnormal()*0.2
gen byte post = (year >= onset) & onset < .
gen double y = 1 + 0.1*year + ue + te + 3.0*post
tempfile d
save `d', replace

prog def shownt
  syntax , tag(string) m(real)
  foreach mm in perm_p perm_p_lo perm_p_hi {
    cap confirm matrix e(`mm')
    if _rc != 0 continue
    tempname X
    matrix `X' = e(`mm')
    local cn : colnames `X'
    local j 0
    foreach c of local cn {
      local ++j
      local v = `X'[1, `j']
      dis "`tag' `mm' `c' = " %20.17f `v' "   p*M = " %12.8f `v'*`m' "   (M=`m')"
    }
  }
  cap confirm matrix e(perm_n)
  if _rc==0 matrix list e(perm_n)
end

dis "===== A: DEFAULT convention + permexact ====="
capture noisily wooldid y unit year onset, att esprelength(2) espostlength(2) cluster(unit) ///
  permute(shuffle) permexact jointtests dots(0)
dis "A_RC=" _rc
dis "A_reps=" e(perm_reps)
shownt , tag(A) m(60)
dis "A_joint_pretrend_p=" e(joint_pretrend_p)
foreach s in joint_pretrend_p joint_horiz joint_horiz_p pretrendjoint_p {
  cap dis "A_`s'=" e(`s')
}

dis "===== B: permabsp + permexact ====="
use `d', clear
capture noisily wooldid y unit year onset, att esprelength(2) espostlength(2) cluster(unit) ///
  permute(shuffle) permexact permabsp jointtests dots(0)
dis "B_RC=" _rc
shownt , tag(B) m(60)

dis "===== C: permside(right) + permexact ====="
use `d', clear
capture noisily wooldid y unit year onset, att esprelength(2) espostlength(2) cluster(unit) ///
  permute(shuffle) permexact permside(right) dots(0)
dis "C_RC=" _rc
shownt , tag(C) m(60)

dis "===== D: permrawp + permexact ====="
use `d', clear
capture noisily wooldid y unit year onset, att esprelength(2) espostlength(2) cluster(unit) ///
  permute(shuffle) permexact permrawp dots(0)
dis "D_RC=" _rc
shownt , tag(D) m(60)
dis "DONE_R5"
