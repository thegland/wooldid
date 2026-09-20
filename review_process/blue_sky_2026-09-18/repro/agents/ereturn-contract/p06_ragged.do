* ================= RAGGED / ASYMMETRIC-SUBGROUP GEOMETRIES =================
* Hunting the F3 mechanism: a namelist grown in one loop while the matrix is grown in
* another, where ONE ROUND APPENDS AND A LATER ROUND DOES NOT.  Stata silently PADS a
* short namelist by repeating the last name, so every later row wears the next
* quantity's label.  The detector is ckc's BDUP / BDIFF / RELYR / ESBDIFF assertions.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vioR.txt"
cap erase "vioR.txt"
include "ckc.doi"
include "ckc2.doi"
cap prog drop nameslog
prog def nameslog
  args tag
  tempname fh
  file open `fh' using "names.txt", write text append
  foreach s in `: e(scalars)' {
    file write `fh' "N|`tag'|SCA|`s'" _n
  }
  foreach s in `: e(macros)' {
    file write `fh' "N|`tag'|MAC|`s'" _n
  }
  foreach s in `: e(matrices)' {
    file write `fh' "N|`tag'|MAT|`s'" _n
  }
  file close `fh'
end
cap prog drop RUN
prog def RUN
  args tag rc
  dis "==================== `tag' (rc=`rc') ===================="
  if "`rc'" != "0" {
    ckcw INF "`tag'" SKIPPED "run failed rc=`rc'; e() is STALE, assertions skipped"
    dis "SKIPPED|`tag'|rc=`rc'"
    exit
  }
  dis "BCOLS|`tag'|`: colfullnames e(b)'"
  cap noi ckc "`tag'"
  dis "CKC_RC=" _rc
  cap noi ckc2 "`tag'"
  dis "CKC2_RC=" _rc
  nameslog "`tag'"
end

* ---------- D5: one subgroup adopts almost immediately (1 pre period), the other late ----
clear
set obs 100
gen long id = _n
gen byte sg = mod(id,2)
gen double ttre = .
replace ttre = 2 if sg==0 & id<=50
replace ttre = 9 if sg==1 & id<=50
gen long clu = ceil(id/5)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .2*t + 1.1*post + rnormal()
save D5.dta, replace

cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(4) espostlength(4) subgroups(sg) cluster(clu)
local lrc = _rc
dis "R01_RC=`lrc'"
RUN R01_shortpre `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(6) espostlength(6) subgroups(sg) cluster(clu)
local lrc = _rc
dis "R02_RC=`lrc'"
RUN R02_shortpre6 `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(4) espostlength(4) subgroups(sg) cluster(clu) jointtests adversarial
local lrc = _rc
dis "R03_RC=`lrc'"
RUN R03_shortprejt `lrc'

* ---------- D6: one subgroup has NO usable pre-treatment outcome at all -------------
use D5.dta, clear
replace y = . if sg==0 & t < ttre & ttre<.
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(4) espostlength(4) subgroups(sg) cluster(clu)
local lrc = _rc
dis "R04_RC=`lrc'"
RUN R04_nopreoutcome `lrc'

* ---------- D7: a subgroup present only in one cohort -------------------------------
clear
set obs 120
gen long id = _n
gen double ttre = .
replace ttre = 4 if inrange(id,1,30)
replace ttre = 7 if inrange(id,31,60)
replace ttre = 10 if inrange(id,61,80)
gen byte sg = 0
replace sg = 1 if inrange(id,1,30)
gen long clu = ceil(id/6)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .2*t + 1.0*post + rnormal()
save D7.dta, replace
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) subgroups(sg) cluster(clu)
local lrc = _rc
dis "R05_RC=`lrc'"
RUN R05_sgonecohort `lrc'

* ---------- D8: the ORIGINAL trigger shape - a TIME-INVARIANT dose with subgroups ----
clear
set obs 36
gen long cty = _n
gen double ttre = .
replace ttre = 4 if inrange(cty,1,10)
replace ttre = 7 if inrange(cty,11,22)
gen byte sg = mod(cty,2)
gen long clu = 1 + mod(cty,9)
gen double cfe = rnormal()
gen double dose0 = 0.5 + runiform()
expand 12
bys cty: gen int t = _n
expand 6
bys cty t: gen byte pid = _n
gen byte post = (t>=ttre) & ttre<.
gen double dose = dose0
gen double y = cfe + .1*t + 0.7*post + 0.4*post*dose + rnormal()
save D8.dta, replace
cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic) subgroups(sg)
local lrc = _rc
dis "R06_RC=`lrc'"
RUN R06_tinvdose_sg `lrc'

cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(untreatedZero) subgroups(sg)
local lrc = _rc
dis "R07_RC=`lrc'"
RUN R07_tinvdose_sg_uz `lrc'

* dose varying only POST onset (dose = dose0*post), with subgroups: the recon F5 DGP
use D8.dta, clear
replace dose = dose0*post
cap noi wooldid y cty t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(untreatedZero) subgroups(sg)
local lrc = _rc
dis "R08_RC=`lrc'"
RUN R08_postdose_sg `lrc'

* ---------- D9: gapped / biennial time, off-grid onsets ------------------------------
clear
set obs 90
gen long id = _n
gen double ttre = .
replace ttre = 2006 if inrange(id,1,25)
replace ttre = 2010 if inrange(id,26,50)
replace ttre = 2014 if inrange(id,51,65)
gen byte sg = mod(id,2)
gen long clu = ceil(id/5)
gen double ufe = rnormal()
expand 9
bys id: gen int t = 2000 + 2*(_n-1)
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .05*t + .9*post + rnormal()
save D9.dta, replace
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sg) cluster(clu)
local lrc = _rc
dis "R09_RC=`lrc'"
RUN R09_biennial `lrc'

* ---------- D10: singleton treated cohorts ------------------------------------------
clear
set obs 80
gen long id = _n
gen double ttre = .
replace ttre = 4 if id==1
replace ttre = 6 if id==2
replace ttre = 8 if inrange(id,3,30)
gen byte sg = mod(id,2)
gen long clu = ceil(id/4)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .2*t + 1.0*post + rnormal()
save D10.dta, replace
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) cluster(clu)
local lrc = _rc
dis "R10_RC=`lrc'"
RUN R10_singletoncoh `lrc'

cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sg) cluster(clu)
local lrc = _rc
dis "R11_RC=`lrc'"
RUN R11_singletoncoh_sg `lrc'

dis "RAGGED_DONE"
