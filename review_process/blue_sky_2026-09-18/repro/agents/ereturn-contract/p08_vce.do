* ====== VCE bookkeeping, e(sample) semantics, and the cfx name surface ======
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid
global CKCOUT "vioV.txt"
cap erase "vioV.txt"
include "ckc.doi"
include "ckc2.doi"

cap prog drop VCEDUMP
prog def VCEDUMP
  args tag rc
  local hascv = 0
  local cv ""
  if "`e(clustvar)'" != "" {
    local hascv = 1
    local cv "`e(clustvar)'"
  }
  dis "VCE|`tag'|rc=`rc'|vce=`e(vce)'|vcetype=`e(vcetype)'|clustvar=`cv'|N_clust=" e(N_clust) "|Nclus=" e(Nclus) "|df_r=" e(df_r) "|wdid_vcespec=`e(wdid_vcespec)'|wt=`e(wtype)'|wexp=`e(wexp)'"
end

local ES "esfixedbaseperiod esprelength(3) espostlength(3)"
use D1.dta, clear
gen long clu2 = 1 + mod(id,7)

* ---- the VCE bookkeeping grid ----
cap noi wooldid y id t ttre, att `ES'
VCEDUMP plain `=_rc'
cap noi wooldid y id t ttre, att `ES' robust
VCEDUMP robust `=_rc'
cap noi wooldid y id t ttre, att `ES' cluster(clu)
VCEDUMP cluster1 `=_rc'
cap noi wooldid y id t ttre, att `ES' cluster(clu clu2)
VCEDUMP cluster2way `=_rc'
cap noi wooldid y id t ttre, att `ES' customvce(cluster clu)
VCEDUMP cvce_cluster `=_rc'
cap noi wooldid y id t ttre, att `ES' customvce(robust)
VCEDUMP cvce_robust `=_rc'
cap noi wooldid y id t ttre, att `ES' customvce(hc3)
VCEDUMP cvce_hc3 `=_rc'
cap noi wooldid y id t ttre, att `ES' unconditionalse
VCEDUMP uncond `=_rc'
cap noi wooldid y id t ttre, att `ES' unconditionalse cluster(clu)
VCEDUMP uncondclu `=_rc'
cap noi wooldid y id t ttre, att `ES' [pw=pwt]
VCEDUMP pweight `=_rc'
cap noi wooldid y id t ttre, att `ES' [aw=awt]
VCEDUMP aweight `=_rc'
cap noi wooldid y id t ttre, att `ES' impvalues cluster(clu)
VCEDUMP impvalues `=_rc'
cap noi wooldid ypos id t ttre, att `ES' poisson
VCEDUMP poisson `=_rc'
cap noi wooldid ypos id t ttre, att `ES' poisson cluster(clu)
VCEDUMP poissonclu `=_rc'
cap noi wooldid y id t ttre, att `ES' fastnose
VCEDUMP fastnose `=_rc'
cap noi wooldid y id t ttre, att `ES' cluster(clu) bootstrap(standard) reps(20) rseed(3)
VCEDUMP boot `=_rc'
cap noi wooldid y id t ttre, att `ES' cluster(clu) permute(shuffle) permreps(20) rseed(3)
VCEDUMP perm `=_rc'

* ---- e(sample) vs the four sample concepts ----
dis "########## e(sample) SEMANTICS ##########"
use D1.dta, clear
qui count
dis "RAW_N=" r(N)
replace y = . if mod(_n,37)==0
qui count if !missing(y)
dis "NONMISS_N=" r(N)
gen byte keepme = (t<=10)
qui count if keepme==1
dis "IFROWS_N=" r(N)
cap noi wooldid y id t ttre if keepme==1, att att_i att_it att_itime `ES'
dis "SAMP_RC=" _rc
qui count if e(sample)
dis "ESAMPLE_N=" r(N) " e(N)=" e(N) " e(N_dropped)=" e(N_dropped)
qui count if e(sample) & keepme!=1
dis "ESAMPLE_OUTSIDE_IF=" r(N)
qui count if e(sample) & missing(y)
dis "ESAMPLE_WITH_MISSING_Y=" r(N)
ckc SAMPIF

* an if that is missing (not 0) on some rows - the BS-37 "exactly 1" rule
use D1.dta, clear
gen byte flag = .
replace flag = 1 if t<=10
replace flag = 0 if t==11
* t==12 stays missing
cap noi wooldid y id t ttre if flag==1, att `ES'
dis "FLAGIF_RC=" _rc " e(N)=" e(N) " e(N_dropped)=" e(N_dropped)
qui count if e(sample)
dis "FLAGIF_ESAMPLE=" r(N)

* ---- the cfx name surface (22-character silent drop) ----
dis "########## CFX NAMES ##########"
use D3.dta, clear
cap noi wooldid y cty t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(basic) makecfxplots
dis "CFX_RC=" _rc
foreach m in `: e(matrices)' {
  if strpos("`m'","cfx") dis "CFXMAT|`m'|" rowsof(e(`m')) "x" colsof(e(`m')) "|len=" strlen("`m'")
}
dis "CFXMATS=`e(wdid_cfxmats)'"
ckc CFXBASIC
* a long continuous-treatment variable name
use D3.dta, clear
rename dose treatment_intensity_measure_long
cap noi wooldid y cty t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) contreat(treatment_intensity_measure_long) contreatcontrols(basic) makecfxplots
dis "CFXLONG_RC=" _rc
foreach m in `: e(matrices)' {
  if strpos("`m'","cfx") dis "CFXMATL|`m'|len=" strlen("`m'")
}
dis "CFXMATSL=`e(wdid_cfxmats)'"

* ---- svy abuse and the properties(svyr) declaration ----
dis "########## SVY ##########"
use D1.dta, clear
gen byte one = 1
cap noi svyset id [pw=pwt], strata(one)
dis "SVYSET_RC=" _rc
cap noi svy: wooldid y id t ttre, att
dis "SVY_RC=" _rc
cap noi bootstrap, reps(3) seed(5): wooldid y id t ttre, att
dis "BSPREFIX_RC=" _rc " e(cmd)=`e(cmd)'"

dis "VCE_DONE"
