*==========================================================================
* probeH.do -- the metamorphic relation an applied user actually assumes:
*   "rows the estimator drops anyway cannot change the answer".
* Variants: a cluster whose outcome is missing on every row (present vs
* physically dropped); a cluster with a MISSING CLUSTER ID (present vs
* dropped); partial missingness inside a cluster (control); and
* never-treated coded as an out-of-window future onset vs plain missing.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

include "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-invariance/_common.do"

mkpanel
save panelH.dta, replace
qui tab ttre, missing

local CH att att_i cluster(state) permute(shuffle) permreps(30) rseed(20260917) esprelength(2) espostlength(2) dots(0) replace

*--- H0 reference
cap noi wooldid y unit t ttre , `CH' permsave(ps_H0)
dis "RC_H0=" _rc
cap dmp H0
cap asgdmp H0

*--- H1: state 5 (a 2006 adopter) has y missing on EVERY row, rows kept
use panelH.dta, clear
qui replace y = . if state == 5
cap noi wooldid y unit t ttre , `CH' permsave(ps_H1ymiss)
dis "RC_H1ymiss=" _rc
cap dmp H1ymiss
cap asgdmp H1ymiss

*--- H2: the same rows physically dropped
use panelH.dta, clear
qui drop if state == 5
cap noi wooldid y unit t ttre , `CH' permsave(ps_H2ydrop)
dis "RC_H2ydrop=" _rc
cap dmp H2ydrop
cap asgdmp H2ydrop

*--- H3: state 5 has a MISSING CLUSTER ID, rows kept
use panelH.dta, clear
qui replace state = . if state == 5
cap noi wooldid y unit t ttre , `CH' permsave(ps_H3cmiss)
dis "RC_H3cmiss=" _rc
cap dmp H3cmiss
cap asgdmp H3cmiss

*--- H4 control: y missing on 3 periods of state 5 only (cluster survives)
use panelH.dta, clear
qui replace y = . if state == 5 & t <= 2002
cap noi wooldid y unit t ttre , `CH' permsave(ps_H4partial)
dis "RC_H4partial=" _rc
cap dmp H4partial
cap asgdmp H4partial

*--- H4b the same rows physically dropped
use panelH.dta, clear
qui drop if state == 5 & t <= 2002
cap noi wooldid y unit t ttre , `CH' permsave(ps_H4bdrop)
dis "RC_H4bdrop=" _rc
cap dmp H4bdrop
cap asgdmp H4bdrop

*--- H5: never-treated coded as an out-of-window future onset (2020)
use panelH.dta, clear
qui replace ttre = 2020 if missing(ttre)
cap noi wooldid y unit t ttre , `CH' permsave(ps_H5future)
dis "RC_H5future=" _rc
cap dmp H5future
cap asgdmp H5future

dis "PROBEH_DONE"
