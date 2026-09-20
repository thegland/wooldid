* ==========================================================================
* r03_names.do -- independent test of claim fpsg-02.
* sthlp:1307: "Under subgroups() each subgroup's rows take the suffix _s#."
* Does a user who follows that sentence get a working lincom/test name?
* ==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 24
gen long id = _n
gen double ttre = .
replace ttre = 7 if id <= 16
gen byte sg = .
replace sg = 1  if inrange(id, 1, 4)
replace sg = 3  if inrange(id, 5, 8)
replace sg = 7  if inrange(id, 9, 12)
replace sg = 99 if inrange(id, 13, 16)
replace sg = 1  if id > 16
gen double eff = cond(sg==1,1,cond(sg==3,2,cond(sg==7,3,4)))
expand 12
bysort id: gen int t = _n
set seed 4242
gen double y = 2*id + 0.4*t + eff*(t >= ttre & ttre < .) + rnormal(0,0.3)

capture noisily wooldid y id t ttre, att att_it att_i att_itime subgroups(sg) esprelength(2) espostlength(2)
dis "RUN_RC=" _rc

dis "=== e(b) COLNAMES (the naming contract of sthlp:1307) ==="
matrix bb = e(b)
local cn : colfullnames bb
dis "COLNAMES: `cn'"
dis "NCOLS=" colsof(bb)

dis "=== what a user who followed sthlp:1307 types ==="
foreach nm in main_att_s1 main_att_s2 main_att_it_s1 pretreatment_att_s1 ES_att_post0_s1 main_att_s99 {
  capture noisily lincom `nm'
  dis "HELPNAME `nm' rc=" _rc
}
dis "=== the names that actually work ==="
foreach nm in main_SG1_att main_SG3_att main_SG99_att ES_att_SG1_post0 pretreatment_SG1_att {
  capture noisily lincom `nm'
  dis "TRUENAME `nm' rc=" _rc
}
dis "=== does _s<value> appear anywhere as an e() NAME? ==="
foreach en in b_es_att_s1 b_es_att_s99 escomp_s3 numpre_att_s99 wdidsgresults_1 {
  capture noisily matrix zz = e(`en')
  dis "ENAME `en' rc=" _rc
}
dis "=== literal grep of the help for the true convention ==="
dis "ALLDONE"
