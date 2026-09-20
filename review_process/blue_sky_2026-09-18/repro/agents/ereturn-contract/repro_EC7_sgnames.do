* EC-7: the help's subgroup coefficient-name suffix does not exist in e(b).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 120
gen long id = _n
gen byte grp = mod(id,5)
gen byte sg = mod(id,2)
gen double ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
replace ttre = 11 if grp==4
gen long clu = ceil(id/6)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + 0.25*t + 1.2*post + rnormal()
wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod ///
    esprelength(3) espostlength(3) subgroups(sg) cluster(clu)
dis "e(b) colnames:"
dis "`: colfullnames e(b)'"
dis "--- what sthlp:1306 tells the user to type ---"
cap noi lincom main_att_s1
dis "lincom main_att_s1  -> rc " _rc
dis "--- what the code produces (ado:4015) ---"
cap noi lincom main_SG1_att
dis "lincom main_SG1_att -> rc " _rc
dis "--- every OTHER subgroup-keyed name does use the documented _s# suffix ---"
foreach n in b_es_att_s1 escomp_s1 wdidsgresults_1 {
  cap confirm matrix e(`n')
  dis "  e(`n') exists as a matrix? " (_rc==0)
}
cap confirm scalar e(numpre_att_s1)
dis "  e(numpre_att_s1) exists as a scalar? " (_rc==0)
