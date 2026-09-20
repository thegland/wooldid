clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
clear
set obs 45
gen long id = _n
gen byte cohortno = ceil(id/5) if id <= 30
gen int ttre = 2003 + cohortno if cohortno < .
gen double effect = 7 - cohortno if cohortno < .
replace effect = 0 if effect >= .
gen double unitfe = 0.37*id - 3
expand 12
bysort id: gen int year = 1999 + _n
gen double yearfe = 0.21*(year-2000)
gen byte post = (ttre < .) & (year >= ttre)
gen double y = unitfe + yearfe + effect*post
gen byte sector = 0
replace sector = 1 if inlist(id,1,6,11,16,21,26)
replace sector = mod(id,2) if id > 30
gen byte blank = (id==1 & year>=2005) | (id==6 & year>=2007)
count if blank
dis "BLANKED_ROWS=" r(N)
count
dis "TOTAL_ROWS=" r(N)
gen byte sector_miss = sector
replace sector_miss = . if blank
dis "===EREBEGIN==="
qui wooldid y id year ttre, att att_i att_it att_itime subgroups(sector_miss)
ereturn list
dis "===EREEND==="
