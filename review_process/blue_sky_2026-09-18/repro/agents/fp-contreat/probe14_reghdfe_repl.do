clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 121212
set sortseed 121212
local D "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat"

use "`D'/p1like_u.dta", clear
gen double relyr = year - ttre
gen byte tf = (ttre<. & relyr>=-2 & relyr<=1)
egen long cell = group(st relyr) if tf
qui replace cell = 0 if cell==.
gen double ringer = dose_it
qui gquantiles lat = dose_it , xtile nquantiles(3)

dis "===== reghdfe, absorb(st year lat##c.ringer), robust ====="
cap noisily reghdfe y_dose ib0.cell##i.lat##c.dose_it , absorb(i.st i.year i.lat##c.(ringer)) vce(robust)
dis "RC=" _rc
if _rc==0 {
  mat V = e(V)
  mata: st_local("mx", strofreal(max(abs(st_matrix("V")))))
  dis "RGHDFE maxabsV=`mx'  df_r=" e(df_r) " rss=" e(rss) " N=" e(N) " rank=" e(rank) " df_a=" e(df_a)
}

dis "===== reghdfe, SAME but absorb only st year (controls explicit), robust ====="
cap noisily reghdfe y_dose i.lat##c.(ringer) ib0.cell##i.lat##c.dose_it , absorb(i.st i.year) vce(robust)
dis "RC=" _rc
if _rc==0 {
  mat V = e(V)
  mata: st_local("mx", strofreal(max(abs(st_matrix("V")))))
  dis "RGHDFE2 maxabsV=`mx'  df_r=" e(df_r) " rss=" e(rss) " N=" e(N) " rank=" e(rank) " df_a=" e(df_a)
}

dis "===== plain reg, everything explicit, robust ====="
cap noisily qui reg y_dose i.st i.year i.lat##c.(ringer) ib0.cell##i.lat##c.dose_it , vce(robust)
dis "RC=" _rc
if _rc==0 {
  mat V = e(V)
  mata: st_local("mx", strofreal(max(abs(st_matrix("V")))))
  dis "REG maxabsV=`mx'  df_r=" e(df_r) " rss=" e(rss) " N=" e(N) " rank=" e(rank)
}

dis "===== reghdfe WITHOUT the lattice (baseline), robust ====="
cap noisily reghdfe y_dose ib0.cell##c.dose_it , absorb(i.st i.year c.ringer) vce(robust)
dis "RC=" _rc
if _rc==0 {
  mat V = e(V)
  mata: st_local("mx", strofreal(max(abs(st_matrix("V")))))
  dis "RGHDFE3 maxabsV=`mx'  df_r=" e(df_r) " rss=" e(rss) " N=" e(N) " rank=" e(rank) " df_a=" e(df_a)
}
dis "DONE"
