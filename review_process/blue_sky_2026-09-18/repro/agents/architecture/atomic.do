clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local BD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/architecture"

dis "=== A: post-regression failure (1 cluster) ==="
use "`BD'/base.dta", clear
gen byte onecl = 1
sort tt id
cap noi wooldid y id tt ttre , att cluster(onecl)
dis "ATOM A rc=" _rc
qui est dir
dis "ATOM A estdir=[`r(names)']"
dis "ATOM A nvars=" c(k) " nobs=" c(N) " sortedby=`: sortedby'"
cap confirm variable _ppmlhdfe_d
dis "ATOM A ppmld_present=" (_rc==0)

dis "=== B: verbose run retains ___wooldidfullmodel? ==="
use "`BD'/base.dta", clear
cap noi wooldid y id tt ttre , att verbose
dis "ATOM B rc=" _rc
qui est dir
dis "ATOM B estdir=[`r(names)']"

dis "=== C: reghdfe singleton drop, default vs nofastcollapse ==="
use "`BD'/base.dta", clear
* make 3 units appear only once -> reghdfe singleton drop in the i dimension
drop if inlist(id,2,5,9) & tt != 4
qui count
dis "ATOM C rows=" r(N)
cap noi wooldid y id tt ttre , att att_it att_i att_itime
dis "ATOM C1 rc=" _rc " N=" e(N) " Ndrop=" e(N_dropped)
foreach s in att att_it att_i att_itime {
  dis "ATOM C1 `s' = " %21x _b[main_`s']
}
cap noi wooldid y id tt ttre , att att_it att_i att_itime nofastcollapse
dis "ATOM C2 rc=" _rc " N=" e(N) " Ndrop=" e(N_dropped)
foreach s in att att_it att_i att_itime {
  dis "ATOM C2 `s' = " %21x _b[main_`s']
}
dis "ATOMDONE"
