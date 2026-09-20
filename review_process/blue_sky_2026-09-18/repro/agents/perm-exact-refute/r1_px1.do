* r1_px1.do -- SKEPTIC probe of PERMEX-PX-1.
* Self-contained. Builds its own 5-cluster staggered panel, runs
*   (a) permute(shuffle) permexact permstepdown   -> M = 30 enumerated assignments
*   (b) permute(shuffle) permreps(199) permstepdown -> sampled draws
* and asks, for each, whether the step-down kernel's p rows are multiples of 1/n
* (plain share) or of 1/(n+1) (add-one), and what the footer/help say they are.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

* ---------- data: 5 units, T = 1..8, onsets 4,4,6,.,. -> M = 5!/(2!1!2!) = 30
clear
set obs 5
gen long unit = _n
gen double onset = .
replace onset = 4 in 1
replace onset = 4 in 2
replace onset = 6 in 3
gen double ue = rnormal()*0.5
expand 8
bysort unit: gen int year = _n
gen double te = rnormal()*0.3
gen byte post = (year >= onset) & onset < .
gen double y = 1 + 0.15*year + ue + te + 0.9*post
tempfile dat
save `dat', replace

dis "===== (a) permexact ====="
capture noisily wooldid y unit year onset, att espostlength(3) permute(shuffle) permexact permstepdown dots(0)
dis "A_RC=" _rc
dis "A_perm_exact=" e(perm_exact)
dis "A_perm_irebuilt=" e(perm_irebuilt)
matrix A_sd = e(perm_stepdown)
matrix list A_sd, format(%21x)
matrix list A_sd
matrix A_n = e(perm_n)
matrix list A_n
* independent check: is every p row a multiple of 1/30 or of 1/31?
local M = 30
forvalues c = 1/`=colsof(A_sd)' {
  forvalues r = 2/4 {
    local v = A_sd[`r',`c']
    local k30 = `v'*`M'
    local k31 = `v'*(`M'+1)
    dis "A_cell r=`r' c=`c' p=" %20.17f `v' "  p*30=" %12.9f `k30' "  p*31=" %12.9f `k31'
  }
}

dis "===== (b) sampled draws, permreps(199) ====="
use `dat', clear
capture noisily wooldid y unit year onset, att espostlength(3) permute(shuffle) permreps(199) permstepdown dots(0)
dis "B_RC=" _rc
dis "B_perm_exact=" e(perm_exact)
matrix B_sd = e(perm_stepdown)
matrix list B_sd
matrix B_n = e(perm_n)
matrix list B_n
forvalues c = 1/`=colsof(B_sd)' {
  forvalues r = 2/4 {
    local v = B_sd[`r',`c']
    dis "B_cell r=`r' c=`c' p=" %20.17f `v' "  p*199=" %12.9f `v'*199 "  p*200=" %12.9f `v'*200
  }
}
dis "DONE_R1"
