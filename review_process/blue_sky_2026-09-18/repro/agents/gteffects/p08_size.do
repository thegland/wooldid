* p08_size.do -- does the gteffects table survive a LARGE cell count?
* ado:1765 builds the per-cell weight block with `cap mkmat ... , matrix(gtraw)'.
* mkmat is capped by matsize / the matrix row limit; on failure `gtn' stays 0 and the
* user is told "gteffects found no cohort-by-relative-period main-effect coefficients
* to report", which would be a FALSE statement about the fit.
clear all
set more off
set processors 2
set maxvar 30000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
set seed 20260916
set sortseed 20260916
local GD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/gteffects"
include "`GD'/gtlib.do"
local EST att att_it att_i att_itime

dis "MATSIZE=" c(matsize) " MAXMATSIZE=" c(max_matsize)
* direct mkmat probe, no wooldid involved
clear
qui set obs 20000
gen double a = _n
gen double b = _n*2
foreach k in 400 401 800 5000 11000 11001 20000 {
  preserve
  qui keep in 1/`k'
  tempname M
  cap mkmat a b, matrix(`M')
  dis "MKMAT_`k'_RC=" _rc
  restore
}

* --- wooldid with a cell count above 400 ---
mkpanel , n(140) t(9) onsets("2004 2006 2008") ysd(0.4) seed(313)
qui count if ttre < .
dis "P8 treated obs=" r(N)
ereturn clear
cap noi wooldid y unit year ttre, `EST' gteffects esprelength(2) espostlength(2) safety(off)
dis "P8a_RC=" _rc
cap confirm matrix e(gteffects)
dis "P8a_HASTABLE=" (_rc==0)
cap noi dis "P8a_ROWS=" rowsof(e(gteffects))
gtcheck , tag(P8a_bigcells) pre es

dis "DONE_P08"
