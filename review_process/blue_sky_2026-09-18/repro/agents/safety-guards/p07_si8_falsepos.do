* p07_si8_falsepos.do -- the si8 maxvar preflight (ado:1355-1377) is documented as a
* LOWER bound on the regression's variable demand, so "an early refusal can only fire
* when the late r(907) was certain anyway" (ado:1341-1348) and the design claims ZERO
* false positives.  But the demand formula adds si8ilev + si8tlev, and on the DEFAULT
* engine path i.i and i.t are ABSORBED by reghdfe (ado:1476), not materialized.
* A wide, mostly-never-treated panel should therefore be refused although it runs.
clear all
set more off
set processors 2
set maxvar 2048          // Stata/BE's fixed ceiling, and Stata's own minimum
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

* 1900 counties, 8 years, only 20 of them ever treated -- an ordinary county panel
local gcounts 10 10
local gnevern 1880
local gT 8
local gt0 2000
local gonsets 2004 2005
local gysd 0.30
include "`GEN'/gen_core.do"
qui levelsof unit, local(uu)
dis "NUNITS=" wordcount("`uu'")
qui count
dis "NROWS=" r(N)
dis "MAXVAR=" c(maxvar)
save "p07.dta", replace

dis "=== default safety: the preflight ==="
capture noisily wooldid ync unit t ttre, att cluster(unit)
dis "P7_DEF_RC=" _rc

dis "=== safety(off): does the run the preflight called impossible actually run? ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) safety(off)
dis "P7_OFF_RC=" _rc
cap mat list e(wdidmainresults), format(%10.7f)
dis "P7_OFF_N=" e(N)
dis "P7_DONE"
