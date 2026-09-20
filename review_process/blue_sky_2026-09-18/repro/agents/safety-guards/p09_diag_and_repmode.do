* p09 -- (i) pin the si8 diagnosis: i.i is a REAL variable demand only on the reg path;
* NOTE FOR A VERIFIER: gen_core's `gcounts'/`gnevern' locals persist past `clear all',
*   so block (ii) below ran on the SAME 1,900-unit panel as block (i), not the 60-unit
*   panel the knobs there ask for.  The F13 results are valid on that data regardless.
*        (ii) the 2026-09-16 cost audit made F13 and si8 skip when wdrepmode==1.  Does
*             the OBSERVED run still fire under each engine?
clear all
set more off
set processors 2
set maxvar 2048
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

local gcounts 10 10
local gnevern 1880
local gT 8
local gt0 2000
local gonsets 2004 2005
local gysd 0.30
include "`GEN'/gen_core.do"
dis "MAXVAR=" c(maxvar)
dis "=== i: reghdfe path, safety(off), maxvar 2048 ==="
capture noisily wooldid ync unit t ttre, att cluster(unit) safety(off)
dis "I_REGHDFE_RC=" _rc
dis "=== i: reg path (customvce) -- i.i IS materialized there ==="
capture noisily wooldid ync unit t ttre, att customvce(cluster unit) safety(off)
dis "I_REGPATH_RC=" _rc

*---------------- (ii) wdrepmode skip of F13 / si8 ----------------
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
local gN 60
local gT 10
local gt0 2000
local gonsets 2004 2006
local gysd 0.30
include "`GEN'/gen_core.do"
* onset recorded only from the onset period onward: the classic F13 violation
gen double ttrebad = ttre
qui replace ttrebad = . if ttre < . & t < ttre
dis "=== F13 observed run, no engine ==="
capture noisily wooldid ync unit t ttrebad, att cluster(unit)
dis "F13_PLAIN_RC=" _rc
dis "=== F13 under bootstrap(standard) ==="
capture noisily wooldid ync unit t ttrebad, att cluster(unit) bootstrap(standard) reps(5)
dis "F13_BOOT_RC=" _rc
dis "=== F13 under permute() ==="
capture noisily wooldid ync unit t ttrebad, att cluster(unit) permute(ttre) permreps(5)
dis "F13_PERM_RC=" _rc
dis "=== F13 with safety(off): documented bypass, what number? ==="
capture noisily wooldid ync unit t ttrebad, att cluster(unit) safety(off)
dis "F13_OFF_RC=" _rc
mat F1 = e(wdidmainresults)
mat list F1, format(%10.7f)
capture noisily wooldid ync unit t ttre, att cluster(unit)
dis "F13_GOOD_RC=" _rc
mat F2 = e(wdidmainresults)
mat list F2, format(%10.7f)
dis "P09_DONE"
