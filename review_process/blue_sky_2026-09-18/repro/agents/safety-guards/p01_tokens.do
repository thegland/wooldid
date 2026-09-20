* p01_tokens.do -- assert the pinned build, then characterise safety() token handling
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

which wooldid
dis "BUILDCHECK_DONE"

* ---------- clean panel, exact DGP ----------
local gN 60
local gT 10
local gt0 2000
local gonsets 2004 2006
local gysd 0
local gaisd 0
local gtrendsd 0
include "`GEN'/gen_core.do"
save "p01_clean.dta", replace

* baseline: default safety
capture noisily wooldid ync unit t ttre, att cluster(unit)
dis "T_DEFAULT_RC=" _rc
mat b0 = e(b)
dis "T_DEFAULT_ATT=" %21x e(wdid_att)

foreach tok in on off ON OFF Off oFF "off " " off" 0 1 yes no {
  capture noisily wooldid ync unit t ttre, att cluster(unit) safety(`tok')
  dis "TOKEN[`tok']_RC=" _rc
}
dis "P01_TOKENS_DONE"

* ---------- now a panel that TRIPS a guard, to see which message wins ----------
* no never-treated: every unit treated inside the window
clear
local gN 40
local gT 8
local gt0 2000
local gonsets 2003 2005
local gnever 0
local gysd 0
local gaisd 0
include "`GEN'/gen_core.do"
count if ttre >= .
dis "NEVERTREATED_ROWS=" r(N)
save "p01_nonever.dta", replace

capture noisily wooldid ync unit t ttre, att cluster(unit)
dis "NN_DEFAULT_RC=" _rc
capture noisily wooldid ync unit t ttre, att cluster(unit) safety(OFF)
dis "NN_OFFCAPS_RC=" _rc
capture noisily wooldid ync unit t ttre, att cluster(unit) safety(off)
dis "NN_OFF_RC=" _rc
dis "P01_DONE"
