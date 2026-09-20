* probe2.do -- reps() caution reachability; reps vs permreps; save-file name handling;
*              string cluster on the re-identification leg.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 222
local gN 40
local gT 7
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
* a realistic STRING cluster variable (state abbreviations)
qui gen str3 stabb = "S" + string(mod(unit,7)+1)
qui gen long state = mod(unit,7)+1
save base2.dta, replace

dis "=============== P5: reps() resolution caution, clustered vs not ==============="
use base2.dta, clear
dis "--- P5a CLUSTERED reps(12) (caution expected) ---"
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(12) rseed(9) dots(0)
dis "P5a_RC=" _rc
dis "--- P5b UNCLUSTERED reps(12) (is any caution printed?) ---"
use base2.dta, clear
cap noi wooldid ync unit t ttre, att bootstrap(standard) reps(12) rseed(9) dots(0)
dis "P5b_RC=" _rc
dis "--- P5c UNCLUSTERED reps(4) boottstat (grid 2/5 = .4) ---"
use base2.dta, clear
cap noi wooldid ync unit t ttre, att bootstrap(standard) boottstat reps(4) rseed(9) dots(0)
dis "P5c_RC=" _rc
dis "--- P5d UNCLUSTERED reps(4) bootdisplay(wald) ---"
use base2.dta, clear
cap noi wooldid ync unit t ttre, att bootstrap(standard) bootdisplay(wald) reps(4) rseed(9) dots(0)
dis "P5d_RC=" _rc

dis "=============== P8: reps()/permreps() resolution ==============="
use base2.dta, clear
dis "--- P8a bootstrap with permreps() too ---"
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(20) permreps(50) rseed(9) dots(0)
dis "P8a_RC=" _rc
dis "--- P8b permute with reps() too ---"
use base2.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) permute(shuffle) permreps(25) reps(20) rseed(9) dots(0)
dis "P8b_RC=" _rc
dis "--- P8c reps() with neither engine ---"
use base2.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) reps(20)
dis "P8c_RC=" _rc

dis "=============== P6: save-file name handling ==============="
cap erase "junkb.dta"
cap erase "junkb, replace.dta"
cap erase "junkp.dta"
cap erase "junkp, replace.dta"
use base2.dta, clear
dis "--- P6a bootsave(junkb, replace) : comma check expected (RF-6) ---"
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(9) dots(0) bootsave(junkb, replace)
dis "P6a_RC=" _rc
use base2.dta, clear
dis "--- P6b permsave(junkp, replace) : is there a comma check? ---"
cap noi wooldid ync unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(9) dots(0) permsave(junkp, replace)
dis "P6b_RC=" _rc
dis "--- files now on disk ---"
cap noi dir "junk*"
dis "P6_DIR_RC=" _rc
cap confirm file "junkp, replace.dta"
dis "P6b_BADFILE_EXISTS_RC=" _rc
cap confirm file "junkb, replace.dta"
dis "P6a_BADFILE_EXISTS_RC=" _rc

dis "=============== P9 (self-found): string cluster variable on the reid leg ==============="
use base2.dta, clear
dis "--- P9a analytic path, string cluster ---"
cap noi wooldid ync unit t ttre, att cluster(stabb)
dis "P9a_RC=" _rc
use base2.dta, clear
dis "--- P9b bootstrap COUNT leg, string cluster ---"
cap noi wooldid ync unit t ttre, att cluster(stabb) bootstrap(standard) reps(15) rseed(9) dots(0)
dis "P9b_RC=" _rc
use base2.dta, clear
dis "--- P9c bootstrap REID leg (att_i), string cluster ---"
cap noi wooldid ync unit t ttre, att att_i cluster(stabb) bootstrap(standard) reps(15) rseed(9) dots(0)
dis "P9c_RC=" _rc
use base2.dta, clear
dis "--- P9d permute, string cluster ---"
cap noi wooldid ync unit t ttre, att cluster(stabb) permute(shuffle) permreps(20) rseed(9) dots(0)
dis "P9d_RC=" _rc

dis "PROBE2_DONE"
