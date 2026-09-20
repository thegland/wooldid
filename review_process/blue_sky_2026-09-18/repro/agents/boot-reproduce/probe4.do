* probe4.do -- EXTERNAL ORACLE: sthlp:532 claims the count leg "matches Stata's own
* bootstrap prefix replication for replication under a common seed". Test it.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 444
local gN 60
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
include "`GEN'/gen_core.do"
save base4.dta, replace

local R  40
local SD 2468

cap erase wdraws.dta
cap erase prefdraws.dta

* --- wooldid's own count leg
use base4.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(wdraws) replace
dis "WD_RC=" _rc
mat bW = e(b)
dis "WD_colnames: " "`: colfullnames bW'"

* --- Stata's own bootstrap prefix, help's recipe verbatim
use base4.dta, clear
cap noi bootstrap attb = _b[main_att], reps(`R') seed(`SD') cluster(unit) idcluster(newc) saving(prefdraws, replace) nodots: wooldid ync unit t ttre, att fastnose
dis "PF_RC=" _rc
cap noi ereturn list
dis "PF_LIST_RC=" _rc

* --- compare replicate for replicate
cap noi {
  use prefdraws.dta, clear
  qui gen long _rep = _n
  qui keep _rep attb
  tempfile tp
  qui save `tp'
  use wdraws.dta, clear
  qui keep if _rep > 0
  qui merge 1:1 _rep using `tp', nogen
  qui gen double d = reldif(main_att, attb)
  qui sum d, meanonly
  dis "PREFIX_VS_COUNTLEG MAXRELDIF = " %12.4e r(max)
  qui count if d > 1e-10
  dis "PREFIX_VS_COUNTLEG n_differ(>1e-10) = " r(N) " of " _N
  list _rep main_att attb in 1/10, noobs
}
dis "CMP_RC=" _rc
dis "PROBE4_DONE"
