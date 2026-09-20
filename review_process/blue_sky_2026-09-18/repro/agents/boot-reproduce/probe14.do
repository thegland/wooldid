* probe14.do -- nail the cause of the probe13 mismatch: is it the phantom cluster?
* Also: does e(boot_df) (default bootdf(clusters) = G-1) count clusters the estimation dropped?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"
which wooldid

local gseed 1414
local gN 60
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
local gymiss 0.15
include "`GEN'/gen_core.do"
qui replace ync = . if unit == 5
save base14.dta, replace

local R  40
local SD 1717

dis "=========== ARM A: the phantom unit is IN the data ==========="
use base14.dta, clear
cap erase a1.dta
cap erase prefA.dta
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(a1) replace
dis "A_RC=" _rc
dis "ARMA|e(N)=" e(N) "|e(N_clust)=" e(N_clust) "|e(boot_df)=" e(boot_df)
use base14.dta, clear
cap noi bootstrap attb = _b[main_att], reps(`R') seed(`SD') cluster(unit) idcluster(newc) saving(prefA, replace) nodots: wooldid ync unit t ttre, att fastnose
dis "ARMA_PREFIX|e(N_clust)=" e(N_clust)
cap noi {
  use prefA.dta, clear
  qui gen long _rep = _n
  qui keep _rep attb
  tempfile tp
  qui save `tp'
  use a1.dta, clear
  qui keep if _rep > 0
  qui merge 1:1 _rep using `tp', nogen
  qui gen double d = reldif(main_att, attb)
  qui sum d, meanonly
  dis "ARMA|PREFIX_VS_COUNTLEG MAXRELDIF = " %12.4e r(max)
}
dis "ARMA_CMP_RC=" _rc

dis "=========== ARM B: the phantom unit is DROPPED from the data ==========="
use base14.dta, clear
qui drop if unit == 5
cap erase b1.dta
cap erase prefB.dta
save base14b.dta, replace
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(b1) replace
dis "B_RC=" _rc
dis "ARMB|e(N)=" e(N) "|e(N_clust)=" e(N_clust) "|e(boot_df)=" e(boot_df)
use base14b.dta, clear
cap noi bootstrap attb = _b[main_att], reps(`R') seed(`SD') cluster(unit) idcluster(newc) saving(prefB, replace) nodots: wooldid ync unit t ttre, att fastnose
dis "ARMB_PREFIX|e(N_clust)=" e(N_clust)
cap noi {
  use prefB.dta, clear
  qui gen long _rep = _n
  qui keep _rep attb
  tempfile tp
  qui save `tp'
  use b1.dta, clear
  qui keep if _rep > 0
  qui merge 1:1 _rep using `tp', nogen
  qui gen double d = reldif(main_att, attb)
  qui sum d, meanonly
  dis "ARMB|PREFIX_VS_COUNTLEG MAXRELDIF = " %12.4e r(max)
}
dis "ARMB_CMP_RC=" _rc

dis "=========== ARM C: a cluster variable with MISSING values ==========="
use base14.dta, clear
qui gen long clx = unit
qui replace clx = . if unit <= 3
cap noi wooldid ync unit t ttre, att cluster(clx)
dis "ARMC_PLAIN|e(N)=" e(N) "|e(N_clust)=" e(N_clust)
use base14.dta, clear
qui gen long clx = unit
qui replace clx = . if unit <= 3
cap noi wooldid ync unit t ttre, att cluster(clx) bootstrap(standard) reps(`R') rseed(`SD') dots(0)
dis "ARMC_BOOT|e(N)=" e(N) "|e(N_clust)=" e(N_clust) "|e(boot_df)=" e(boot_df)

dis "PROBE14_DONE"
