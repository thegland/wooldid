* probe13.do -- SELF-FOUND: the bootstrap's resampling universe is `smp' (user-if & nonmissing
* cluster), NOT e(sample).  On a panel with missing outcomes -- an ordinary administrative
* shape -- rows and whole units sit in `smp' but not in e(sample).  Does the count leg still
* match Stata's own bootstrap prefix "replication for replication under a common seed"
* (sthlp:532)?  And does bootdf(clusters) use the right G?
* Also: bootsave()/permsave() must refuse a pre-existing file BEFORE the replications run.
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
* one unit whose outcome is entirely missing: a cluster that exists in `smp' and cannot
* exist in e(sample)
qui replace ync = . if unit == 5
qui count
dis "NOBS=" r(N)
qui count if !missing(ync)
dis "NOBS_NONMISS_Y=" r(N)
save base13.dta, replace

local R  40
local SD 1717
cap erase d1.dta
cap erase d2.dta
cap erase pref13.dta

dis "=========== plain analytic run: how many clusters does the ESTIMATION see? ==========="
use base13.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit)
dis "PL_RC=" _rc
dis "PLAIN|N=" e(N) "|N_clust=" e(N_clust)
qui count if e(sample)
dis "PLAIN|e(sample) rows=" r(N)

dis "=========== bootstrap count leg ==========="
use base13.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(d1) replace
dis "B_RC=" _rc
dis "BOOT|N=" e(N) "|boot_df=" e(boot_df) "|boot_completed=" e(boot_completed) "|boot_failed=" e(boot_failed)
mat bB = e(b)
mat vB = e(V)
dis "BOOT|b_att=" %21x bB[1,colnumb(bB,"main_att")] "|se_att=" %21x sqrt(vB[colnumb(bB,"main_att"),colnumb(bB,"main_att")])

dis "=========== Stata's own bootstrap prefix, same seed ==========="
use base13.dta, clear
cap noi bootstrap attb = _b[main_att], reps(`R') seed(`SD') cluster(unit) idcluster(newc) saving(pref13, replace) nodots: wooldid ync unit t ttre, att fastnose
dis "PF_RC=" _rc
dis "PREFIX|N_reps=" e(N_reps) "|N_misreps=" e(N_misreps) "|N_clust=" e(N_clust)

dis "=========== compare the draws ==========="
cap noi {
  use pref13.dta, clear
  qui gen long _rep = _n
  qui keep _rep attb
  tempfile tp
  qui save `tp'
  use d1.dta, clear
  qui keep if _rep > 0
  qui merge 1:1 _rep using `tp', nogen
  qui gen double d = reldif(main_att, attb)
  qui sum d, meanonly
  dis "P13|PREFIX_VS_COUNTLEG MAXRELDIF = " %12.4e r(max)
  qui count if d > 1e-6
  dis "P13|n_differ(>1e-6) = " r(N) " of " _N
  list _rep main_att attb in 1/12, noobs
}
dis "P13CMP_RC=" _rc

dis "=========== how many w_ columns did bootsave write? ==========="
use d1.dta, clear
qui ds w_*
dis "P13|w_ columns = " wordcount("`r(varlist)'")
cap confirm variable w_5
dis "P13|w_5 (the all-missing-y unit) present rc=" _rc
qui sum w_5 if _rep > 0
dis "P13|w_5 mean multiplicity = " r(mean) " (a phantom cluster if it is drawn at all)"

dis "=========== UNCLUSTERED leg on the same data ==========="
use base13.dta, clear
cap noi wooldid ync unit t ttre, att bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(d2) replace
dis "U_RC=" _rc
dis "UBOOT|N=" e(N) "|boot_completed=" e(boot_completed)
mat bU = e(b)
mat vU = e(V)
dis "UBOOT|se_att=" %21x sqrt(vU[colnumb(bU,"main_att"),colnumb(bU,"main_att")])

dis "=========== pre-existing save file without replace must be r(602), BEFORE the loop ==========="
use base13.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(200) rseed(1) dots(1) bootsave(d1)
dis "P13|BOOTSAVE_EXISTING_RC=" _rc
use base13.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) permute(shuffle) permreps(200) rseed(1) dots(1) permsave(d1)
dis "P13|PERMSAVE_EXISTING_RC=" _rc

dis "PROBE13_DONE"
