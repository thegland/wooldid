* probe20: the bootstrap's re-id leg REFINES the partition, so a cohort-constant
* covariate stays collinear there -- the contrast that isolates FP-1 to the rebuild.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm"
use "dataB3.dta", clear
cap noi qui wooldid y coh t onset, att_i cluster(coh) bootstrap(standard) reps(99) rseed(4242)
dis "BOOT_PLAIN rc=" _rc " b=" %21x _b[main_att_i] " se=" %21x _se[main_att_i]
use "dataB3.dta", clear
cap noi qui wooldid y coh t onset, att_i cluster(coh) bootstrap(standard) reps(99) rseed(4242) fe(wave)
dis "BOOT_FEWAVE rc=" _rc " b=" %21x _b[main_att_i] " se=" %21x _se[main_att_i]
use "dataB3.dta", clear
cap noi qui wooldid y coh t onset, att_i cluster(coh) bootstrap(standard) reps(99) rseed(4242) controls(cohx)
dis "BOOT_CTRL rc=" _rc " b=" %21x _b[main_att_i] " se=" %21x _se[main_att_i]
