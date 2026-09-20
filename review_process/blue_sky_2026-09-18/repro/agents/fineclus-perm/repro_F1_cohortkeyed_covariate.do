*==============================================================================
* SELF-CONTAINED REPRO  --  FINDING F1
* permonsetcohorts: a covariate keyed to the OBSERVED cohort is silently dropped
* from the observed fit (perfectly collinear with the i fixed effects) yet is LIVE
* inside every replication, because wooldid_permrun substitutes the rebuilt i into
* the VARLIST only (ado:7639) while repopts still names the original variables
* (ado:125-163).  Result: an option that provably cannot move the point estimate
* halves the permutation p-value.
*
* Run:  stata-mp -e do repro_F1_cohortkeyed_covariate.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

*---- a ragged administrative-shaped panel: 6 states, unequal caseloads and windows
set seed 909090
set sortseed 909090
set obs 6
gen int state = _n
gen double sfx = rnormal()*2
gen int nper = cond(state==1,2,cond(state==2,3,cond(state==3,1,cond(state==4,2,cond(state==5,3,1)))))
gen int tlo  = cond(state==2,2,cond(state==5,3,1))
gen int thi  = cond(state==3,7,8)
expand 8
bysort state: gen int t = _n
drop if t < tlo | t > thi
expand nper
bysort state t: gen int unit = _n
gen double ufx = rnormal()
gen double eps = rnormal()
* i IS the onset grouping, as permonsetcohorts asserts; cluster() = state is FINER
gen int coh   = cond(inlist(state,1,2),1,cond(inlist(state,3,4),2,3))
gen double onset = cond(coh==1,4,cond(coh==2,6,.))
gen byte D = (onset<. & t>=onset)
gen double y = 1 + sfx + ufx + 0.3*t + 1.5*D + eps
*---- the contaminating covariate: an "adoption wave" index. It is constant within
*---- the observed cohort, so it is EXACTLY collinear with the cohort fixed effects
gen double wave = 10*coh + 0.5*coh*coh
*---- and the same idea as a fixed effect: an "adoption wave" grouping (early vs rest)
gen byte wavefe = (coh <= 2)
keep state unit t coh onset y wave wavefe
save "_F1.dta", replace

*---- (1) the covariate provably cannot move the point estimate --------------
use "_F1.dta", clear
qui wooldid y coh t onset, att cluster(state) fastnose
scalar b0 = _b[main_att]
use "_F1.dta", clear
qui wooldid y coh t onset, att cluster(state) fastnose controls(wave)
scalar b1 = _b[main_att]
dis "F1 point estimate WITHOUT controls(wave): " %21x b0
dis "F1 point estimate WITH    controls(wave): " %21x b1
dis "F1 point estimates identical?  " cond(b0 == b1, "YES (bitwise)", "NO")

*---- (2) yet the permutation p-value halves ---------------------------------
use "_F1.dta", clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts
dis "F1 BASELINE   p = " e(perm_p_att) "   null sd = " e(perm_sd_att) "   b = " %21x _b[main_att]
scalar pA  = e(perm_p_att)
scalar sdA = e(perm_sd_att)

use "_F1.dta", clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts controls(wave)
dis "F1 +ctrl(wave) p = " e(perm_p_att) "   null sd = " e(perm_sd_att) "   b = " %21x _b[main_att]
scalar pB  = e(perm_p_att)
scalar sdB = e(perm_sd_att)

dis " "
dis "==> p BASELINE = " pA "   p WITH A PROVABLY-INERT OPTION = " pB
dis "==> null-distribution SD shrinks from " sdA " to " sdB
dis "==> ratio p(B)/p(A) = " pB/pA

*---- (2b) the same thing as an adoption-wave FIXED EFFECT -------------------
use "_F1.dta", clear
qui wooldid y coh t onset, att cluster(state) fastnose fe(wavefe)
dis "F1b point estimate WITH fe(wavefe): " %21x _b[main_att] "   (baseline " %21x b0 ")"
use "_F1.dta", clear
cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permexact permonsetcohorts fe(wavefe)
dis "F1b +fe(wavefe) p = " e(perm_p_att) "   null sd = " e(perm_sd_att)

*---- (3) CONTROL ARM: cluster() AT i, so there is no rebuild. The same inert
*----     option must then leave the p-value alone -- and it does.
use "_F1.dta", clear
cap noi wooldid y coh t onset, att cluster(coh) permute(shuffle) permexact
dis "F1 CONTROL (cluster at i) baseline    p = " e(perm_p_att)
scalar pC = e(perm_p_att)
use "_F1.dta", clear
cap noi wooldid y coh t onset, att cluster(coh) permute(shuffle) permexact controls(wave)
dis "F1 CONTROL (cluster at i) +ctrl(wave) p = " e(perm_p_att)
scalar pD = e(perm_p_att)
dis "==> control arm p identical? " cond(pC == pD, "YES", "NO")

*---- (4) mechanism: inside a replication the covariate is no longer collinear.
*----     Hand-build one replication (states 1,3 -> onset 4; 2,5 -> onset 6; 4,6 NT)
*----     exactly as the ado does (rebuilt i takes the observed cohort's code).
use "_F1.dta", clear
gen double repons = cond(inlist(state,1,3),4,cond(inlist(state,2,5),6,.))
gen double repi   = cond(repons==4,1,cond(repons==6,2,3))
qui wooldid y repi t repons, att fastnose cluster(state) esrelativeto(-1)
scalar r_plain = _b[main_att]
qui wooldid y repi t repons, att fastnose cluster(state) esrelativeto(-1) controls(wave)
scalar r_wave = _b[main_att]
dis "F1 one replication, no covariate    : " %21x r_plain
dis "F1 one replication, controls(wave)  : " %21x r_wave
dis "F1 replication estimates differ by  : " r_wave - r_plain
dis "    (in the OBSERVED run the same covariate changes nothing at all)"
