* repro_br3.do -- SELF-CONTAINED.
* Without cluster(), `bootstrap(standard) att' and `bootstrap(standard) att att_i' at ONE
* rseed() draw DIFFERENT resamples, so att's reported standard error moves when a second
* estimand is added.  With cluster() the two legs draw the SAME resamples and att does not move.
* Mechanism: the count leg (ado:9295) and the physical-expansion leg (ado:9238) both draw with
*   bsample, cluster(<rowid>) idcluster(...)
* but the re-identification leg (ado:9345) uses a PLAIN `bsample' when cluster() is absent.
* The code comment at ado:9340 claims "same-seed draws equal the count path's".
clear all
set more off
set processors 1
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 1
set sortseed 1
* --- a tiny staggered panel, built here, no repo data ---
clear
set obs 40
gen long unit = _n
gen long ttre = cond(unit <= 12, 2003, cond(unit <= 24, 2005, .))
expand 8
bysort unit: gen int t = 1999 + _n
set seed 7
gen double y = 0.4*unit/40 + 0.1*t + 1.5*(t >= ttre) + rnormal()
tempfile d
save `d'

cap prog drop shse
prog def shse
  syntax , tag(string)
  tempname b v
  mat `b' = e(b)
  mat `v' = e(V)
  local k = colnumb(`b',"main_att")
  dis "`tag' : b_att = " %21x `b'[1,`k'] "   se_att = " %21x sqrt(`v'[`k',`k']) "  = " %9.6f sqrt(`v'[`k',`k'])
end

foreach f in r1 r2 r3 r4 {
  cap erase `f'.dta
}

dis "================ NO cluster(): the two legs, ONE seed ================"
use `d', clear
qui wooldid y unit t ttre, att bootstrap(standard) reps(200) rseed(12345) dots(0) bootsave(r1) replace
shse , tag("UNCLUSTERED  att            (count leg)")
use `d', clear
qui wooldid y unit t ttre, att att_i bootstrap(standard) reps(200) rseed(12345) dots(0) bootsave(r2) replace
shse , tag("UNCLUSTERED  att att_i      (reid leg) ")

dis "================ WITH cluster(unit): the two legs, ONE seed ================"
use `d', clear
qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(200) rseed(12345) dots(0) bootsave(r3) replace
shse , tag("CLUSTERED    att            (count leg)")
use `d', clear
qui wooldid y unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(200) rseed(12345) dots(0) bootsave(r4) replace
shse , tag("CLUSTERED    att att_i      (reid leg) ")

dis "================ are the DRAWS the same? ================"
use r1.dta, clear
qui keep _rep main_att
qui rename main_att att_countleg
tempfile a
qui save `a'
use r2.dta, clear
qui keep _rep main_att
qui rename main_att att_reidleg
qui merge 1:1 _rep using `a', nogen
qui gen double d = reldif(att_countleg, att_reidleg)
qui sum d if _rep > 0, meanonly
dis "UNCLUSTERED: max relative difference between the two legs' replicate draws = " %12.4e r(max)
list _rep att_countleg att_reidleg in 1/6, noobs

use r3.dta, clear
qui keep _rep main_att
qui rename main_att att_countleg
qui save `a', replace
use r4.dta, clear
qui keep _rep main_att
qui rename main_att att_reidleg
qui merge 1:1 _rep using `a', nogen
qui gen double d = reldif(att_countleg, att_reidleg)
qui sum d if _rep > 0, meanonly
dis "CLUSTERED  : max relative difference between the two legs' replicate draws = " %12.4e r(max)
list _rep att_countleg att_reidleg in 1/6, noobs
dis "REPRO_BR3_DONE"
