* ============================================================================
* SKEPTIC FOLLOW-UP for hv-noparent-F3 (+ two F1 loose ends)
* Part 1: what does `margins, eydx()` actually compute when fitted values
*         straddle zero -- mean of ratios, or ratio of means? (standalone, no
*         wooldid; settles whether a returned number is "complete but
*         meaningless" or merely "correct arithmetic on a bad estimand")
* Part 2: wooldid on a panel whose TREATED-POST fitted values straddle zero,
*         with the weight syntax CORRECT this time ([aw=w] before the comma).
* Part 3: F1 loose ends -- the cfxplottypes validator with a full parent, and
*         ccnoabsorb under contreat() with no contreatcontrols().
* run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r_f3b_semantics.do
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150

which wooldid

* ================= PART 1: margins eydx semantics, no wooldid ===============
clear
set obs 400
gen int g = mod(_n,20)
gen byte x = (_n>200)
gen double u = -3 + 6*runiform()
gen double y = u + 1.0*x + rnormal()*0.05
qui reg y x i.g
predict double yh, xb
qui su yh if x==1
scalar mean_yh = r(mean)
scalar min_yh  = r(min)
scalar max_yh  = r(max)
gen double rat = _b[x]/yh
qui su rat if x==1
scalar mean_of_ratios = r(mean)
scalar ratio_of_means = _b[x]/mean_yh
dis "ORACLE_MINYHAT=" min_yh
dis "ORACLE_MAXYHAT=" max_yh
dis "ORACLE_MEANOFRATIOS=" mean_of_ratios
dis "ORACLE_RATIOOFMEANS=" ratio_of_means
capture noisily margins, eydx(x) subpop(if x==1) noestimcheck
dis "ORACLE_MARGINS_RC=" _rc
capture noisily matrix list r(table)

* a positive-yhat control, same model shape
replace y = y + 10
qui reg y x i.g
capture drop yh
predict double yh, xb
capture drop rat
gen double rat = _b[x]/yh
qui su rat if x==1
dis "ORACLE_POS_MEANOFRATIOS=" r(mean)
qui su yh if x==1
dis "ORACLE_POS_RATIOOFMEANS=" _b[x]/r(mean)
dis "ORACLE_POS_MINYHAT=" r(min)
capture noisily margins, eydx(x) subpop(if x==1) noestimcheck
dis "ORACLE_POS_RC=" _rc
capture noisily matrix list r(table)

* ================= PART 2: wooldid, straddling panel, weights right =========
clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = -4 + 8*runiform()
expand 10
bysort unit: gen int t = _n
gen double y = ufe + 0.05*t + rnormal()*0.2
replace y = y + 0.8*(t>=ttre) if ttre<.
gen double w = 1 + runiform()
tempfile Q
save `Q', replace

use `Q', clear
qui su y if t>=ttre & ttre<.
dis "Q_MEANY_TP=" r(mean)
dis "Q_MINY_TP=" r(min)
gen double rr = .8/y
qui su rr if t>=ttre & ttre<.
dis "Q_MEANOFRATIOS_approx=" r(mean)

use `Q', clear
capture noisily wooldid y unit t ttre, att robust
dis "QB0_RC=" _rc

use `Q', clear
capture noisily wooldid y unit t ttre, att robust semielasticity
dis "QB1_RC=" _rc

use `Q', clear
capture noisily wooldid y unit t ttre [aw=w], att robust semielasticity
dis "QB2_AW_RC=" _rc

use `Q', clear
capture noisily wooldid y unit t ttre [pw=w], att robust semielasticity
dis "QB3_PW_RC=" _rc

use `Q', clear
capture noisily wooldid y unit t ttre [aw=w], att robust
dis "QB4_CONTROL_AW_RC=" _rc

* a tighter straddle -- fitted values right on top of zero
clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = rnormal()*0.3
expand 10
bysort unit: gen int t = _n
gen double y = ufe + 0.01*t + rnormal()*0.05
replace y = y + 0.2*(t>=ttre) if ttre<.
gen double w = 1 + runiform()
tempfile R
save `R', replace
use `R', clear
qui su y if t>=ttre & ttre<.
dis "R_MINY_TP=" r(min)
dis "R_MEANY_TP=" r(mean)
use `R', clear
capture noisily wooldid y unit t ttre, att robust semielasticity
dis "RB1_RC=" _rc
use `R', clear
capture noisily wooldid y unit t ttre [aw=w], att robust semielasticity
dis "RB2_AW_RC=" _rc

* ================= PART 3: F1 loose ends ====================================
clear
set obs 80
gen long unit = _n
gen byte grp = mod(_n,4)
gen ttre = .
replace ttre = 5  if grp==1
replace ttre = 7  if grp==2
replace ttre = 9  if grp==3
gen double ufe = rnormal()
gen double dose0 = runiform()*3
expand 10
bysort unit: gen int t = _n
gen double y = 2 + ufe + 0.1*t + rnormal()
replace y = y + 0.8*(t>=ttre) if ttre<.
gen double dose = dose0*(t>=ttre & ttre<.)
tempfile C
save `C', replace

* cfxplottypes validator WITH the full parent (makecfxplots + a lattice)
use `C', clear
capture noisily wooldid y unit t ttre, att robust contreat(dose) contreatcontrols(untreatedZero) ///
    makecfxplots cfxlattice(dose0) cfxplottypes(psot) safety(off)
dis "VG6b_fullparent_RC=" _rc

* same request, parent ABSENT, same data
use `C', clear
capture noisily wooldid y unit t ttre, att robust contreat(dose) contreatcontrols(untreatedZero) ///
    cfxlattice(dose0) cfxplottypes(psot) safety(off)
dis "VG7b_noparent_RC=" _rc

* ccnoabsorb under contreat() but with NO contreatcontrols(): sthlp:1121 says inert
use `C', clear
capture noisily wooldid y unit t ttre, att robust contreat(dose) safety(off)
dis "CC0_RC=" _rc
capture matrix CC0 = e(wdidmainresults)
use `C', clear
capture noisily wooldid y unit t ttre, att robust contreat(dose) ccnoabsorb safety(off)
dis "CC1_RC=" _rc
capture matrix CC1 = e(wdidmainresults)
capture noisily dis "CC_MRELDIF=" mreldif(CC1, CC0)

dis "DONE_F3B"
