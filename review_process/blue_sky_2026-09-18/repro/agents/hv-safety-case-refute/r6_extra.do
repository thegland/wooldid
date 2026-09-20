* r6_extra.do -- SELF-CONTAINED.
* (a) Does the float-merged panel that lost the BS-34 caution under contreatwithin
*     actually report a different number from the panel the user meant?
* (b) Does contreatwithin + contreat() also swallow the esfixedbaseperiod r(459)
*     "cannot find sufficient pre- and post-treatment periods" refusal (ado:908)?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

which wooldid

******************** (a) float-merged identifiers ********************
clear
set obs 60
gen float id  = 20000000 + _n
gen double idd = 20000000 + _n
egen double g = group(id)
gen double ttre = .
replace ttre = 4 if mod(g,3) == 1
replace ttre = 6 if mod(g,3) == 2
gen double ue = rnormal()
expand 8
bys id (ttre): gen double t = _n
replace t = mod(_n-1,8)+1
gen double dose = 0.5 + runiform()
gen byte post = (t >= ttre) & ttre < .
gen double y = 2 + 0.4*t + ue + 1.1*post*dose + rnormal()*0.4
tempfile F
save "`F'"

use "`F'", clear
dis "LEGSTART_Afloat"
cap noi wooldid y id t ttre, att contreat(dose) contreatcontrols(byCohort) contreatwithin cluster(id)
dis "LEG_Afloat_RC=" _rc
cap noi dis "Afloat_N=" e(N)
cap noi mat li e(wdidmainresults), format(%21x)
dis "LEGEND_Afloat"

use "`F'", clear
dis "LEGSTART_Bdouble"
cap noi wooldid y idd t ttre, att contreat(dose) contreatcontrols(byCohort) contreatwithin cluster(idd)
dis "LEG_Bdouble_RC=" _rc
cap noi dis "Bdouble_N=" e(N)
cap noi mat li e(wdidmainresults), format(%21x)
dis "LEGEND_Bdouble"

******************** (b) no usable pre-period for any treated unit ********************
clear
set obs 40
gen double id = _n
gen double ttre = .
replace ttre = 1 if _n <= 30          // treated from the very first period: no pre-period
gen double ue = rnormal()
expand 8
bys id: gen double t = _n
gen double dose = 0.5 + runiform()
gen byte post = (t >= ttre) & ttre < .
gen double y = 2 + 0.4*t + ue + 1.1*post*dose + rnormal()*0.4
tempfile Q
save "`Q'"

use "`Q'", clear
dis "LEGSTART_noctw_esf"
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(byCohort)
dis "LEG_noctw_esf_RC=" _rc

use "`Q'", clear
dis "LEGSTART_ctw_esf"
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(2) espostlength(2) contreat(dose) contreatcontrols(byCohort) contreatwithin
dis "LEG_ctw_esf_RC=" _rc

use "`Q'", clear
dis "LEGSTART_noctw_pool"
cap noi wooldid y id t ttre, att esprelength(2) espostlength(2) contreat(dose) contreatcontrols(byCohort)
dis "LEG_noctw_pool_RC=" _rc

use "`Q'", clear
dis "LEGSTART_ctw_pool"
cap noi wooldid y id t ttre, att esprelength(2) espostlength(2) contreat(dose) contreatcontrols(byCohort) contreatwithin
dis "LEG_ctw_pool_RC=" _rc

dis "ALLDONE"
