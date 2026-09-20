* r11_boottstat.do -- SELF-CONTAINED.
* (1) Is the printed "joint pre-trend test" of pretrendpower() reproduced EXACTLY by a
*     hand `test` on e(V) after a bootstrap fit?  sthlp:532 says of the jointtests family
*     "Neither is the Wald test on e(V) that a hand test reproduces; that one is a
*     normal-theory test on the bootstrap covariance and is not printed."
* (2) Under boottstat, where the run's own caution says "the joint pre-trend tests ...
*     are bootstrapped as test statistics", is pretrendpower()'s printed joint pre-trend
*     test bootstrapped, or is it still the normal-theory Wald on F(k, G-1)?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 30
gen state = _n
gen onset = .
replace onset = 7 + mod(_n, 6) if mod(_n,3)==0
expand 2
bysort state: gen cty = _n
gen county = state*10 + cty
expand 16
bysort county: gen t = _n
bysort county (t): gen double ushock = rnormal(0,.3) if _n==1
bysort county (t): replace ushock = ushock[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .8 + .12*rel, 0)
gen double y = 2 + .05*t + ushock + te + e
drop rel te

local OPT "cluster(state) esfixedbaseperiod esprelength(6) espostlength(4) jointtests pretrendpower(80)"
local PRE "ES_att_pre7 ES_att_pre6 ES_att_pre5 ES_att_pre4 ES_att_pre3 ES_att_pre2"

dis "==== analytic"
qui cap noi wooldid y county t onset, `OPT'
cap noi test `PRE'
dis "HAND|ANALYTIC|handp=" %14.10g r(p) "|handdf=" r(df) "|drop=" r(drop) "|roth_prep=" %14.10g e(roth_prep_att) "|rankjoint=" %14.10g e(joint_horiz) "|dfr=" e(df_r)

dis "==== bootstrap, default display"
qui cap noi wooldid y county t onset, `OPT' bootstrap(standard) reps(60) rseed(2468)
local rp = e(roth_prep_att)
local rj = e(joint_horiz)
local df = e(df_r)
local jf "`e(boot_jt_family)'"
cap noi test `PRE'
dis "HAND|BOOT60|handp=" %14.10g r(p) "|handdf=" r(df) "|drop=" r(drop) "|roth_prep=" %14.10g `rp' "|rankjoint=" %14.10g `rj' "|dfr=`df'|jt_family=`jf'"

dis "==== bootstrap + boottstat, same seed and reps"
qui cap noi wooldid y county t onset, `OPT' bootstrap(standard) reps(60) rseed(2468) boottstat
local rp2 = e(roth_prep_att)
local rj2 = e(joint_horiz)
local df2 = e(df_r)
local jf2 "`e(boot_jt_family)'"
local sl2 = e(roth_slope_att)
local se2 = e(roth_slope_se_att)
cap noi test `PRE'
dis "HAND|BOOTTSTAT60|handp=" %14.10g r(p) "|handdf=" r(df) "|drop=" r(drop) "|roth_prep=" %14.10g `rp2' "|rankjoint=" %14.10g `rj2' "|dfr=`df2'|jt_family=`jf2'|slope=" %14.10g `sl2' "|se=" %14.10g `se2'
dis "DONE"
