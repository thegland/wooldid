*=======================================================================
* r6_varabbrev.do -- sharpener for sk-prefix-01: is the silent
* redirection of _b[main_att] -> main_att_i governed by Stata's
* `set varabbrev' switch?  If yes, the finding can name a user-side
* remedy; if no, there is none.  SELF-CONTAINED.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 192
gen long unit = _n
gen int grp = .
replace grp = ceil(_n/40)            if _n <= 120
replace grp = 3 + ceil((_n-120)/6)   if _n > 120 & _n <= 132
replace grp = 5 + ceil((_n-132)/30)  if _n > 132
gen int ttre = .
replace ttre = 2004 if inrange(grp,1,3)
replace ttre = 2007 if inrange(grp,4,5)
gen double truetau = cond(inrange(grp,1,3), 0.4, cond(inrange(grp,4,5), 2.5, 0))
gen double ai = rnormal()
expand 10
bysort unit: gen int year = 2000 + _n - 1
gen byte post = (ttre < .) & (year >= ttre)
gen double y = ai + 0.05*(year-2000) + truetau*post + 0.20*rnormal()
compress

* prove no dataset variable is named main_att / main_att_i
cap confirm variable main_att
noi dis "VAR_main_att_RC=" _rc
cap confirm variable main_att_i
noi dis "VAR_main_att_i_RC=" _rc

qui wooldid y grp year ttre, att esprelength(2) espostlength(2) robust
noi dis "TRUE_ATT=" _b[main_att]

qui wooldid y grp year ttre, att_i esprelength(2) espostlength(2) robust
noi dis "ATT_I=" _b[main_att_i]

set varabbrev on
cap noi dis "VA_ON_main_att = " _b[main_att]
noi dis "VA_ON_RC=" _rc

set varabbrev off
cap noi dis "VA_OFF_main_att = " _b[main_att]
noi dis "VA_OFF_RC=" _rc
cap noi lincom main_att
noi dis "VA_OFF_LINCOM_RC=" _rc
cap noi dis "VA_OFF_main_att_i = " _b[main_att_i]
noi dis "VA_OFF_exact_RC=" _rc
set varabbrev on
dis "R6_DONE"
