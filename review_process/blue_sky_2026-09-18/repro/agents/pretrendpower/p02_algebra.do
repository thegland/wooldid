clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

local N 80
local T 10
set obs `N'
gen long unit = _n
gen double a_i = rnormal()
gen int onset = .
replace onset = 5 if mod(unit,5)==1
replace onset = 7 if mod(unit,5)==2
replace onset = 9 if mod(unit,5)==3
expand `T'
bysort unit: gen int t = _n
gen double f_t = 0.1*t
gen byte D = (onset<. & t>=onset)
gen double y = a_i + f_t + 1.0*D + rnormal()

timer clear
timer on 2
qui wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2)
timer off 2
dis "PLAIN_ONLY"
timer list 2

timer on 3
qui wooldid y unit t onset, esfixedbaseperiod esprelength(3) espostlength(2) pretrendpower(80)
timer off 3
dis "WITH_PTP"
timer list 3

* ---------- hand GLS on the posted b/V ----------
mat B = e(b_es_att)
mat V = e(V_es_att)
local dfr = e(df_r)
mat S = V[1..3,1..3]
mat S = (S + S')/2
mat bp = (B[1,1] \ B[1,2] \ B[1,3])
mat m  = (-3 \ -2 \ -1)
mat Si = invsym(S)
mat W1 = m'*Si*m
scalar c2 = W1[1,1]
mat W2 = m'*Si*bp
scalar hslope = W2[1,1]/c2
scalar hse = 1/sqrt(c2)
mat W3 = bp'*Si*bp
scalar hW = W3[1,1]
scalar hp = Ftail(3,`dfr',hW/3)
scalar lam80 = npnF(3,`dfr',invF(3,`dfr',.95),1-.80)
scalar hmdt80 = sqrt(lam80/c2)
scalar hbias80 = hmdt80*(e(roth_hbar_att) - (-1))

dis "HAND slope  = " %21x hslope   "   E: " %21x e(roth_slope_att)
dis "HAND se     = " %21x hse      "   E: " %21x e(roth_slope_se_att)
dis "HAND W      = " %21x hW       "   E: " %21x e(roth_W_att)
dis "HAND p      = " %21x hp       "   E: " %21x e(roth_prep_att)
mat RP = e(roth_power_att)
dis "HAND mdt80  = " %21x hmdt80   "   E: " %21x RP[1,2]
dis "HAND bias80 = " %21x hbias80  "   E: " %21x RP[1,4]
dis "HAND nis80  = ?                        E: " %21x RP[1,3]

* ---------- external oracle: Roth's own pretrends package ----------
tempname bb vv
mat `bb' = e(b_es_att)
mat `vv' = e(V_es_att)
mat colnames `bb' = c1 c2 c3 c4 c5
mat colnames `vv' = c1 c2 c3 c4 c5
mat rownames `vv' = c1 c2 c3 c4 c5
cap noi pretrends power 0.80, b(`bb') v(`vv') time(-4 -3 -2 0 1) ref(-1)
dis "PRET80_RC=" _rc
dis "PRETRENDS slope80 = " %21x r(slope)
cap noi pretrends power 0.50, b(`bb') v(`vv') time(-4 -3 -2 0 1) ref(-1)
dis "PRET50_RC=" _rc
dis "PRETRENDS slope50 = " %21x r(slope)

* what power does pretrends attribute to wooldid's own mdt_nis(80)?
scalar wnis80 = RP[1,3]
cap noi pretrends , b(`bb') v(`vv') time(-4 -3 -2 0 1) ref(-1) slope(`=wnis80') nocoefplot
dis "PT_AT_WNIS_RC=" _rc
dis "PRETRENDS power at wooldid mdt_nis(80) = " r(Power)
scalar wwald80 = RP[1,2]
cap noi pretrends , b(`bb') v(`vv') time(-4 -3 -2 0 1) ref(-1) slope(`=wwald80') nocoefplot
dis "PT_AT_WWALD_RC=" _rc
dis "PRETRENDS power at wooldid mdt_wald(80) = " r(Power)
