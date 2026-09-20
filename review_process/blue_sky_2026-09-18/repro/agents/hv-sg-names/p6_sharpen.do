*==========================================================================
* p6_sharpen.do -- two sharpening checks for the name findings:
*  (1) ES_att vs ES_att_it on a design where cell weighting BITES, so the
*      "att_i has no ES_att_i_* row" issue has a wrong-number consequence.
*  (2) rigorous count of what coefplot actually plots when handed the
*      help's non-matching drop(*_s1*) pattern.
*  SELF-CONTAINED.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop dumpb
prog def dumpb
  args tag
  local nn : colnames e(b)
  foreach n of local nn {
    di "BVAL|`tag'|`n'|" %21x _b[`n']
  }
end

* ---- i = GROUP (coarser than unit), cohorts of very unequal size ----
clear
set obs 900
gen long unit = _n
gen byte grp = .
* group 1..3 : onset 2003, sizes 300 / 60 / 60
replace grp = 1 if _n <=300
replace grp = 2 if inrange(_n,301,360)
replace grp = 3 if inrange(_n,361,420)
* group 4..5 : onset 2006, sizes 30 / 30
replace grp = 4 if inrange(_n,421,450)
replace grp = 5 if inrange(_n,451,480)
* group 6..8 : never treated
replace grp = 6 if inrange(_n,481,620)
replace grp = 7 if inrange(_n,621,760)
replace grp = 8 if inrange(_n,761,900)
gen int ttre = .
replace ttre = 2003 if inlist(grp,1,2,3)
replace ttre = 2006 if inlist(grp,4,5)
gen double coheff = 0
replace coheff = 0.4 if ttre == 2003
replace coheff = 2.5 if ttre == 2006
gen double g_fe = .
forval g = 1/8 {
  qui replace g_fe = rnormal() if grp == `g' & mi(g_fe)
}
bysort grp (unit): replace g_fe = g_fe[1]
expand 10
bysort unit: gen int year = 1999 + _n
gen double post = (year >= ttre) & ttre < .
gen double y = g_fe + 0.2*(year-2000) + post*coheff + 0.25*rnormal()
save "p6_data.dta", replace

di "=== RUN U: att + att_i together, i = grp ==="
cap noi wooldid y grp year ttre, att att_i esprelength(1) espostlength(2) robust
di "U_RC=" _rc
if _rc == 0 dumpb "U"

di "=== RUN V: coefplot rigour on the help's drop(*_s1*) ==="
use "p6_data.dta", clear
gen byte sgA = 1 + mod(grp, 2)
cap noi wooldid y grp year ttre, att subgroups(sgA) esprelength(1) espostlength(1) robust
di "V_RC=" _rc
local nb = colsof(e(b))
di "V_NCOEF=" `nb'

preserve
clear
cap noi coefplot ., drop(*_s1*) nodraw generate(cpA_)
di "V_DROP_HELPNAME_RC=" _rc
cap qui count if !mi(cpA__b)
di "V_DROP_HELPNAME_NPLOTTED=" r(N)
restore

preserve
clear
cap noi coefplot ., drop(*SG1*) nodraw generate(cpB_)
di "V_DROP_TRUENAME_RC=" _rc
cap qui count if !mi(cpB__b)
di "V_DROP_TRUENAME_NPLOTTED=" r(N)
restore

preserve
clear
cap noi coefplot ., nodraw generate(cpC_)
di "V_NOOPT_RC=" _rc
cap qui count if !mi(cpC__b)
di "V_NOOPT_NPLOTTED=" r(N)
restore

di "PROBE6_DONE"
