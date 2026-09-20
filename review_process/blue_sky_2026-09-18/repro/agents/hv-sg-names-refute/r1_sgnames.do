*=======================================================================
* r1_sgnames.do -- SKEPTIC's independent repro of hv-sgname-01.
*
* CLAIM UNDER TEST: sthlp:1307 says "Under subgroups() each subgroup's
* rows take the suffix _s#" for e(b)/e(V), but the code actually posts
* an INFIX _SG<value>_ after the family stem (main_SG1_att), and the
* string "_SG" appears nowhere in wooldid.sthlp.
*
* This do-file is SELF-CONTAINED: it generates its own panel, uses no
* repo data and no repo test-suite file, and prints tagged sentinels.
* Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do r1_sgnames.do
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

*--- HAZARD (b): assert we are on the pinned build, not the PLUS v1.3 ---
which wooldid

*--- helper: dump a matrix's column names one per line (80-col log safe) --
cap prog drop dumpcols
prog def dumpcols
  args mname tag
  local cn : colnames `mname'
  local ce : coleq   `mname'
  dis "DUMPCOLS_`tag'_NCOL=" colsof(`mname')
  forval j = 1/`=colsof(`mname')' {
    local nj : word `j' of `cn'
    local ej : word `j' of `ce'
    dis "DC_`tag' `j' eq=<`ej'> nm=<`nj'>"
  }
end
cap prog drop dumprows
prog def dumprows
  args mname tag
  local rn : rownames `mname'
  local re : roweq   `mname'
  dis "DUMPROWS_`tag'_NROW=" rowsof(`mname')
  forval j = 1/`=rowsof(`mname')' {
    local nj : word `j' of `rn'
    local ej : word `j' of `re'
    dis "DR_`tag' `j' eq=<`ej'> nm=<`nj'>"
  }
end

*=================== OWN DATA: 90 units x 10 years =====================
* 5 unit classes by mod(unit,5): classes 1,2,3 are treated cohorts with
* onsets 2004/2006/2008; classes 0,4 are never-treated (36 units).
* Subgroup sg in {1,2} alternates by unit, so each cohort carries both.
clear
set obs 90
gen long unit  = _n
gen byte sg    = 1 + mod(_n, 2)
gen byte cls   = mod(_n, 5)
gen int  ttre  = .
replace  ttre  = 2004 if cls == 1
replace  ttre  = 2006 if cls == 2
replace  ttre  = 2008 if cls == 3
gen double ai  = rnormal()
expand 10
bysort unit: gen int year = 2000 + _n - 1
gen byte post  = (ttre < .) & (year >= ttre)
gen double tau = 1.0 + 0.5*(sg == 2)
gen double y   = ai + 0.10*(year-2000) + tau*post + 0.30*rnormal()
compress
dis "DATA_N=" _N

*=================== RUN A: four estimands x subgroups =================
cap noi wooldid y unit year ttre, att att_it att_i att_itime ///
    subgroups(sg) esprelength(2) espostlength(2) robust
dis "RUNA_RC=" _rc

tempname B V
matrix `B' = e(b)
matrix `V' = e(V)
dumpcols `B' EB
dumprows `V' EV

*--- Q1: is there ANY equation structure that could rescue "_s#"? -------
* (a flat colname list with eq=<_> means [main_s1]att etc. cannot work)
local ceq : coleq `B'
dis "EB_COLEQ_RAW=<`ceq'>"

*=================== Q2: the help's names vs the real names ============
* sthlp:1303-1307 declares e(b)/e(V) colnames to be what esttab,
* coefplot, lincom and test address a coefficient by.
foreach nm in main_att_s1 main_att_s2 main_SG1_att main_SG2_att ///
              main_att_it_s1 main_SG1_att_it main_att_itime_s1 ///
              main_SG1_att_itime pretreatment_att_s1 pretreatment_SG1_att ///
              ES_att_post0_s1 ES_att_s1_post0 ES_att_SG1_post0 {
    cap di _b[`nm']
    dis "BREF `nm' RC=" _rc
}

cap noi test main_att_s1 = main_att_s2
dis "TEST_HELPNAME_RC=" _rc
cap noi test main_SG1_att = main_SG2_att
dis "TEST_TRUENAME_RC=" _rc
cap noi lincom main_att_s1
dis "LINCOM_HELPNAME_RC=" _rc
cap noi lincom main_SG1_att
dis "LINCOM_TRUENAME_RC=" _rc

*--- estout / esttab, the other consumer the help names ----------------
cap noi esttab ., keep(*_s1*) nostar
dis "ESTTAB_KEEP_HELPNAME_RC=" _rc
cap noi esttab ., keep(*SG1*) nostar
dis "ESTTAB_KEEP_TRUENAME_RC=" _rc
cap noi esttab ., drop(*_s1*) nostar
dis "ESTTAB_DROP_HELPNAME_RC=" _rc

*--- coefplot: is drop(*_s1*) a SILENT no-op? --------------------------
* count plotted coefficients via coefplot's generate() output
cap drop cpA*
cap noi coefplot, generate(cpA) nodraw
dis "CP_PLAIN_RC=" _rc
cap count if cpA_b < .
dis "CP_PLAIN_NCOEF=" r(N)
cap drop cpB*
cap noi coefplot, drop(*_s1*) generate(cpB) nodraw
dis "CP_DROPHELP_RC=" _rc
cap count if cpB_b < .
dis "CP_DROPHELP_NCOEF=" r(N)
cap drop cpC*
cap noi coefplot, drop(*SG1*) generate(cpC) nodraw
dis "CP_DROPTRUE_RC=" _rc
cap count if cpC_b < .
dis "CP_DROPTRUE_NCOEF=" r(N)
cap drop cpD*
cap noi coefplot, keep(*_s1*) generate(cpD) nodraw
dis "CP_KEEPHELP_RC=" _rc
cap count if cpD_b < .
dis "CP_KEEPHELP_NCOEF=" r(N)

*=================== Q3: the OTHER two subgroup key regimes ============
foreach en in b_es_att_s1 b_es_att_SG1 b_es_att_it_s1 V_es_att_s1 ///
              escomp_s1 wdidsgresults_1 wdidsgresults_s1 wdidsgresults_SG1 {
    cap confirm matrix e(`en')
    dis "EMAT `en' RC=" _rc
}

*=================== Q4: does the SG ES covariance SPLICE land? ========
* ado:10910/11061 say the SG event-study block is spliced into e(V)
* under the ES_<base>_SG<sg>_ prefix. Verify the numbers actually match
* the standalone e(V_es_att_s1).
tempname Ves
cap matrix `Ves' = e(V_es_att_s1)
dis "VES_S1_RC=" _rc
cap dumprows `Ves' VESS1
cap dumpcols `Ves' VESS1
foreach pr in pre2 pre1 post0 post1 {
  foreach pc in pre2 pre1 post0 post1 {
    tempname a b
    cap matrix `a' = `Ves'["`pr'","`pc'"]
    local ra = _rc
    cap matrix `b' = `V'["ES_att_SG1_`pr'","ES_att_SG1_`pc'"]
    local rb = _rc
    if `ra' == 0 & `rb' == 0 {
      dis "SPLICE `pr' `pc' ves=" %21x `a'[1,1] " eV=" %21x `b'[1,1] ///
          " same=" (`a'[1,1] == `b'[1,1])
    }
    else dis "SPLICE `pr' `pc' RC_ves=`ra' RC_eV=`rb'"
  }
}

*=================== Q5: what the PRINTED table calls the rows =========
dis "--- e(wdidsgresults_1) rownames ---"
tempname SG1
cap matrix `SG1' = e(wdidsgresults_1)
cap dumprows `SG1' SGRES1

*=================== Q6: adversarial scalar name (sthlp:1398) ==========
cap noi wooldid y unit year ttre, att subgroups(sg) esprelength(2) ///
    espostlength(2) robust esfixedbaseperiod adversarial
dis "RUNB_RC=" _rc
dis "--- ereturn list, all (grep adversarial) ---"
cap noi ereturn list, all

dis "R1_DONE"
