* r1_core.do -- SKEPTIC's own repro of hv-cwm-01.
* Self-contained: generates its own panel, no project data, no project test suite.
*
* NOTE TO A VERIFIER: the `dumparm' estimand loop below is INERT.  e(b)'s column names
* carry an equation prefix (`main:att'), so colnumb(e(b),"att") returns missing and the
* `<TAG>_att = ...' lines never print.  Only the `<TAG>_gtrows' line from that program is
* live.  The four estimands per arm are read from wooldid's own displayed table in the log
* (ARM A at log line ~175, ARM B ~320, ARM S ~460, ARM D ~598), and every claim in
* FINDINGS.md that depends on them is cross-checked against the `<TAG>_wtdmean_attit' /
* `<TAG>_simplemean' lines, which ARE live and are computed from e(gteffects).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771102
set sortseed 771102

which wooldid

*----------------------------------------------------------------------
* Panel: 14 counties x 10 years x 4 establishments per county-year.
* Staggered adoption, 6 never-treated counties.
* sw  = establishment employment (sample weight)
* pop = county population, ~40:1 spread, CORRELATED with the county's effect
* popg= same thing growing over time (a time-varying multiplier)
*----------------------------------------------------------------------
clear
set obs 14
gen int county = _n
gen int onset = .
replace onset = 2011 if inrange(county,1,3)
replace onset = 2013 if inrange(county,4,6)
replace onset = 2015 if inrange(county,7,8)
* counties 9-14 never treated
gen double pop = 1
replace pop = 40 if inlist(county,1,5)
replace pop = 18 if inlist(county,2,7)
replace pop = 6  if inlist(county,4,10,12)
replace pop = 3  if inlist(county,3,9)
gen double cfe = rnormal(0,1)
* heterogeneous effect correlated with pop -- so weighting genuinely matters
gen double teff = 1.2 + 0.055*pop + rnormal(0,0.3)
expand 10
bysort county: gen int year = 2007 + _n
expand 4
bysort county year: gen int estab = _n
gen double sw = exp(rnormal(1.1,0.55))
gen double yfe = 0.12*(year-2007)
gen byte post = (onset<. & year>=onset)
gen double y = 3 + cfe + yfe + teff*post + rnormal(0,0.8)
gen double popg = pop*(1.15^(year-2007))
compress
save "panel_r1.dta", replace

*----------------------------------------------------------------------
* helper: dump the four estimands + the gteffects table for one arm
*----------------------------------------------------------------------
capture program drop dumparm
program define dumparm
  syntax , TAG(string)
  foreach e in att att_it att_i att_itime {
    local c = colnumb(e(b),"`e'")
    if `c' < . {
      local v = el(e(b),1,`c')
      dis "`tag'_`e' = " %21x `v' "   (" %18.10f `v' ")"
    }
  }
  matrix GT_`tag' = e(gteffects)
  dis "`tag'_gtrows = " rowsof(GT_`tag')
end

*======================================================================
* ARM A -- no cwm
*======================================================================
use "panel_r1.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it att_i att_itime gteffects
dis "ARMA_RC=" _rc
dumparm , tag(A)
gen byte esmpA = e(sample)
qui count if esmpA
dis "ARMA_NESAMP=" r(N) "  ARMA_NROWS=" _N

*======================================================================
* ARM B -- cwm(pop), pop constant within county
*======================================================================
use "panel_r1.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it att_i att_itime gteffects customweightmultiplier(pop)
dis "ARMB_RC=" _rc
dumparm , tag(B)

*======================================================================
* ARM S -- cwm(three), a CONSTANT multiplier. Must be an exact no-op.
*======================================================================
use "panel_r1.dta", clear
gen double three = 3
cap noi wooldid y county year onset [aw=sw], att att_it att_i att_itime gteffects customweightmultiplier(three)
dis "ARMS_RC=" _rc
dumparm , tag(S)

*======================================================================
* ARM D -- cwm(popg), TIME-VARYING multiplier
*======================================================================
use "panel_r1.dta", clear
cap noi wooldid y county year onset [aw=sw], att att_it att_i att_itime gteffects customweightmultiplier(popg)
dis "ARMD_RC=" _rc
dumparm , tag(D)

*======================================================================
* Build the RAW-DATA ORACLE for the cell weights, and the two candidate
* att_it numbers (documented simple mean vs. reported).
*======================================================================
use "panel_r1.dta", clear
gegen double swsum  = total(sw),     by(county year)
gegen double swpsum = total(sw*pop), by(county year)
gegen double swgsum = total(sw*popg),by(county year)
gen double oracle_pop  = swpsum/swsum
gen double oracle_popg = swgsum/swsum
gegen byte tagc = tag(county year)
keep if tagc
keep county year oracle_pop oracle_popg
rename county i
rename year caltime
save "oracle_r1.dta", replace

foreach tag in A B S D {
  clear
  svmat double GT_`tag', names(col)
  gen int _rowid = _n
  merge 1:1 i caltime using "oracle_r1.dta", keep(master match) nogen
  * (1) reported aggregate reproduced by the w_att_it column (sthlp:1136 identity)
  qui sum estimate [aw=w_att_it] if relyr >= 0
  dis "`tag'_wtdmean_attit = " %21x r(mean) "  (" %18.10f r(mean) ")"
  * (2) the DOCUMENTED estimand: the simple, unweighted mean of the cells
  qui sum estimate if relyr >= 0
  dis "`tag'_simplemean    = " %21x r(mean) "  (" %18.10f r(mean) ")  ncells=" r(N)
  * (3) oracle check on the cell weights
  qui sum w_att_it
  dis "`tag'_w_att_it min=" %12.8f r(min) " max=" %12.8f r(max)
  if "`tag'" == "B" {
    gen double dev = reldif(w_att_it, oracle_pop)
    qui sum dev
    dis "B_ORACLE_maxreldif_cellw = " %12.4e r(max)
  }
  if "`tag'" == "D" {
    gen double dev = reldif(w_att_it, oracle_popg)
    qui sum dev
    dis "D_ORACLE_maxreldif_cellw = " %12.4e r(max)
  }
  * (4) per-cohort w_att_i totals, POST block and PRE block separately
  qui levelsof onset if onset < ., local(ons)
  foreach o of local ons {
    qui sum w_att_i if onset == `o' & relyr >= 0
    local tp = r(sum)
    qui sum w_att_i if onset == `o' & relyr < 0
    local tq = r(sum)
    dis "`tag'_watti onset=`o'  POSTtot=" %12.8f `tp' "  PREtot=" %12.8f `tq'
  }
}

*======================================================================
* SECONDARY CLAIM CONTROL: is "pre row and main row weight the cohorts
* differently" cwm-specific?  Arm T: NO cwm at all, but a TIME-VARYING
* sample weight, att only.  If the pre/main cohort shares differ here too,
* the phenomenon is generic to time-varying weights, not a cwm defect.
*======================================================================
use "panel_r1.dta", clear
gen double swt = sw*(1.15^(year-2007))
cap noi wooldid y county year onset [aw=swt], att gteffects
dis "ARMT_RC=" _rc
matrix GT_T = e(gteffects)
clear
svmat double GT_T, names(col)
qui sum w_att if relyr >= 0
local Tpost = r(sum)
qui sum w_att if relyr < 0
local Tpre = r(sum)
qui levelsof onset if onset < ., local(ons)
foreach o of local ons {
  qui sum w_att if onset == `o' & relyr >= 0
  local sp = r(sum)/`Tpost'
  qui sum w_att if onset == `o' & relyr < 0
  local sq = r(sum)/`Tpre'
  dis "T_share onset=`o'  POSTshare=" %10.6f `sp' "  PREshare=" %10.6f `sq' "  reldif=" %10.6f reldif(`sp',`sq')
}

dis "R1_DONE"
