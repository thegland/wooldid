* ==========================================================================
* r04_realistic.do -- how big is the fpsg-01 distortion on data a real applied
* economist would actually have?  Noisy, staggered, unbalanced, with two
* realistic item-nonresponse patterns for the subgroup variable.
*
* ARM 1  clean sg (benchmark)
* ARM 2  MCAR: 15% of post rows have sg unobserved
* ARM 3  AS BUILT this arm is NOT effect-correlated: the rank below is computed
*        AFTER expand, so it ranks ROWS, and `effrank > 20` is true for almost
*        every unit -- the realized arm is "the subgroup code stops being
*        collected from relative year 2 on for nearly all units" (mean per-unit
*        weight mass 0.292).  Still a realistic pattern, but the correctly
*        built effect-correlated arm is in r05_contrast.do, which also has a
*        known truth.  Left as-is because the arm is informative as realized.
*
* ORACLE.  For each arm, e(gteffects) is pulled down and the SAME cell
* estimates are aggregated two ways over the cells that are IN subgroup 1:
*   EQUALUNIT  = simple mean over units of (that unit's w_att_it-weighted mean)
*   SHAREWTD   = sum(w_att_i * est) / sum(w_att_i)
* If wooldid's reported att_i for sg==1 equals SHAREWTD and not EQUALUNIT, the
* claim holds on realistic data.  The gap is then reported in SE units.
*
* NOTE: e(histogramestimates) is NOT a valid oracle here -- the margins call
* that builds it (ado:12405) carries no subgroup term in its subpop, so it
* averages each unit over ALL its post cells, in-subgroup or not.
* ==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---------------------------- panel ---------------------------------------
clear
set obs 90
gen long id = _n
gen byte sg = .
replace sg = 1 if id <= 30
replace sg = 0 if id > 30 & id <= 60
replace sg = 0 if id > 60
gen double ttre = .
replace ttre = 5 + 2*mod(id-1,3) if id <= 60      /* onsets 5, 7, 9 */
gen double eff = 0
replace eff = 1 + 2*runiform() if id <= 60        /* heterogeneous effects */
gen double u_i = rnormal(0,1)
expand 14
bysort id: gen int t = _n
* unbalanced entry/exit, as administrative panels are
gen byte keeprow = 1
replace keeprow = 0 if mod(id,7) == 0 & t <= 2
replace keeprow = 0 if mod(id,11) == 0 & t >= 13
drop if keeprow == 0
gen double y = u_i + 0.3*t + eff*(t >= ttre & ttre < .) + rnormal(0,0.6)

* ---- ARM 2: MCAR 15% of post rows lose sg
gen byte sg_mcar = sg
replace sg_mcar = . if ttre < . & t >= ttre & runiform() < 0.15

* ---- ARM 3: the biggest-effect units lose sg in their later years
* NB: computed after expand, so this ranks ROWS not units -- see header.
egen double effrank = rank(eff) if id <= 30
gen byte sg_corr = sg
replace sg_corr = . if id <= 30 & effrank > 20 & ttre < . & t >= ttre + 2

tempfile cellsg
preserve
  keep id t sg sg_mcar sg_corr
  duplicates drop
  rename t caltime
  save `cellsg', replace
restore
tempfile panel
save `panel', replace

* ---------------------------- the three arms -------------------------------
foreach arm in sg sg_mcar sg_corr {
  use `panel', clear
  dis "############### ARM `arm' ###############"
  dis "CAUTIONSCAN_`arm'_START"
  capture noisily wooldid y id t ttre, att att_i att_it att_itime subgroups(`arm') gteffects
  dis "ARM_`arm'_RC=" _rc
  dis "CAUTIONSCAN_`arm'_END"
  matrix R = e(wdidsgresults_1)
  local rn : rowfullnames R
  dis "ARM_`arm'_ROWNAMES: `rn'"
  scalar rep_att   = R[1,1]
  scalar rep_attit = R[2,1]
  scalar rep_atti  = R[3,1]
  scalar rep_attim = R[4,1]
  scalar se_atti   = R[3,2]
  matrix G = e(gteffects)

  clear
  svmat double G, names(col)
  * gteffects columns: i onset caltime relyr estimate se p lb ub inwindow nobs
  *                    w_att w_att_it w_att_i w_att_itime
  rename i id
  merge m:1 id caltime using `cellsg', keep(master match) nogen
  keep if relyr >= 0
  keep if `arm' == 1
  * SHARE-weighted (what the claim says wooldid computes)
  qui summ estimate [aw = w_att_i], meanonly
  scalar SHAREWTD = r(mean)
  * EQUAL-unit (what att_i is documented to be)
  preserve
    collapse (mean) ubar = estimate [aw = w_att_it], by(id)
    qui summ ubar, meanonly
    scalar EQUALUNIT = r(mean)
    qui count
    scalar NUNITS = r(N)
  restore
  * per-unit total w_att_i mass -- should be 1.0 for a fully-observed unit
  preserve
    collapse (sum) wmass = w_att_i, by(id)
    qui summ wmass
    scalar WMIN = r(min)
    scalar WMAX = r(max)
    scalar WMEAN = r(mean)
    scalar WTOT = r(sum)
  restore

  dis "---- ARM `arm' ----"
  dis "REPORTED att_i (sg==1)  = " %18.10f rep_atti
  dis "REPORTED att_itime      = " %18.10f rep_attim
  dis "REPORTED att            = " %18.10f rep_att
  dis "ORACLE   SHAREWTD       = " %18.10f SHAREWTD
  dis "ORACLE   EQUALUNIT      = " %18.10f EQUALUNIT
  dis "N units in sg1          = " NUNITS
  dis "per-unit w_att_i mass: min=" %9.6f WMIN " max=" %9.6f WMAX " mean=" %9.6f WMEAN " total=" %9.4f WTOT
  dis "GAP_reported_vs_SHARE   = " %18.10e (rep_atti - SHAREWTD)
  dis "GAP_reported_vs_EQUAL   = " %18.10e (rep_atti - EQUALUNIT)
  dis "SE of att_i             = " %18.10f se_atti
  dis "GAP_EQUAL_in_SE_units   = " %18.6f ((rep_atti - EQUALUNIT)/se_atti)
  dis "RELGAP_pct              = " %18.6f (100*(rep_atti - EQUALUNIT)/EQUALUNIT)
}

dis "ALLDONE"
