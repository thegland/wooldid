* ===========================================================================
* repro_cwm_estimand.do   -- SELF-CONTAINED repro, no project data used.
*
* CLAIM: customweightmultiplier() silently stops att_it / att_i / att_itime
* being the estimand their help entries describe ("the simple, unweighted mean
* of all of the cohort-period-specific estimates", sthlp:220-221; "the simple
* average of cohort-specific treatment effects", sthlp:229).  cwm sits in the
* NUMERATOR of all four estimand weights (ado:1581,1630,1637,1644,1650) and in
* the DENOMINATOR of none (W_it/W_i are built from the sample weight alone,
* ado:1619-1620), so each treatment cell carries total weight
*        sum_cell(w*cwm) / sum_cell(w)                 (= cwm_cell if cell-constant)
* instead of 1.
*
* SECOND HALF: the bootctl cell-weighted gate (ado:6430-6442) asserts that "a
* within-cluster weight multiplier cancels out of them".  That is a DIFFERENT
* path - bootstrap(bayes) passes its Dirichlet draw as the SAMPLE WEIGHT
* (repw = basew*drawv, ado:9449-9451), which enters W_it too and therefore does
* cancel.  Arm C below shows the cancellation.  Both statements are true at once.
*
* Run:  stata-mp -e do repro_cwm_estimand.do   (then read repro_cwm_estimand.log)
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242
which wooldid

* ------------------------------------------------------------------ data ---
* 8 counties (4 adopt a policy in staggered years, 4 never), 10 years,
* 5 establishments per county-year carrying establishment sample weights.
* pop is the county's population: 1 for four counties, 50 for two of the
* treated ones - a perfectly ordinary 50:1 population spread.
clear
set obs 8
gen int county = _n
gen int adopt = .
replace adopt = 4 if county==1
replace adopt = 5 if county==2
replace adopt = 6 if county==3
replace adopt = 7 if county==4
gen double pop = cond(inlist(county,3,4), 50, 1)
expand 10
bys county: gen int year = _n
expand 5
bys county year: gen byte estab = _n
gen double sw = 0.5 + runiform()                     // establishment weight
gen byte post = (adopt < .) & (year >= adopt)
* effects are heterogeneous across counties: 1, 2, 3, 4 (+0.2 per year since adoption)
gen double tau = cond(county==1,1,cond(county==2,2,cond(county==3,3,cond(county==4,4,0))))
replace tau = tau + 0.2*(year-adopt) if post==1
replace tau = 0 if post==0
gen double y = 0.3*county + 0.1*year + tau*post + 0.05*rnormal()
gen double sw_x_pop = sw*pop                          // arm C
gen double popgrow  = pop*(1.15^year)                 // arm D: population that grows
compress
tempfile PANEL
save `PANEL'

* ---------------------------------------------------------------- ARM A ----
* no customweightmultiplier: the documented estimands
dis _n "=================== ARM A: no cwm ==================="
use `PANEL', clear
wooldid y county year adopt [aw=sw], att att_it att_i att_itime gteffects esprelength(3)
mat MA = e(wdidmainresults)
mat GA = e(gteffects)

* ---------------------------------------------------------------- ARM B ----
* the ordinary request: "population-weight my treatment effects"
dis _n "=================== ARM B: customweightmultiplier(pop) ==================="
use `PANEL', clear
wooldid y county year adopt [aw=sw], att att_it att_i att_itime gteffects esprelength(3) customweightmultiplier(pop)
mat MB = e(wdidmainresults)
mat GB = e(gteffects)

* ---------------------------------------------------------------- ARM C ----
* the SAME multiplier folded into the sample weight instead (what bootstrap(bayes)
* does with its Dirichlet draw): the cell weights are unchanged at 1
dis _n "=================== ARM C: [aw=sw*pop], no cwm ==================="
use `PANEL', clear
wooldid y county year adopt [aw=sw_x_pop], att att_it att_i att_itime gteffects esprelength(3)
mat GC = e(gteffects)

* ------------------------------------------------------------- verdicts ----
mata:
  function grab(string scalar m, string scalar cname) {
     return(st_matrix(m)[., selectindex(st_matrixcolstripe(m)[.,2]:==cname)[1]])
  }
  function mainb(string scalar m, string scalar nm) {
     real matrix M ; string matrix rs ; real colvector ii
     M = st_matrix(m) ; rs = st_matrixrowstripe(m)
     ii = selectindex((rs[.,1]:=="main"):*(rs[.,2]:==nm))
     if (length(ii)==0) return(.)
     return(M[ii[1], selectindex(st_matrixcolstripe(m)[.,2]:=="estimate")[1]])
  }
  rA = grab("GA","relyr") ; eA = grab("GA","estimate") ; wA = grab("GA","w_att_it")
  rB = grab("GB","relyr") ; eB = grab("GB","estimate") ; wB = grab("GB","w_att_it")
  rC = grab("GC","relyr") ;                              wC = grab("GC","w_att_it")
  pA = selectindex(rA:>=0) ; pB = selectindex(rB:>=0) ; pC = selectindex(rC:>=0)
  printf("\n--------------------------- VERDICT ---------------------------\n")
  printf("ARM A  att_it reported              = %12.8f\n", mainb("MA","att_it"))
  printf("ARM A  simple mean of gteffects     = %12.8f   <- they agree\n", mean(eA[pA]))
  printf("ARM A  cell weight w_att_it min/max = %g / %g\n", min(wA[pA]), max(wA[pA]))
  printf("\n")
  printf("ARM B  att_it reported              = %12.8f\n", mainb("MB","att_it"))
  printf("ARM B  simple mean of gteffects     = %12.8f   <- DOCUMENTED estimand\n", mean(eB[pB]))
  printf("ARM B  pop-weighted mean of cells   = %12.8f   <- what is reported\n", sum(eB[pB]:*wB[pB])/sum(wB[pB]))
  printf("ARM B  cell weight w_att_it min/max = %g / %g   <- should be 1/1\n", min(wB[pB]), max(wB[pB]))
  printf("ARM B  relative error vs documented = %g\n", abs(mainb("MB","att_it")-mean(eB[pB]))/abs(mean(eB[pB])))
  printf("ARM B  att reported                 = %12.8f   <- att_it has collapsed onto att\n", mainb("MB","att"))
  printf("\n")
  printf("ARM C  cell weight w_att_it min/max = %g / %g   <- multiplier CANCELS\n", min(wC[pC]), max(wC[pC]))
end

* ----- per-cohort weight totals, pre block vs post block, under cwm ---------
dis _n "=== ARM B: per-cohort TOTAL estimand weight, pre rows vs post rows ==="
dis "(with no cwm every one of these is exactly 1: arm A table below)"
foreach g in B A {
  preserve
    clear
    svmat double G`g', names(col)
    gen byte pst = relyr >= 0
    collapse (sum) w_att_i w_att_itime w_att_it (count) ncells=relyr, by(i pst)
    dis as result "---- ARM `g' ----"
    list, noobs sepby(i)
  restore
}
* --------------------------------------------------------------- ARM D ----
* a cwm that VARIES OVER TIME inside the cohort: W_i is total(w) by
* (posttreatment, i) (ado:1620), so the pretreatment att_i row and the main
* att_i row end up weighting the cohorts differently.
dis _n "=================== ARM D: customweightmultiplier(popgrow) ==================="
use `PANEL', clear
wooldid y county year adopt [aw=sw], att_i att_itime gteffects esprelength(3) customweightmultiplier(popgrow)
mat GD = e(gteffects)
dis _n "=== ARM D: per-cohort TOTAL estimand weight, pre block (pst=0) vs post block (pst=1) ==="
dis "(arm A has 1 in every cell of this table; here the pre and post blocks differ)"
preserve
  clear
  svmat double GD, names(col)
  gen byte pst = relyr >= 0
  collapse (sum) w_att_i w_att_itime (count) ncells=relyr, by(i pst)
  list, noobs sepby(i)
  reshape wide w_att_i w_att_itime ncells, i(i) j(pst)
  gen double ratio_prepost = w_att_i0/w_att_i1
  dis as result "per-cohort pre/post weight RATIO (all 1 with no cwm, and equal across cohorts iff comparable)"
  list i w_att_i0 w_att_i1 ratio_prepost, noobs
restore

dis "DONE_REPRO"
