* ===========================================================================
* p10 -- SELF-CONTAINED REPRO: a subgroup variable that is MISSING on some of a
*        unit's rows silently turns subgroup att_i / att_itime from the
*        equal-unit average they are defined to be into a SHARE-weighted average.
*
* wooldid REFUSES (r459, default safety) a unit whose subgroup value CHANGES
* over time, and its own error text gives the reason:
*     "att_i and att_itime enter a split unit into each subgroup weighted by its
*      share of sample weight there, rather than with the equal unit weighting
*      those estimands normally imply"
* A unit whose subgroup value is merely MISSING on some rows produces exactly
* the same share weighting -- and is accepted with no error and no caution.
*
* The DGP is noiseless, so every estimate is exact and the two candidate
* definitions are far apart and unmistakable:
*      unit 1 has effect 5.0, units 2-5 have effect 1.0, all in subgroup 1
*      equal-unit att_i(sg==1)  = (5 + 1 + 1 + 1 + 1)/5            = 1.80
*      share-weighted att_i     = (s*5 + 4*1)/(s + 4), s = 1/6     = 1.16
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---------------- a small noiseless panel ----------------
qui set obs 60
gen long unit = _n
gen double ttre = .
qui replace ttre = 2004 if unit <= 30          // 30 treated units, common onset
gen double a_i = 0.1*unit
expand 10
bysort unit: gen int t = 2000 + _n - 1
gen double g_t = 0.3*t
gen byte post = (ttre < . & t >= ttre)

* unit-level treatment effect: unit 1 is an outlier
gen double tau_u = 1.0
qui replace tau_u = 5.0 if unit == 1
gen double truetau = post * tau_u
gen double y = a_i + g_t + truetau

* subgroup: units 1-5 in subgroup 1, units 6-30 in subgroup 0
gen byte sg = .
qui replace sg = 1 if unit <= 5
qui replace sg = 0 if unit > 5

* ---- the ONLY manipulation: blank unit 1's subgroup in 5 of its 6 post years
gen byte sg_miss = sg
qui replace sg_miss = . if unit == 1 & t >= 2005

tab sg_miss if post==1 & unit<=5, m
tab unit if post==1 & sg_miss==1 & unit<=5

* ---------------- run A: clean subgroup variable ----------------
dis "@@@@@@@@@@ A: clean sg @@@@@@@@@@"
cap noi wooldid y unit t ttre , att att_it att_i att_itime cluster(unit) subgroups(sg)
dis "A_RC=" _rc
matrix A_S1 = e(wdidsgresults_1)
dis "A: att(sg1)=" %14.10f A_S1[1,1] " att_it=" %14.10f A_S1[2,1] " att_i=" %14.10f A_S1[3,1] " att_itime=" %14.10f A_S1[4,1]

* ---------------- run B: subgroup missing on some of unit 1's rows -------
dis "@@@@@@@@@@ B: sg missing on 5 of unit 1's 6 post rows @@@@@@@@@@"
cap noi wooldid y unit t ttre , att att_it att_i att_itime cluster(unit) subgroups(sg_miss)
dis "B_RC=" _rc
matrix B_S1 = e(wdidsgresults_1)
dis "B: att(sg1)=" %14.10f B_S1[1,1] " att_it=" %14.10f B_S1[2,1] " att_i=" %14.10f B_S1[3,1] " att_itime=" %14.10f B_S1[4,1]

* ---------------- the two candidate definitions ----------------
* unit 1 has 6 post rows (2004..2009), of which 1 carries sg==1
local s = 1/6
local equalunit = (5.0 + 1.0 + 1.0 + 1.0 + 1.0)/5
local sharewtd   = (`s'*5.0 + 4*1.0)/(`s' + 4)
dis ""
dis "EQUAL-UNIT att_i for sg1 (the documented definition) = " %14.10f `equalunit'
dis "SHARE-WEIGHTED att_i for sg1                          = " %14.10f `sharewtd'
dis "wooldid run B att_i                                   = " %14.10f B_S1[3,1]
dis "wooldid run B att_itime                               = " %14.10f B_S1[4,1]
dis "d(equal-unit)  = " %12.3e (B_S1[3,1]-`equalunit')
dis "d(share-wtd)   = " %12.3e (B_S1[3,1]-`sharewtd')
dis ""
dis "For comparison, run A (nothing missing) att_i = " %14.10f A_S1[3,1] " (must be 1.80)"
dis ""
dis "att (cell-weighted) is NOT affected by this: it is a row-weighted mean by"
dis "definition, so dropping rows is the right behaviour there."
dis "A att = " %14.10f A_S1[1,1] "   B att = " %14.10f B_S1[1,1]

* ---------------- run C: the case wooldid DOES refuse, for contrast -------
dis "@@@@@@@@@@ C: the adjacent case wooldid refuses @@@@@@@@@@"
gen byte sg_flip = sg
qui replace sg_flip = 0 if unit == 1 & t >= 2005
cap noi wooldid y unit t ttre , att att_it att_i att_itime cluster(unit) subgroups(sg_flip)
dis "C_RC=" _rc

* ===========================================================================
* run D -- the same thing on NOISY data, with no reliance on a known truth.
* histogramcohorteffects returns the per-unit (per-i) effect estimates, which are
* exactly the quantities att_i averages.  Build both candidate averages from
* THOSE and see which one the subgroup table reports.
* ===========================================================================
dis "@@@@@@@@@@ D: noisy data, oracle = e(histogramestimates) @@@@@@@@@@"
qui replace y = y + rnormal()*0.7
* per-unit effects from a run with NO subgroups
cap noi wooldid y unit t ttre , att_i cluster(unit) histogramcohorteffects
dis "D0_RC=" _rc
matrix HE = e(histogramestimates)
dis "HE dims: " rowsof(HE) " x " colsof(HE)
local hen : colnames HE
dis "HE colnames: `hen'"
matrix list HE

* subgroup run on the SAME noisy y
cap noi wooldid y unit t ttre , att att_it att_i att_itime cluster(unit) subgroups(sg_miss)
dis "D_RC=" _rc
matrix D_S1 = e(wdidsgresults_1)
dis "D: wooldid att_i(sg1) = " %18.14f D_S1[3,1]

mata:
  HE = st_matrix("HE")
  rs = st_matrixrowstripe("HE")
  cs = st_matrixcolstripe("HE")
  printf("rows=%f cols=%f\n", rows(HE), cols(HE))
  for (c=1; c<=cols(HE); c++) printf("col %f = %s / %s\n", c, cs[c,1], cs[c,2])
  // locate the per-unit rows for units 1..5
  b = J(5,1,.)
  for (u=1; u<=5; u++) {
    for (r=1; r<=rows(HE); r++) {
      if (strtoreal(rs[r,2])==u | rs[r,2]==strofreal(u)) {
        if (b[u]>=.) b[u] = HE[r,1]
      }
    }
  }
  b
  s = (1/6 \ 1 \ 1 \ 1 \ 1)
  equalunit = mean(b)
  sharewtd  = (s'*b)/sum(s)
  w = st_matrix("D_S1")[3,1]
  printf("\nEQUAL-UNIT  att_i(sg1) = %20.16f\n", equalunit)
  printf("SHARE-WTD   att_i(sg1) = %20.16f\n", sharewtd)
  printf("wooldid     att_i(sg1) = %20.16f\n", w)
  printf("d(equal-unit) = %10.3e   d(share-wtd) = %10.3e\n", w-equalunit, w-sharewtd)
end
