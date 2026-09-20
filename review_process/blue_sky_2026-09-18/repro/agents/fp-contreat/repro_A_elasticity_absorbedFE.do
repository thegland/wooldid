* ===========================================================================
* FINDING A  (self-contained repro; run:  stata-mp -e do <thisfile>)
*
* CLAIM: on wooldid's DEFAULT estimation path (reghdfe), the average marginal
* effects requested with contreatelasticitytype(eydx) and (eyex) divide by a
* prediction that OMITS the absorbed fixed effects (i.i, i.t and the absorbed
* dose control).  The reported number is therefore not the "% change in y"
* the help file defines (sthlp:905-907), and it differs from the same
* estimand computed on the plain-regress path (customvce()/unconditionalse),
* where the fixed effects are explicit regressors and the answer is exact.
*
* ORACLE: the panel below is EXACT-FIT - y is produced with zero residual by a
* model nested in wooldid's continuous-treatment specification - so
* E[y|x] == y and the truths are closed form over the treated post rows:
*     dydx = mean(beta_it)          eydx = mean(beta_it / y)
*     dyex = mean(beta_it*dose)     eyex = mean(beta_it*dose / y)
* 20 states (the cohort id i), 5 counties each, 10 periods, staggered onsets.
* ===========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 60601
set sortseed 60601

clear
set obs 1000                       /* 20 states x 5 counties x 10 periods */
gen long unit = ceil(_n/10)
gen int  t    = _n - (unit-1)*10
gen byte i    = ceil(unit/5)       /* 20 states, 5 counties each          */
gen double ttre = .
qui replace ttre = 4 if inrange(i, 1, 5)
qui replace ttre = 6 if inrange(i, 6,10)
qui replace ttre = 8 if inrange(i,11,15)
* states 16-20 never treated

* state-level level differences -- exactly what reghdfe absorbs as i.i
gen double sfe = 40 + 3*i
gen double dose = 1 + 2*runiform()          /* policy intensity, strictly > 0 */

gen double kappatrue = 0
gen double betatrue  = 0
qui replace kappatrue = 1.0 + 0.20*(t-ttre) if ttre < . & t >= ttre
qui replace betatrue  = 2.0 - 0.30*(t-ttre) if ttre < . & t >= ttre

* EXACT FIT for:  absorb(i.i i.t) + contreatcontrols(basic) + cell##c.dose
gen double y = sfe + 0.5*t + 0.9*dose + kappatrue + betatrue*dose

* ---------------- closed-form truths ---------------------------------------
gen double q_dydx = betatrue
gen double q_dyex = betatrue*dose
gen double q_eydx = betatrue/y
gen double q_eyex = betatrue*dose/y
foreach k in dydx dyex eydx eyex {
  qui sum q_`k' if ttre < . & t >= ttre
  local T_`k' = r(mean)
  dis as result "TRUTH `k' = " %18.12f `T_`k''
}

foreach et in dydx dyex eydx eyex {
  foreach path in "cluster(unit)" "customvce(cluster unit)" {
    cap qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) ///
        contreatelasticitytype(`et') `path' regtol(14)
    if _rc == 0 {
      mat M = e(wdidmainresults)
      dis as result "RESULT|`et'|`path'| ame = " %18.12f M[2,1] ///
        " | truth = " %18.12f `T_`et'' " | rel.err = " ///
        %10.3e abs(M[2,1]-`T_`et'')/abs(`T_`et'') " | R2 = " %8.6f e(r2)
    }
    else dis as error "RESULT|`et'|`path'| RC=" _rc
  }
}
dis "DONE"
