*=========================================================================
* repro_S4_imunit.do -- SELF-CONTAINED, ~8 minutes.
*
* CLAIM DEMONSTRATED: the `impvalues' under-coverage is NOT confined to the
* geometry where a cluster is finer than i.  With i = unit and cluster(unit)
* -- 45 treated clusters, each owning its own cohort-by-relative-period
* coefficients, which is the geometry in which Ibragimov-Muller is textbook-
* legitimate -- the reported IM standard error is still only about 2/3 of the
* estimator's actual sampling standard deviation and the nominal 95% interval
* covers about 83%.  On the same draws the analytic clustered interval covers
* about 89% (itself the documented at-or-above-i understatement).
*
* This is a DIFFERENT mechanism from the finer-than-i one: here the 45 reads
* really are 45 distinct numbers, but they are positively correlated, because
* every treated unit's effect is measured against the same never-treated
* comparison group and the same estimated period effects, and IM's
* sd/sqrt(G) assumes independence across clusters.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_S4_imunit.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local R 80

cap prog drop mkdat
prog def mkdat
  args nunit sd
  clear
  set seed `sd'
  set sortseed `sd'
  set obs `nunit'
  gen long unit = _n
  gen byte cg = mod(_n-1,4)+1
  gen double ttre = .
  replace ttre = 2004 if cg==1
  replace ttre = 2006 if cg==2
  replace ttre = 2008 if cg==3
  gen double ui   = rnormal(0,0.70)
  gen double hmul = 1 + rnormal(0,0.25)
  expand 12
  bys unit: gen int t = 2000 + _n - 1
  gen double e = .
  bys unit (t): replace e = rnormal(0,0.55) if _n==1
  bys unit (t): replace e = 0.4*e[_n-1] + rnormal(0,0.55) if _n>1
  gen double truetau = 0
  replace truetau = 1.00*hmul if cg==1 & t>=ttre
  replace truetau = 1.15*hmul if cg==2 & t>=ttre
  replace truetau = 0.85*hmul if cg==3 & t>=ttre
  gen double y = ui + 0.25*(t-2000) + truetau + e
end

tempname M
matrix `M' = J(`R', 10, .)
forvalues r = 1/`R' {
  mkdat 60 `=70000+`r''
  qui su truetau if !missing(ttre) & t>=ttre
  local tatt = r(mean)
  tempvar um
  qui egen double `um' = mean(cond(!missing(ttre) & t>=ttre, truetau, .)), by(unit)
  qui bys unit: gen byte __f = (_n==1)
  qui su `um' if !missing(ttre) & __f==1, meanonly
  local tceq = r(mean)
  drop __f
  tempfile dd
  qui save "`dd'"
  use "`dd'", clear
  cap qui wooldid y unit t ttre, att cluster(unit) impvalues level(95)
  if _rc==0 {
    mat T = e(wdidmainresults)
    matrix `M'[`r',1]=T[1,1]
    matrix `M'[`r',2]=T[1,2]
    matrix `M'[`r',3]=T[1,5]
    matrix `M'[`r',4]=T[1,6]
  }
  use "`dd'", clear
  cap qui wooldid y unit t ttre, att cluster(unit) level(95)
  if _rc==0 {
    mat T = e(wdidmainresults)
    matrix `M'[`r',5]=T[1,1]
    matrix `M'[`r',6]=T[1,2]
    matrix `M'[`r',7]=T[1,5]
    matrix `M'[`r',8]=T[1,6]
  }
  matrix `M'[`r',9]=`tatt'
  matrix `M'[`r',10]=`tceq'
}
clear
qui svmat double `M'
rename (`M'1 `M'2 `M'3 `M'4 `M'5 `M'6 `M'7 `M'8 `M'9 `M'10) ///
       (bI seI lbI ubI bC seC lbC ubC tatt tceq)
qui su bI
local mI=r(mean)
local sdI=r(sd)
qui su seI
local msI=r(mean)
qui su bC
local sdC=r(sd)
qui su seC
local msC=r(mean)
gen byte cIa = (lbI<=tatt & tatt<=ubI)
gen byte cIc = (lbI<=tceq & tceq<=ubI)
gen byte cIm = (lbI<=`mI'  & `mI' <=ubI)
gen byte cCa = (lbC<=tatt & tatt<=ubC)
di as txt _n(2) "============ repro_S4 SUMMARY ( `R' draws, i = unit, cluster(unit) ) ============"
di as txt "impvalues : mean(b)=" %8.4f `mI' "  sd(b)=" %8.4f `sdI' "  mean(se)=" %8.4f `msI' "   se/sd=" %6.3f `msI'/`sdI'
di as txt "analytic  :                   sd(b)=" %8.4f `sdC' "  mean(se)=" %8.4f `msC' "   se/sd=" %6.3f `msC'/`sdC'
qui su cIa
di as txt "impvalues 95% coverage vs att truth                = " %6.3f r(mean)
qui su cIc
di as txt "impvalues 95% coverage vs cluster-equal-weight truth= " %6.3f r(mean)
qui su cIm
di as txt "impvalues 95% coverage around its own MC mean       = " %6.3f r(mean)
qui su cCa
di as txt "analytic  95% coverage vs att truth                 = " %6.3f r(mean)
di as txt "3 binomial SE at p=.95, n=`R' = " %6.3f 3*sqrt(.95*.05/`R')
di as txt _n "ALLDONE_S4"
