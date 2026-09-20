*=========================================================================
* repro_S1_impvalues.do  -- SELF-CONTAINED.
*
* CLAIM DEMONSTRATED: under `impvalues' with a cluster() FINER than i, the
* reported Ibragimov-Muller standard error is the dispersion of only as many
* DISTINCT per-cluster reads as there are treated cohorts, but divided by the
* square root of the number of CLUSTERS and referred to a t on (clusters - 1)
* degrees of freedom.  The standard error is therefore deflated by
* sqrt(clusters / distinct reads) and the 95% interval covers about 25%.
*
* Three parts, ~4 minutes total:
*   PART 1  three arms on ONE dataset + an exact arithmetic identity that
*           proves the per-cluster reads take only 3 distinct values.
*   PART 2  60-draw Monte Carlo: coverage and mean(se) vs sd(b).
*   PART 3  the crossed-cluster corollary: every se/t/p/CI blanked with a
*           message blaming "deterministic or noiseless data", on data whose
*           idiosyncratic sd is 0.55.
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_S1_impvalues.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

*---- generator: 60 units, i = cg (3 onset cohorts + never-treated),
*     12 periods, unit random intercept, AR(1) errors, within-cohort
*     effect heterogeneity.  Nothing about this panel is degenerate.
cap prog drop mkdat
prog def mkdat
  args nunit sd
  clear
  set seed `sd'
  set sortseed `sd'
  set obs `nunit'
  gen long unit = _n
  gen byte cg  = mod(_n-1,4)+1
  gen byte reg = mod(_n-1,5)+1
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

*========================== PART 1 ======================================
di as txt _n(2) "################ PART 1: three arms, one dataset ################"
mkdat 60 9021
qui su e
di as txt "sd of the idiosyncratic error in this panel = " as res r(sd)
tempfile d1
qui save "`d1'"

di as txt _n "--- ARM 1: impvalues, cluster(unit)  [cluster FINER than i=cg] ---"
use "`d1'", clear
wooldid y cg t ttre, att cluster(unit) impvalues level(95)
scalar se1 = e(wdidmainresults)[1,2]
scalar b1  = e(wdidmainresults)[1,1]
scalar df1 = e(df_r)

di as txt _n "--- ARM 2: impvalues, cluster(cg)  [cluster == i: textbook IM] ---"
use "`d1'", clear
wooldid y cg t ttre, att cluster(cg) impvalues level(95)
scalar se2 = e(wdidmainresults)[1,2]
scalar b2  = e(wdidmainresults)[1,1]
scalar df2 = e(df_r)

di as txt _n "--- ARM 3: ordinary analytic clustered run, cluster(unit) ---"
use "`d1'", clear
wooldid y cg t ttre, att cluster(unit) level(95)
scalar se3 = e(wdidmainresults)[1,2]
scalar b3  = e(wdidmainresults)[1,1]
scalar df3 = e(df_r)

di as txt _n "================= PART 1 SUMMARY ================="
di as txt "impvalues cluster(unit) : b = " %12.9f b1 "  se = " %12.9f se1 "  df_r = " df1
di as txt "impvalues cluster(cg)   : b = " %12.9f b2 "  se = " %12.9f se2 "  df_r = " df2
di as txt "analytic  cluster(unit) : b = " %12.9f b3 "  se = " %12.9f se3 "  df_r = " df3
di as txt "se ratio  cluster(cg)/cluster(unit) = " %8.4f se2/se1
di as txt "se ratio  analytic/impvalues-unit   = " %8.4f se3/se1

*---- the arithmetic identity ------------------------------------------
* If the 45 per-unit reads take only 3 distinct values (one per treated
* cohort) with multiplicity 15, then
*     sd over 45 reads = sd over 3 reads * sqrt(2*15/44)
* and the reported IM se is that divided by sqrt(45).  se2 gives the sd
* over the 3 reads directly (se2 = sd3/sqrt(3)).
scalar sd3  = se2*sqrt(3)
scalar pred = sd3*sqrt(2*15/44)/sqrt(45)
di as txt _n "PREDICTED impvalues cluster(unit) se under the "
di as txt "  '3 distinct reads, each repeated 15 times' hypothesis = " %14.11f pred
di as txt "OBSERVED                                                = " %14.11f se1
di as txt "ratio observed/predicted                                = " %18.15f se1/pred
di as txt "(a ratio of 1 to 8+ digits is only possible if the 45 reads really"
di as txt " do take exactly 3 distinct values.)"

*========================== PART 2 ======================================
di as txt _n(2) "################ PART 2: 60-draw Monte Carlo ################"
tempname MB
matrix `MB' = J(60,10,.)
forvalues r = 1/60 {
  mkdat 60 `=91000+`r''
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
  qui wooldid y cg t ttre, att cluster(unit) impvalues level(95)
  mat T = e(wdidmainresults)
  matrix `MB'[`r',1] = T[1,1]
  matrix `MB'[`r',2] = T[1,2]
  matrix `MB'[`r',3] = T[1,5]
  matrix `MB'[`r',4] = T[1,6]
  use "`dd'", clear
  qui wooldid y cg t ttre, att cluster(unit) level(95)
  mat T = e(wdidmainresults)
  matrix `MB'[`r',5] = T[1,1]
  matrix `MB'[`r',6] = T[1,2]
  matrix `MB'[`r',7] = T[1,5]
  matrix `MB'[`r',8] = T[1,6]
  matrix `MB'[`r',9]  = `tatt'
  matrix `MB'[`r',10] = `tceq'
}
clear
qui svmat double `MB'
rename (`MB'1 `MB'2 `MB'3 `MB'4 `MB'5 `MB'6 `MB'7 `MB'8 `MB'9 `MB'10) ///
       (bI seI lbI ubI bC seC lbC ubC tatt tceq)
gen byte covI_att = (lbI<=tatt & tatt<=ubI)
gen byte covI_ceq = (lbI<=tceq & tceq<=ubI)
gen byte covC_att = (lbC<=tatt & tatt<=ubC)
qui su bI
local sdI = r(sd)
local mI  = r(mean)
qui su seI
local mseI = r(mean)
qui su bC
local sdC = r(sd)
qui su seC
local mseC = r(mean)
gen byte covI_mc = (lbI<=`mI' & `mI'<=ubI)
di as txt _n "================= PART 2 SUMMARY (60 draws) ================="
di as txt "impvalues cluster(unit): mean(b) = " %9.5f `mI' "  sd(b) = " %9.5f `sdI' ///
          "  mean(se) = " %9.5f `mseI' "   se/sd = " %6.3f `mseI'/`sdI'
di as txt "analytic  cluster(unit):                  sd(b) = " %9.5f `sdC' ///
          "  mean(se) = " %9.5f `mseC' "   se/sd = " %6.3f `mseC'/`sdC'
qui su covI_att
di as txt "impvalues 95% coverage vs the att truth                 = " %6.3f r(mean)
qui su covI_ceq
di as txt "impvalues 95% coverage vs the cluster-equal-weight truth= " %6.3f r(mean)
qui su covI_mc
di as txt "impvalues 95% coverage around its OWN Monte Carlo mean  = " %6.3f r(mean)
qui su covC_att
di as txt "analytic  95% coverage vs the att truth                 = " %6.3f r(mean)
di as txt "3 binomial SE at p=.95, n=60 = " %6.3f 3*sqrt(.95*.05/60)

*========================== PART 3 ======================================
di as txt _n(2) "########## PART 3: crossed cluster -> everything blanked ##########"
mkdat 60 9021
* reg has 5 levels, each holding units from every cohort in equal numbers,
* so every cluster's weight profile over the (i x relative period) cells is
* identical and the per-cluster reads collapse to a single value.
wooldid y cg t ttre, att cluster(reg) impvalues level(95)
di as txt "^ the message blames 'deterministic or noiseless data'; the"
di as txt "  idiosyncratic sd of this panel is 0.55 and the unit sd is 0.70."

di as txt _n "ALLDONE_S1"
