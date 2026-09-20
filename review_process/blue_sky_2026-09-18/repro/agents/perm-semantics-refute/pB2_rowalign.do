*! SKEPTIC repro B2 -- INDEPENDENT test of PERMSEM-01, with replication failures
*! produced by a REALISTIC unbalanced-panel geometry rather than a probabilities scheme.
*! 18 counties, staggered adoption with FOUR DISTINCT onset years (one county each), and
*! an administrative file in which four counties report only from 2016 and three stop
*! reporting after 2014.  A shuffled onset landing on a late reporter leaves that
*! single-county cohort with no pre-treatment period; landing on an early leaver leaves
*! it with no post-treatment period.  Both make the replication's fit fail.
*! Claims under test:
*!  (i)   e(permdist_b)/e(permdist) are compacted to completed replications
*!  (ii)  e(perm_jointW)/e(perm_advW) keep full requested length with missing rows
*!  (iii) "neither is indexed by permsave's _rep"      <-- the part I doubt
*!  (iv)  the "sharper case": a naive (1+sum(W:>=Wobs))/(n+1) on the posted full-length
*!        e(perm_jointW) over-counts by exactly e(perm_failed)
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do pB2_rowalign.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 31337
set sortseed 31337
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-refute"

which wooldid

clear
set obs 18
gen long cty = _n
gen double onset = .
replace onset = 2013 if cty == 2
replace onset = 2015 if cty == 6
replace onset = 2016 if cty == 10
replace onset = 2017 if cty == 15
gen byte late  = inlist(cty,4,8,12,17)
gen byte early = inlist(cty,5,13,18)
gen double a_i = rnormal()*2
expand 10
bysort cty: gen int year = 2010 + _n - 1
drop if late  & year < 2016
drop if early & year > 2014
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.25*(year-2010) + 1.0*post + rnormal()
count

cap noi wooldid y cty year onset, cluster(cty) att esprelength(2) espostlength(2) ///
    jointtests adversarial esfixedbaseperiod permute(shuffle) permreps(80) ///
    rseed(2468) permsave(pB2_save) replace
dis "PB_RC=" _rc
dis "PB_reps="  e(perm_reps) "  completed=" e(perm_completed) "  failed=" e(perm_failed)

foreach m in permdist_b permdist perm_jointW perm_advW perm_jointWobs perm_advWobs {
  cap confirm matrix e(`m')
  if _rc == 0 {
    matrix TMP = e(`m')
    dis "PB_DIM e(`m') = " rowsof(TMP) " x " colsof(TMP)
  }
  else dis "PB_DIM e(`m') = ABSENT"
}
dis "PB_JOINT_HORIZ p = " %21x e(joint_horiz)
dis "PB_ADV p = " %21x e(adversarial_pval_att)

matrix _PDB = e(permdist_b)
cap confirm matrix e(perm_jointW)
local hasjw = (_rc == 0)
if `hasjw' {
  matrix _JW  = e(perm_jointW)
  matrix _JWO = e(perm_jointWobs)
}
scalar _NFAIL = e(perm_failed)

preserve
use pB2_save.dta, clear
dis "PB_SAVE_ROWS = " _N
qui ds _rep, not
local coefvars `r(varlist)'
qui gen byte allmiss = 1
foreach v of local coefvars {
  qui replace allmiss = 0 if !missing(`v')
}
qui count if _rep >= 1 & allmiss == 1
dis "PB_SAVE_FAILEDROWS = " r(N)
qui levelsof _rep if _rep >= 1 & allmiss == 1, local(failreps)
dis "PB_SAVE_FAILED_REPIDS = `failreps'"
sort _rep
qui drop if _rep == 0
mkmat allmiss, matrix(_FILEMISS)
mkmat `coefvars', matrix(_FILEALL)
restore

mata:
D  = st_matrix("_PDB")
FM = st_matrix("_FILEMISS")
FA = st_matrix("_FILEALL")
printf("PB_MATA permdist_b rows=%g ; permsave rep rows=%g ; file-failed=%g\n",
       rows(D), rows(FM), sum(FM))
if (st_local("hasjw") == "1") {
  JW  = st_matrix("_JW")
  JWO = st_matrix("_JWO")
  printf("PB_MATA perm_jointW rows=%g\n", rows(JW))
  if (rows(JW) == rows(FM)) {
    jwmiss = (JW[.,1] :>= .)
    printf("PB_ALIGN agree=%g of %g ; JWmiss=%g FILEmiss=%g\n",
           sum(jwmiss :== FM[.,1]), rows(FM), sum(jwmiss), sum(FM))
    if (sum(jwmiss :== FM[.,1]) == rows(FM))
      printf("PB_VERDICT3 perm_jointW IS row-indexed by permsave _rep\n")
    else
      printf("PB_VERDICT3 perm_jointW is NOT row-indexed by permsave _rep\n")
  }
  else printf("PB_VERDICT3 perm_jointW row count %g != permsave rep rows %g\n", rows(JW), rows(FM))
  // (iv) the naive vs correct reconstruction of the joint p from the posted pair
  nn    = sum(JW[.,1] :< .)
  nex_c = sum(select(JW[.,1], JW[.,1] :< .) :>= JWO[1,1])
  nex_n = sum(JW[.,1] :>= JWO[1,1])
  printf("PB_RECON correct: n=%g nex=%g p=%s\n", nn, nex_c, strofreal((1+nex_c)/(nn+1), "%21x"))
  printf("PB_RECON naive  : N=%g nex=%g p=%s\n", rows(JW), nex_n, strofreal((1+nex_n)/(rows(JW)+1), "%21x"))
  printf("PB_RECON posted joint_horiz = %s ; overcount = %g ; e(perm_failed) = %g\n",
         strofreal(st_numscalar("e(joint_horiz)"), "%21x"), nex_n - nex_c, st_numscalar("_NFAIL"))
}
// does permdist_b equal the surviving permsave rows, in order?
ok = (rowsum(FA :< .) :> 0)
FO = select(FA, ok)
if (rows(FO) == rows(D) & cols(FO) == cols(D))
  printf("PB_EQ permdist_b == compacted permsave block, maxabsdiff=%g\n",
         max(abs(editmissing(FO,0) - editmissing(D,0))))
else
  printf("PB_EQ shape mismatch: file %g x %g vs permdist_b %g x %g\n", rows(FO), cols(FO), rows(D), cols(D))
end

dis "PB_DONE"
