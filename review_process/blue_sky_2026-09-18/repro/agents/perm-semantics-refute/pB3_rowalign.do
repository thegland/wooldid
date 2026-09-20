*! SKEPTIC repro B3 -- INDEPENDENT test of PERMSEM-01 with real replication failures.
*! Failure mechanism: a panel with attrition.  Ten districts, two adopt in 2017; three
*! districts stop reporting after 2014.  A shuffled assignment that puts BOTH onsets on
*! the three early leavers leaves the replication with no treated observation at all, so
*! that replication's fit fails.  P(both on leavers) = C(3,2)/C(10,2) = 3/45.
*! (No permprobs, no zero-treated-by-construction probability scheme.)
*!
*! Claims under test:
*!  (i)   e(permdist_b)/e(permdist) are compacted to the completed replications
*!  (ii)  e(perm_jointW)/e(perm_advW) keep the full requested length with missing rows
*!  (iii) "neither is indexed by permsave's _rep"       <-- the part I doubt
*!  (iv)  a naive (1+#{W_r>=Wobs})/(N+1) on the FULL-length e(perm_jointW) over-counts
*!        by exactly e(perm_failed)
*! NOTE on the permsave file: with permsave() the file ALSO carries one onset_<id> column
*! per cluster, and those are populated on failed replications too, so a "failed row" must
*! be detected on the COEFFICIENT columns only (the columns of e(permdist_b)).
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do pB3_rowalign.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-refute"

which wooldid

clear
set obs 10
gen long dist = _n
gen double onset = cond(inlist(dist,4,9), 2017, .)
gen byte leaver = inlist(dist,2,6,10)
gen double a_i = rnormal()*2
expand 10
bysort dist: gen int year = 2010 + _n - 1
drop if leaver & year > 2014
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.25*(year-2010) + 1.0*post + rnormal()
count
tab leaver

cap noi wooldid y dist year onset, cluster(dist) att esprelength(2) espostlength(2) ///
    jointtests adversarial esfixedbaseperiod permute(shuffle) permreps(100) ///
    rseed(99) permsave(pB3_save) replace
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

matrix _PDB = e(permdist_b)
local pdbnames : colnames _PDB
dis "PB_PDBNAMES = `pdbnames'"
cap confirm matrix e(perm_jointW)
local hasjw = (_rc == 0)
if `hasjw' {
  matrix _JW  = e(perm_jointW)
  matrix _JWO = e(perm_jointWobs)
}
scalar _NFAIL = e(perm_failed)

preserve
use pB3_save.dta, clear
dis "PB_SAVE_ROWS = " _N
qui ds onset_*
local onsetvars `r(varlist)'
dis "PB_SAVE_ONSETCOLS = `: word count `onsetvars''"
qui ds _rep onset_*, not
local coefvars `r(varlist)'
dis "PB_SAVE_COEFCOLS = `: word count `coefvars''"
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
printf("PB_MATA permdist_b rows=%g cols=%g ; permsave rep rows=%g cols=%g ; file-failed=%g\n",
       rows(D), cols(D), rows(FM), cols(FA), sum(FM))
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
  else printf("PB_VERDICT3 perm_jointW rows %g != permsave rep rows %g\n", rows(JW), rows(FM))
  nn    = sum(JW[.,1] :< .)
  nex_c = sum(select(JW[.,1], JW[.,1] :< .) :>= JWO[1,1])
  nex_n = sum(JW[.,1] :>= JWO[1,1])
  printf("PB_RECON correct: n=%g nex=%g p=%s\n", nn, nex_c, strofreal((1+nex_c)/(nn+1), "%21x"))
  printf("PB_RECON naive  : N=%g nex=%g p=%s\n", rows(JW), nex_n, strofreal((1+nex_n)/(rows(JW)+1), "%21x"))
  printf("PB_RECON posted joint_horiz=%s overcount=%g e(perm_failed)=%g\n",
         strofreal(st_numscalar("e(joint_horiz)"), "%21x"), nex_n - nex_c, st_numscalar("_NFAIL"))
}
ok = (rowsum(FA :< .) :> 0)
FO = select(FA, ok)
if (rows(FO) == rows(D) & cols(FO) == cols(D))
  printf("PB_EQ permdist_b == compacted permsave coefficient block, maxabsdiff=%g\n",
         max(abs(editmissing(FO,0) - editmissing(D,0))))
else
  printf("PB_EQ shape mismatch: file %g x %g vs permdist_b %g x %g\n", rows(FO), cols(FO), rows(D), cols(D))
end

dis "PB_DONE"
