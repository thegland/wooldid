*! SKEPTIC repro B -- INDEPENDENT test of PERMSEM-01.
*! Claims under test, taken one at a time:
*!  (i)  e(permdist_b) / e(permdist) are compacted to completed replications
*!  (ii) e(perm_jointW) / e(perm_advW) keep the full requested length with missing rows
*!  (iii) "neither is indexed by permsave's _rep"      <-- the part I doubt
*! Own data: an unbalanced administrative panel (late-entering units), permute(shuffle),
*! so replication failure arises from the DESIGN, not from a probabilities scheme that
*! can draw zero treated clusters.
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do pB_rowalign.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 90210
set sortseed 90210
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-refute"

which wooldid

* ================= DGP: 16 counties, 2010-2019, 4 treated in 2015 ==============
*  four counties enter the administrative file only in 2016 (a realistic
*  late-reporting shape).  Under a shuffled assignment an onset of 2015 landing on
*  one of those counties leaves that cohort with no pre-treatment period.
clear
set obs 16
gen long cty = _n
gen double onset = cond(inlist(cty,3,7,11,14), 2015, .)
gen byte late = inlist(cty,5,9,13,16)
gen double a_i = rnormal()*2
expand 10
bysort cty: gen int year = 2010 + _n - 1
drop if late & year < 2016
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.25*(year-2010) + 1.0*post + rnormal()
tab late
count

cap noi wooldid y cty year onset, cluster(cty) att esprelength(2) espostlength(2) ///
    jointtests adversarial esfixedbaseperiod permute(shuffle) permreps(100) ///
    rseed(777) permsave(pB_save) replace
dis "PB_RC=" _rc
dis "PB_reps="  e(perm_reps) "  completed=" e(perm_completed) "  failed=" e(perm_failed)

foreach m in permdist_b permdist perm_jointW perm_advW {
  cap confirm matrix e(`m')
  if _rc == 0 {
    matrix TMP = e(`m')
    dis "PB_DIM e(`m') = " rowsof(TMP) " x " colsof(TMP)
  }
  else dis "PB_DIM e(`m') = ABSENT"
}

cap confirm matrix e(perm_jointW)
local hasjw = (_rc == 0)
if `hasjw' matrix _JW = e(perm_jointW)
cap confirm matrix e(perm_advW)
if _rc == 0 matrix _AW = e(perm_advW)
matrix _PDB = e(permdist_b)

* ------------- the alignment test the claim turns on -----------------
preserve
use pB_save.dta, clear
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
* put the file's failure indicator into a Mata vector, in _rep order
sort _rep
qui drop if _rep == 0
mkmat allmiss, matrix(_FILEMISS)
restore

mata:
D  = st_matrix("_PDB")
FM = st_matrix("_FILEMISS")
printf("PB_MATA permdist_b rows = %g ; permsave replication rows = %g ; file-failed = %g\n",
       rows(D), rows(FM), sum(FM))
if (st_local("hasjw") == "1") {
  JW = st_matrix("_JW")
  printf("PB_MATA perm_jointW rows = %g\n", rows(JW))
  if (rows(JW) == rows(FM)) {
    jwmiss = (JW[.,1] :>= .)
    printf("PB_ALIGN jointW-missing rows vs permsave-failed rows: agree=%g of %g ; JWmiss=%g FILEmiss=%g\n",
           sum(jwmiss :== FM[.,1]), rows(FM), sum(jwmiss), sum(FM))
    if (sum(jwmiss :== FM[.,1]) == rows(FM)) printf("PB_VERDICT perm_jointW IS row-indexed by permsave _rep\n")
    else printf("PB_VERDICT perm_jointW is NOT row-indexed by permsave _rep\n")
  }
  else printf("PB_VERDICT perm_jointW row count differs from permsave replication rows\n")
}
* does permdist_b line up with the SURVIVING permsave rows?
if (rows(D) == rows(FM) - sum(FM)) printf("PB_PDB permdist_b rows == surviving permsave rows (compacted)\n")
else printf("PB_PDB permdist_b rows = %g, surviving permsave rows = %g (MISMATCH)\n", rows(D), rows(FM)-sum(FM))
end

* ---- does the compaction of permdist_b actually equal the surviving file rows? ----
preserve
use pB_save.dta, clear
qui drop if _rep == 0
qui ds _rep, not
local coefvars `r(varlist)'
qui gen byte anyok = 0
foreach v of local coefvars {
  qui replace anyok = 1 if !missing(`v')
}
qui keep if anyok == 1
sort _rep
mkmat `coefvars', matrix(_FILEOK)
restore
mata:
FO = st_matrix("_FILEOK")
D  = st_matrix("_PDB")
if (rows(FO) == rows(D) & cols(FO) == cols(D)) {
  printf("PB_EQ max abs diff permdist_b vs compacted permsave block = %g\n", max(abs(editmissing(FO,0) - editmissing(D,0))))
}
else printf("PB_EQ shape mismatch: file %g x %g vs permdist_b %g x %g\n", rows(FO), cols(FO), rows(D), cols(D))
end

dis "PB_DONE"
