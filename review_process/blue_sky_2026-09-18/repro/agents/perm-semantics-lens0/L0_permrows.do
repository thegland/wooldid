*! LENS-0 (independent reproduction) of PERMSEM-01.
*! Third, independent failure mechanism: a ROLLING ADMINISTRATIVE FRAME.
*! 30 districts observed 2010-2019.  12 report the whole window; 10 enter the frame in
*! 2015; 8 enter in 2017.  Three districts adopt in 2017.  The observed assignment puts
*! all three onsets on full-window districts, so the observed fit carries every horizon.
*! A shuffled assignment that lands all three onsets on late entrants leaves the deepest
*! pre-treatment horizon (relyr -3) with no contributing unit, so that replication
*! COMPLETES but its -3 coefficient column is missing; one that lands them all on the
*! 2017 entrants leaves no pre-period at all.
*! This gives both row types the two earlier repros did not separate:
*!   (a) fully failed replications  (all coefficient columns missing)
*!   (b) PARTIAL replications       (completed, but missing on the JOINTLY TESTED columns)
*! Claims under test, stated from the finding alone:
*!   C1 e(permdist_b)/e(permdist) are compacted to completed replications
*!   C2 e(perm_jointW)/e(perm_advW) keep the full REQUESTED length with missing rows
*!   C3 "neither is indexed by permsave's _rep"
*!   C4 naive (1+#{W_r>=Wobs})/(N+1) on full-length e(perm_jointW) over-counts by
*!      EXACTLY e(perm_failed)
*!   C5 no printed number is wrong (posted p equals the correct hand computation)
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do L0_permrows.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 771103
set sortseed 771103
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-lens0"

which wooldid

* ---------------------------------------------------------------- data
clear
set obs 30
gen long dist = _n
* entry cohort of the administrative frame
gen int firstyr = 2010
replace firstyr = 2015 if inrange(dist, 13, 22)
replace firstyr = 2017 if inrange(dist, 23, 30)
* three adopters, all full-window districts in the OBSERVED assignment
gen double onset = cond(inlist(dist, 3, 7, 11), 2017, .)
gen double a_i = rnormal()*1.5
expand 10
bysort dist: gen int year = 2010 + _n - 1
drop if year < firstyr
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.20*(year - 2010) + 0.8*post + rnormal()
qui count
dis "L0_NOBS=" r(N)
qui levelsof dist if firstyr == 2015, local(l15)
qui levelsof dist if firstyr == 2017, local(l17)
dis "L0_LATE2015=`l15'"
dis "L0_LATE2017=`l17'"

* ---------------------------------------------------------------- the run
cap noi wooldid y dist year onset, cluster(dist) att esprelength(3) espostlength(2) ///
    jointtests adversarial esfixedbaseperiod permute(shuffle) permreps(120) ///
    rseed(4041) permsave(L0_save) replace
dis "L0_RC=" _rc
dis "L0_reps=" e(perm_reps) " completed=" e(perm_completed) " failed=" e(perm_failed)
dis "L0_pconv=" e(perm_pconv) " exact=" e(perm_exact) " npossible=" e(perm_npossible_str)

foreach m in permdist_b permdist perm_jointW perm_jointWobs perm_advW perm_advWobs perm_n {
  cap confirm matrix e(`m')
  if _rc == 0 {
    matrix TMP = e(`m')
    dis "L0_DIM e(`m') = " rowsof(TMP) " x " colsof(TMP)
  }
  else dis "L0_DIM e(`m') = ABSENT"
}

* per-element replication counts: if these all equal e(perm_completed) there are no
* partial rows in this run
cap confirm matrix e(perm_n)
if _rc == 0 {
  matrix _PN = e(perm_n)
  mata: st_local("pnmin", strofreal(min(st_matrix("_PN")))); st_local("pnmax", strofreal(max(st_matrix("_PN"))))
  dis "L0_PERMN min=`pnmin' max=`pnmax' completed=" e(perm_completed)
}

* names of the joint / adversarial tests actually posted
cap confirm matrix e(perm_jointW)
local hasjw = (_rc == 0)
if `hasjw' {
  matrix _JW = e(perm_jointW)
  matrix _JWO = e(perm_jointWobs)
  local jwn : colnames _JW
  dis "L0_JWNAMES=`jwn'"
}
cap confirm matrix e(perm_advW)
local hasaw = (_rc == 0)
if `hasaw' {
  matrix _AW = e(perm_advW)
  matrix _AWO = e(perm_advWobs)
  local awn : colnames _AW
  dis "L0_AWNAMES=`awn'"
}
matrix _PDB = e(permdist_b)
scalar _NREPS = e(perm_reps)
scalar _NCOMP = e(perm_completed)
scalar _NFAIL = e(perm_failed)

* the posted p's, in hex, for every joint / adversarial name
qui ereturn list
foreach s in joint_horiz joint_horiz_att joint_pre {
  cap dis "L0_POSTED `s' = " %21x e(`s')
}

* ---------------------------------------------------------------- the saved file
preserve
cap noi use L0_save.dta, clear
dis "L0_USE_RC=" _rc
dis "L0_SAVE_ROWS=" _N
qui ds onset_*
local onsetvars `r(varlist)'
dis "L0_SAVE_ONSETCOLS=`: word count `onsetvars''"
qui ds _rep onset_*, not
local coefvars `r(varlist)'
local ncoef : word count `coefvars'
dis "L0_SAVE_COEFCOLS=`ncoef'"
qui gen int nmiss = 0
foreach v of local coefvars {
  qui replace nmiss = nmiss + missing(`v')
}
qui gen byte allmiss = (nmiss == `ncoef')
qui gen byte partial = (nmiss > 0 & nmiss < `ncoef')
qui count if _rep >= 1 & allmiss
dis "L0_FILE_FAILEDROWS=" r(N)
qui count if _rep >= 1 & partial
dis "L0_FILE_PARTIALROWS=" r(N)
qui levelsof _rep if _rep >= 1 & allmiss, local(fr)
dis "L0_FILE_FAILED_REPIDS=`fr'"
qui levelsof _rep if _rep >= 1 & partial, local(pr)
dis "L0_FILE_PARTIAL_REPIDS=`pr'"
qui sum _rep
dis "L0_FILE_MAXREP=" r(max)
sort _rep
qui drop if _rep == 0
mkmat allmiss, matrix(_FMISS)
mkmat `coefvars', matrix(_FALL)
restore

* ---------------------------------------------------------------- the arithmetic
mata:
D  = st_matrix("_PDB")
FM = st_matrix("_FMISS")
FA = st_matrix("_FALL")
NR = st_numscalar("_NREPS"); NC = st_numscalar("_NCOMP"); NF = st_numscalar("_NFAIL")
printf("L0_BASES permdist_b %gx%g | permsave rep rows %g | reps=%g completed=%g failed=%g\n",
       rows(D), cols(D), rows(FM), NR, NC, NF)

// C1: permdist_b == the file's non-all-missing coefficient block, in order
ok = (rowsum(FA :< .) :> 0)
FC = select(FA, ok)
printf("L0_C1 compacted-file rows=%g permdist_b rows=%g\n", rows(FC), rows(D))
if (rows(FC) == rows(D) & cols(FC) == cols(D)) {
  dd = FC - D
  printf("L0_C1 maxabsdiff(permdist_b, compacted file block)=%g\n", max(abs(select(dd, dd :< .))))
}

// C2/C3 on the jointW half
if (st_local("hasjw") == "1") {
  JW = st_matrix("_JW"); JWO = st_matrix("_JWO")
  printf("L0_C2 rows(perm_jointW)=%g  reps=%g  filereprows=%g\n", rows(JW), NR, rows(FM))
  jwm = (JW[.,1] :>= .)
  // necessary condition for _rep indexing: every all-missing FILE row is missing in JW
  bad = sum(jwm :< FM[.,1])          // file failed but JW present -> impossible if aligned
  extra = sum(jwm :> FM[.,1])        // JW missing on a row the file kept -> PARTIAL row
  printf("L0_C3 rows=%g  JWmiss=%g  FILEallmiss=%g  violations(JW present on failed rep)=%g  extra(JW missing on kept rep)=%g\n",
         rows(JW), sum(jwm), sum(FM), bad, extra)
  if (bad == 0) printf("L0_C3_VERDICT perm_jointW IS row-indexed by permsave _rep (every failed rep is a missing row)\n")
  else          printf("L0_C3_VERDICT perm_jointW is NOT row-indexed by permsave _rep\n")

  njw   = sum(JW[.,1] :< .)
  nex_c = sum(select(JW[.,1], JW[.,1] :< .) :>= JWO[1,1])
  nex_n = sum(JW[.,1] :>= JWO[1,1])
  p_c = (1 + nex_c)/(njw + 1)
  p_n = (1 + nex_n)/(rows(JW) + 1)
  printf("L0_C4 n_joint=%g (completed=%g, reps=%g) nex_correct=%g nex_naive=%g\n", njw, NC, NR, nex_c, nex_n)
  printf("L0_C4 p_correct=%s  p_naive=%s  overcount=%g  e(perm_failed)=%g  reps-n_joint=%g\n",
         strofreal(p_c, "%21x"), strofreal(p_n, "%21x"), nex_n - nex_c, NF, NR - njw)
  if (nex_n - nex_c == NF) printf("L0_C4_VERDICT overcount EQUALS e(perm_failed)\n")
  else printf("L0_C4_VERDICT overcount (%g) does NOT equal e(perm_failed) (%g); it equals reps - n_joint = %g\n",
              nex_n - nex_c, NF, NR - njw)
  printf("L0_C4 p_correct=%12.8f p_naive=%12.8f\n", p_c, p_n)
  // C5: does the posted p equal the correct hand computation?
  jn = st_matrixcolstripe("_JW")
  for (k = 1; k <= cols(JW); k++) {
    nk = sum(JW[.,k] :< .)
    ek = sum(select(JW[.,k], JW[.,k] :< .) :>= JWO[1,k])
    pk = (1 + ek)/(nk + 1)
    pk_naive = (1 + sum(JW[.,k] :>= JWO[1,k]))/(rows(JW) + 1)
    printf("L0_C5 test=%s posted=%s hand_correct=%s hand_naive=%s n=%g\n",
           jn[k,2], strofreal(st_numscalar("e(" + jn[k,2] + ")"), "%21x"),
           strofreal(pk, "%21x"), strofreal(pk_naive, "%21x"), nk)
  }
}

// the adversarial half
if (st_local("hasaw") == "1") {
  AW = st_matrix("_AW"); AWO = st_matrix("_AWO")
  an = st_matrixcolstripe("_AW")
  printf("L0_ADV rows(perm_advW)=%g reps=%g\n", rows(AW), NR)
  for (k = 1; k <= cols(AW); k++) {
    nk = sum(AW[.,k] :< .)
    ek = sum(select(AW[.,k], AW[.,k] :< .) :>= AWO[1,k])
    pk = (1 + ek)/(nk + 1)
    pn = (1 + sum(AW[.,k] :>= AWO[1,k]))/(rows(AW) + 1)
    printf("L0_ADV test=%s n=%g hand_correct=%s hand_naive=%s\n",
           an[k,2], nk, strofreal(pk, "%21x"), strofreal(pn, "%21x"))
    printf("L0_ADV %s correct=%12.8f naive=%12.8f overcount=%g\n",
           an[k,2], pk, pn, sum(AW[.,k] :>= AWO[1,k]) - ek)
  }
}
end

* posted adversarial p's by name, for the C5 comparison
qui ereturn list
dis "L0_DONE"
