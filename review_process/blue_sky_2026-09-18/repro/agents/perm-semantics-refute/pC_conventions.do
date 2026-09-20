*! SKEPTIC probe C -- structural checks on the permute() p-value surface that nobody
*! seeded.  Five independent questions, each with its own sentinel:
*!  C1  multi-onset enumeration count: with onsets {2014, 2016, 2016} on 7 clusters
*!      (4 never treated) the number of DISTINCT assignments is 7!/(1! 2! 4!) = 105.
*!      Does e(perm_npossible) say 105, and does the run complete 105 replications?
*!      (The count and the enumerator are two independent algorithms in the source;
*!       nothing reconciles them - ado:13773 says so.)
*!  C2  e(perm_p_ci) containment on a SAMPLED run: is e(perm_p) inside [lo, hi]?
*!  C3  the help says at most one of permabsp / permrawp / permside() may be given.
*!      Does the guard fire?
*!  C4  e(perm_p) == min(1, 2*min(e(perm_p_lo), e(perm_p_hi))) under the default.
*!  C5  under permexact the Monte-Carlo se / Clopper-Pearson rows must be WITHHELD for
*!      every element whose observed row is a member. Are they missing?
*! Run: /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do pC_conventions.do
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 1234567
set sortseed 1234567
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-semantics-refute"

which wooldid

* ===================== panel: 7 hospitals, 2010-2019 ==========================
clear
set obs 7
gen long hosp = _n
gen double onset = .
replace onset = 2014 if hosp == 2
replace onset = 2016 if hosp == 3
replace onset = 2016 if hosp == 6
gen double a_i = rnormal()*2
expand 10
bysort hosp: gen int year = 2010 + _n - 1
gen byte post = (year >= onset) & !missing(onset)
gen double y = a_i + 0.25*(year-2010) + 1.0*post + rnormal()
save pC_panel.dta, replace

* ---------------------------- C1 + C4 + C5 -----------------------------------
cap noi wooldid y hosp year onset, cluster(hosp) att esprelength(1) espostlength(1) ///
    permute(shuffle) permexact
dis "PC1_RC=" _rc
dis "PC1_NPOSSIBLE=" e(perm_npossible) "  EXPECT 105"
dis "PC1_COMPLETED=" e(perm_completed) "  FAILED=" e(perm_failed) "  REPS=" e(perm_reps)
dis "PC1_EXACTFLAG=" e(perm_exact)
dis "PC1_ADDONE=|`e(perm_addone)'|"
dis "PC1_INCOMPLETE=|`e(perm_incomplete)'|"

matrix _P  = e(perm_p)
matrix _PL = e(perm_p_lo)
matrix _PH = e(perm_p_hi)
matrix _PS = e(perm_p_se)
matrix _PC = e(perm_p_ci)
matrix _NR  = e(perm_n)
mata:
P=st_matrix("_P"); PL=st_matrix("_PL"); PH=st_matrix("_PH")
PS=st_matrix("_PS"); PC=st_matrix("_PC"); NN=st_matrix("_NR")
nm=st_matrixcolstripe("_P")
bad4=0; bad5=0
for (j=1;j<=cols(P);j++) {
  if (P[1,j] >= .) continue
  tgt = min((1, 2*min((PL[1,j], PH[1,j]))))
  if (abs(P[1,j]-tgt) > 1e-15*max((abs(tgt),1))) {
    bad4++
    printf("PC4_FAIL %s p=%s vs 2min=%s\n", nm[j,2], strofreal(P[1,j],"%21x"), strofreal(tgt,"%21x"))
  }
  if (PS[1,j] < . | PC[1,j] < . | PC[2,j] < .) {
    bad5++
    printf("PC5_FAIL %s n=%g se_p=%s cp=[%s,%s] (should all be missing under permexact+member)\n",
           nm[j,2], NN[1,j], strofreal(PS[1,j],"%21x"), strofreal(PC[1,j],"%21x"), strofreal(PC[2,j],"%21x"))
  }
}
printf("PC4_RESULT badcols=%g of %g\n", bad4, cols(P))
printf("PC5_RESULT badcols=%g of %g\n", bad5, cols(P))
end

* ------------------------------- C2 ------------------------------------------
use pC_panel.dta, clear
cap noi wooldid y hosp year onset, cluster(hosp) att esprelength(1) espostlength(1) ///
    permute(shuffle) permreps(200) rseed(4242)
dis "PC2_RC=" _rc
dis "PC2_COMPLETED=" e(perm_completed) " FAILED=" e(perm_failed)
matrix _P2  = e(perm_p)
matrix _PS2 = e(perm_p_se)
matrix _PC2 = e(perm_p_ci)
matrix _PL2 = e(perm_p_lo)
matrix _PH2 = e(perm_p_hi)
mata:
P=st_matrix("_P2"); PC=st_matrix("_PC2"); PS=st_matrix("_PS2")
PL=st_matrix("_PL2"); PH=st_matrix("_PH2")
nm=st_matrixcolstripe("_P2")
bad=0
for (j=1;j<=cols(P);j++) {
  if (P[1,j] >= . | PC[1,j] >= . | PC[2,j] >= .) continue
  if (P[1,j] < PC[1,j] - 1e-12 | P[1,j] > PC[2,j] + 1e-12) {
    bad++
    printf("PC2_FAIL %s p=%12.9f ci=[%12.9f,%12.9f] se=%12.9f\n", nm[j,2], P[1,j], PC[1,j], PC[2,j], PS[1,j])
  }
  else printf("PC2_OK   %s p=%12.9f ci=[%12.9f,%12.9f]\n", nm[j,2], P[1,j], PC[1,j], PC[2,j])
}
printf("PC2_RESULT outside=%g of %g\n", bad, cols(P))
end

* ------------------------------- C3 ------------------------------------------
use pC_panel.dta, clear
cap noi wooldid y hosp year onset, cluster(hosp) att permute(shuffle) permreps(20) permabsp permside(right)
dis "PC3a_RC=" _rc
cap noi wooldid y hosp year onset, cluster(hosp) att permute(shuffle) permreps(20) permabsp permrawp
dis "PC3b_RC=" _rc
cap noi wooldid y hosp year onset, cluster(hosp) att permute(shuffle) permreps(20) permrawp permside(left)
dis "PC3c_RC=" _rc

dis "PC_DONE"
