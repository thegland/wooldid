* probe9: the aw==duplication identity on three surfaces probe2 could not reach.
*   L: continuous treatment with a GENERATED LATTICE.  The lattice cut points are
*      WEIGHTED quantiles (sthlp:922); weighted quantiles are the classic place where
*      aweight and fweight semantics diverge, so [aw=k] must still cut where the
*      physically duplicated frame cuts.
*   F: cluster() FINER than i under bootstrap(standard) with a cell-weighted estimand -
*      the two-regime unit rule shipped 2026-09-16 and not yet assessed.
*   O: permute(shuffle) + permonsetcohorts - the finer-than-i permutation machinery,
*      also shipped 2026-09-16.
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

mata:
void wd_rel(string scalar a, string scalar b, string scalar tag) {
  real matrix A, B, D
  A = st_matrix(a); B = st_matrix(b)
  if (rows(A)==0 | rows(B)==0) {
    printf("REL_%s=EMPTY\n", tag)
  }
  else if (rows(A)!=rows(B) | cols(A)!=cols(B)) {
    printf("REL_%s=SHAPEMISMATCH %g x %g vs %g x %g\n", tag, rows(A),cols(A),rows(B),cols(B))
  }
  else {
    D = abs(A:-B) :/ (abs(A):*(abs(A):>1e-12) :+ (abs(A):<=1e-12))
    printf("REL_%s=%10.3e  maxabs=%10.3e\n", tag, colmax(rowmax(D)), colmax(rowmax(abs(A:-B))))
  }
}
end

program define pairrun
  syntax , fa(string) ca(string asis) fb(string) cb(string asis) tag(string)
  tempname b1 b2 m1 m2
  qui use "`fa'", clear
  cap noi `ca'
  local r1 = _rc
  if `r1' == 0 {
    matrix `b1' = e(b)
    cap matrix `m1' = e(wdidmainresults)
  }
  qui use "`fb'", clear
  cap noi `cb'
  local r2 = _rc
  if `r2' == 0 {
    matrix `b2' = e(b)
    cap matrix `m2' = e(wdidmainresults)
  }
  dis "PR_`tag'_RC=" `r1' "," `r2'
  if `r1' == 0 & `r2' == 0 {
    mata: wd_rel("`b1'","`b2'","`tag'_b")
    mata: wd_rel("`m1'","`m2'","`tag'_main")
  }
end

* continuous dose with genuine within-cell variation so the lattice has work to do
use "p2_A.dta", clear
set seed 71717
gen double dose2 = post*(runiform()*3)
save "p9_A.dta", replace
expand k
drop k
save "p9_B.dta", replace

*--- L: continuous treatment, with and without a lattice ------------------------
pairrun , fa(p9_A.dta) ca(wooldid y i t ttre [aw=k], att contreat(dose2) contreatcontrols(basic)) fb(p9_B.dta) cb(wooldid y i t ttre, att contreat(dose2) contreatcontrols(basic)) tag(CONT)
pairrun , fa(p9_A.dta) ca(wooldid y i t ttre [aw=k], att contreat(dose2) contreatcontrols(basic) lattice(all 3)) fb(p9_B.dta) cb(wooldid y i t ttre, att contreat(dose2) contreatcontrols(basic) lattice(all 3)) tag(LATALL)
pairrun , fa(p9_A.dta) ca(wooldid y i t ttre [aw=k], att contreat(dose2) contreatcontrols(basic) lattice(baseTreated 3)) fb(p9_B.dta) cb(wooldid y i t ttre, att contreat(dose2) contreatcontrols(basic) lattice(baseTreated 3)) tag(LATBT)
pairrun , fa(p9_A.dta) ca(wooldid y i t ttre [aw=k], att contreat(dose2) contreatcontrols(basic) lattice(all 3) latticeignoreweights) fb(p9_B.dta) cb(wooldid y i t ttre, att contreat(dose2) contreatcontrols(basic) lattice(all 3) latticeignoreweights) tag(LATIGN)
pairrun , fa(p9_A.dta) ca(wooldid y i t ttre [aw=k], att att_it contreat(dose2) contreatcontrols(basic) esprelength(2) espostlength(2)) fb(p9_B.dta) cb(wooldid y i t ttre, att att_it contreat(dose2) contreatcontrols(basic) esprelength(2) espostlength(2)) tag(CONTES)

*--- F: finer-than-i cluster under bootstrap with a cell-weighted estimand ------
pairrun , fa(p9_A.dta) ca(wooldid y i t ttre [aw=k], att_i cluster(fineclus) bootstrap(standard) reps(60) rseed(31415)) fb(p9_B.dta) cb(wooldid y i t ttre, att_i cluster(fineclus) bootstrap(standard) reps(60) rseed(31415)) tag(FINEBOOT)
pairrun , fa(p9_A.dta) ca(wooldid y i t ttre [aw=k], att cluster(fineclus) bootstrap(standard) reps(60) rseed(31415) boottstat) fb(p9_B.dta) cb(wooldid y i t ttre, att cluster(fineclus) bootstrap(standard) reps(60) rseed(31415) boottstat) tag(FINEBTST)

*--- O: permonsetcohorts -------------------------------------------------------
pairrun , fa(p9_A.dta) ca(wooldid y i t ttre [aw=k], att cluster(fineclus) permute(shuffle) permreps(60) rseed(2718) permonsetcohorts) fb(p9_B.dta) cb(wooldid y i t ttre, att cluster(fineclus) permute(shuffle) permreps(60) rseed(2718) permonsetcohorts) tag(PERMONSET)

dis "PROBE9_DONE"
