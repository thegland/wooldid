* probe15: WEIGHTED QUANTILES.  A generated lattice cuts the dose at WEIGHTED quantiles
* (sthlp:922).  Weighted quantiles are the classic place where aweight and fweight
* semantics come apart, so this is the sharpest remaining test of the identity.
* probe9's lattice arm died r(322) on both frames (3 rows per cell split into 3 bins gives
* singleton treated cells) - a fixture problem, not a wooldid one.  Here: 6 rows per cell,
* 2 lattice bins, and a dose with wide support.
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
    printf("REL_%s=SHAPEMISMATCH\n", tag)
  }
  else {
    D = abs(A:-B) :/ (abs(A):*(abs(A):>1e-12) :+ (abs(A):<=1e-12))
    printf("REL_%s=%10.3e  maxabs=%10.3e\n", tag, colmax(rowmax(D)), colmax(rowmax(abs(A:-B))))
  }
}
end

clear
set obs 24
gen i = _n
gen ttre = .
replace ttre = 4 if i <= 6
replace ttre = 5 if i > 6  & i <= 12
replace ttre = 6 if i > 12 & i <= 18
gen double ufe = rnormal()
expand 8
bysort i: gen t = _n
expand 6
bysort i t: gen r = _n
gen byte post = (t >= ttre) & !missing(ttre)
gen double dose = post*(0.2 + runiform()*4)
gen double y = ufe + 0.15*t + 0.9*dose + rnormal()*0.5
gen int k = 1 + mod(_n, 4)
save "p15_A.dta", replace
expand k
drop k
save "p15_B.dta", replace

program define pairrun
  syntax , fa(string) ca(string asis) fb(string) cb(string asis) tag(string)
  tempname b1 b2
  qui use "`fa'", clear
  cap noi `ca'
  local r1 = _rc
  if `r1' == 0 matrix `b1' = e(b)
  qui use "`fb'", clear
  cap noi `cb'
  local r2 = _rc
  if `r2' == 0 matrix `b2' = e(b)
  dis "PR_`tag'_RC=" `r1' "," `r2'
  if `r1' == 0 & `r2' == 0 mata: wd_rel("`b1'","`b2'","`tag'_b")
end

pairrun , fa(p15_A.dta) ca(wooldid y i t ttre [aw=k], att contreat(dose) contreatcontrols(basic) nofastcollapse) fb(p15_B.dta) cb(wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) nofastcollapse) tag(NOLAT)
pairrun , fa(p15_A.dta) ca(wooldid y i t ttre [aw=k], att contreat(dose) contreatcontrols(basic) lattice(all 2) nofastcollapse) fb(p15_B.dta) cb(wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) lattice(all 2) nofastcollapse) tag(LAT2)
pairrun , fa(p15_A.dta) ca(wooldid y i t ttre [aw=k], att contreat(dose) contreatcontrols(basic) lattice(all 2) latticeignoreweights nofastcollapse) fb(p15_B.dta) cb(wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) lattice(all 2) latticeignoreweights nofastcollapse) tag(LAT2IGN)
pairrun , fa(p15_A.dta) ca(wooldid y i t ttre [aw=k], att contreat(dose) contreatcontrols(basic) lattice(baseTreated 2) nofastcollapse) fb(p15_B.dta) cb(wooldid y i t ttre, att contreat(dose) contreatcontrols(basic) lattice(baseTreated 2) nofastcollapse) tag(LAT2BT)

* the lattice cut points themselves: a weighted quantile must fall where the duplicated
* frame's unweighted quantile falls
use "p15_A.dta", clear
_pctile dose [aw=k] if dose > 0, nq(2)
dis "AW_MEDIAN=" %21x r(r1)
use "p15_B.dta", clear
_pctile dose if dose > 0, nq(2)
dis "DUP_MEDIAN=" %21x r(r1)
use "p15_A.dta", clear
_pctile dose [fw=k] if dose > 0, nq(2)
dis "FW_MEDIAN=" %21x r(r1)

dis "PROBE15_DONE"
