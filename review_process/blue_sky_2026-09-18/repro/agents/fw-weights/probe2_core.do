* probe2: THE METAMORPHIC IDENTITY.  [aw=k] with integer k must reproduce the dataset
* with each row physically repeated k times, for every POINT ESTIMATE and every
* estimand weight.  (SE / N / df legitimately differ under aweight semantics and are
* reported but not asserted.)
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

*--------------------------------------------------------------------
* generator: 24 units, 8 periods, 3 rows per (i,t); weights vary WITHIN cell
*--------------------------------------------------------------------
clear
set obs 24
gen i = _n
gen ttre = .
replace ttre = 4 if i <= 6
replace ttre = 5 if i > 6  & i <= 12
replace ttre = 6 if i > 12 & i <= 18
gen double ufe = rnormal()
gen byte sg   = mod(i,2)
expand 8
bysort i: gen t = _n
expand 3
bysort i t: gen r = _n
gen byte post = (t >= ttre) & !missing(ttre)
gen double x  = rnormal()
gen double dose = post*(0.5 + 0.5*mod(i,3))
gen double y = ufe + 0.15*t + 1.25*post + 0.3*x + rnormal()*0.5
gen int k = 1 + mod(_n, 4)
gen double one = 1
gen long fineclus = 10*i + r          // finer than i
save "p2_A.dta", replace
expand k
drop k
save "p2_B.dta", replace
use "p2_A.dta", clear
count
dis "NA=" r(N)
use "p2_B.dta", clear
count
dis "NB=" r(N)

*--------------------------------------------------------------------
program define runpair
  syntax , opts(string) tag(string) [wopt(string)]
  if "`wopt'" == "" local wopt "[aw=k]"
  tempname bA bB eA eB gA gB
  local escA 0
  local gtA  0
  local nA   .
  local nB   .
  qui use "p2_A.dta", clear
  cap noi wooldid y i t ttre `wopt' , `opts'
  local rcA = _rc
  if `rcA' == 0 {
    matrix `bA' = e(b)
    local nA = e(N)
    cap matrix `eA' = e(escomp)
    if _rc == 0 {
      local escA 1
    }
    cap matrix `gA' = e(gteffects)
    if _rc == 0 {
      local gtA 1
    }
  }
  qui use "p2_B.dta", clear
  cap noi wooldid y i t ttre , `opts'
  local rcB = _rc
  if `rcB' == 0 {
    matrix `bB' = e(b)
    local nB = e(N)
    cap matrix `eB' = e(escomp)
    cap matrix `gB' = e(gteffects)
  }
  dis "PAIR_`tag'_RC=" `rcA' "," `rcB'
  if `rcA' == 0 & `rcB' == 0 {
    dis "PAIR_`tag'_N=" `nA' "," `nB'
    mata: wd_rel("`bA'","`bB'","`tag'_b")
    if `escA' == 1 {
      mata: wd_rel("`eA'","`eB'","`tag'_escomp")
    }
    if `gtA' == 1 {
      mata: wd_rel("`gA'","`gB'","`tag'_gteff")
    }
  }
end

runpair , opts(att) tag(ATT)
runpair , opts(att att_i att_it att_itime) tag(ALL4)
runpair , opts(att att_i att_it att_itime nofastcollapse) tag(ALL4NFC)
runpair , opts(att att_it esprelength(2) espostlength(3)) tag(ES)
runpair , opts(att att_it esprelength(2) espostlength(3) nofastcollapse) tag(ESNFC)
runpair , opts(att gteffects) tag(GTE)
runpair , opts(att att_i subgroups(sg)) tag(SG)
runpair , opts(att cluster(i)) tag(CLUSI)
runpair , opts(att att_i cluster(fineclus)) tag(CLUSFINE)
runpair , opts(att controls(x)) tag(CTRL)
runpair , opts(att fe(sg)) tag(FE)
runpair , opts(att contreat(dose) contreatcontrols(dose)) tag(CONTREAT)
runpair , opts(att summarystats) tag(SUMSTAT)
runpair , opts(att customweightmultiplier(one)) tag(CWM)

* pweight leg : point estimates only
runpair , opts(att att_i att_it att_itime) tag(PWALL4) wopt([pw=k])

dis "PROBE2_DONE"
