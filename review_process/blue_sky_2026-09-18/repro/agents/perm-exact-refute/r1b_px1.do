* r1b_px1.do -- SKEPTIC probe of PERMEX-PX-1 (self-contained).
* 5 clusters, onsets 4,4,6,.,. -> M = 5!/(2!1!2!) = 30 enumerated assignments.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

clear
set obs 5
gen long unit = _n
gen double onset = .
replace onset = 4 in 1
replace onset = 4 in 2
replace onset = 6 in 3
gen double ue = rnormal()*0.5
expand 8
bysort unit: gen int year = _n
gen double te = rnormal()*0.3
gen byte post = (year >= onset) & onset < .
gen double y = 1 + 0.15*year + ue + te + 0.9*post
tempfile dat
save `dat', replace

prog def rptsd
  syntax , tag(string) nn(real)
  cap confirm matrix e(perm_stepdown)
  if _rc != 0 {
    dis "`tag'_STEPDOWN_ABSENT"
    exit
  }
  tempname S
  matrix `S' = e(perm_stepdown)
  matrix list `S'
  local C = colsof(`S')
  forvalues c = 1/`C' {
    forvalues r = 2/4 {
      local v = `S'[`r',`c']
      dis "`tag'_cell r=`r' c=`c' p=" %20.17f `v' "  p*n=" %14.10f `v'*`nn' "  p*(n+1)=" %14.10f `v'*(`nn'+1)
    }
  }
end

dis "===== (a) permexact ====="
capture noisily wooldid y unit year onset, att espostlength(3) cluster(unit) permute(shuffle) permexact permstepdown dots(0)
dis "A_RC=" _rc
dis "A_perm_exact=" e(perm_exact)
dis "A_perm_reps=" e(perm_reps)
cap confirm matrix e(perm_n)
if _rc==0 matrix list e(perm_n)
rptsd , tag(A) nn(30)
* how many rows of permdist_b are listwise complete over the ES_att columns?
cap confirm matrix e(permdist_b)
if _rc==0 {
  matrix PB = e(permdist_b)
  dis "A_permdist_rows=" rowsof(PB) " cols=" colsof(PB)
  local pbn : colnames PB
  dis "A_permdist_cols: `pbn'"
}

dis "===== (b) sampled draws, permreps(199) ====="
use `dat', clear
capture noisily wooldid y unit year onset, att espostlength(3) cluster(unit) permute(shuffle) permreps(199) permstepdown dots(0)
dis "B_RC=" _rc
dis "B_perm_exact=" e(perm_exact)
cap confirm matrix e(perm_n)
if _rc==0 matrix list e(perm_n)
rptsd , tag(B) nn(199)
dis "DONE_R1B"
