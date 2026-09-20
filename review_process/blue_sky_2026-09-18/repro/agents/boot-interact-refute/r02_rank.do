* r02_rank.do -- SELF-CONTAINED. BOOTX-2 (rank guard) + BOOTX-3 (e(V) rank) at
* REALISTIC replication counts, on my own data.
* Q1: does the analytic ptp pre block come back rank-deficient (caution printed)
*     while the BOOTSTRAP block on the same data comes back "rank k of k" with no
*     caution, at reps() a user would actually run?
* Q2: how big is the pure Monte-Carlo (seed-to-seed) scatter of the printed Roth
*     slope, relative to the se wooldid prints beside it?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

clear
set obs 30
gen state = _n
gen onset = .
replace onset = 7 + mod(_n, 6) if mod(_n,3)==0
expand 2
bysort state: gen cty = _n
gen county = state*10 + cty
expand 16
bysort county: gen t = _n
bysort county (t): gen double ushock = rnormal(0,.3) if _n==1
bysort county (t): replace ushock = ushock[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .8 + .12*rel, 0)
gen double y = 2 + .05*t + ushock + te + e
drop rel te

cap prog drop rpt2
prog def rpt2
  args tag
  tempname b V S bp Si
  matrix `b' = e(b_es_att)
  matrix `V' = e(V_es_att)
  local ln : colnames `b'
  local k 0
  local idx
  local j 0
  foreach l of local ln {
    local ++j
    if substr("`l'",1,3)=="pre" {
      local ++k
      local idx `idx' `j'
    }
  }
  matrix `S' = J(`k',`k',0)
  local r 0
  foreach i1 of local idx {
    local ++r
    local c 0
    foreach i2 of local idx {
      local ++c
      matrix `S'[`r',`c'] = `V'[`i1',`i2']
    }
  }
  matrix `S' = (`S' + `S'')/2
  matrix `Si' = invsym(`S')
  local keffhand = `k' - diag0cnt(`Si')
  mata: st_local("eigs", invtokens(strofreal(sort(symeigenvalues(st_matrix("`S'"))', 1)', "%12.4e")))
  * rank of the full e(V) and of e(V_es_att)
  local rkV = .
  cap confirm matrix e(V)
  if _rc==0 {
    tempname VV
    matrix `VV' = e(V)
    mata: st_local("rkV", strofreal(rank(st_matrix("`VV'"))))
    local cV = colsof(`VV')
  }
  else local cV = .
  mata: st_local("rkVes", strofreal(rank(st_matrix("`V'"))))
  dis "RK|`tag'|k=`k'|keff_hand=`keffhand'|keff_stored=" e(roth_keff_att) ///
      "|rank_Ves=`rkVes'/" rowsof(`V') "|rank_V=`rkV'/`cV'" ///
      "|ncomp=" e(boot_completed) "|NV=" e(boot_N_V) ///
      "|slope=" %12.8g e(roth_slope_att) "|se=" %12.8g e(roth_slope_se_att) ///
      "|ll=" %12.8g e(roth_slope_ll_att) "|ul=" %12.8g e(roth_slope_ul_att) ///
      "|mdt80=" %12.8g e(roth_power_att)[1,2] "|prep=" %12.8g e(roth_prep_att) ///
      "|rankjoint=" %12.8g e(joint_horiz) ///
      "|eig=`eigs'"
end

local OPT "cluster(state) esfixedbaseperiod esprelength(6) espostlength(4) jointtests pretrendpower(80)"

* ---- analytic: does the rank caution fire?
cap noi wooldid y county t onset, `OPT'
dis "A_RC=" _rc
rpt2 ANALYTIC

* ---- bootstrap at realistic reps, several seeds (seed-to-seed = pure Monte-Carlo noise)
foreach R in 50 200 {
  foreach SD in 101 202 303 404 {
    qui cap noi wooldid y county t onset, `OPT' bootstrap(standard) reps(`R') rseed(`SD')
    dis "R`R'S`SD'_RC=" _rc
    rpt2 BOOT`R'_S`SD'
  }
}
dis "DONE"
