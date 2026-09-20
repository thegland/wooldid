*! L2A - LENS2 (realism/impact) on FCW-03.  ANALYTIC arms only: the sampling distribution
*! of att_i under (A1) the published convention where the cohort effects are HELD FIXED
*! across replications, and (A2) the same design where the cohort effects are REDRAWN, i.e.
*! the cohorts themselves are part of the sampling experiment.  A3 repeats A1 at LOW cohort
*! heterogeneity.  The bootstrap arms are in L2B.do.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 20260917
set sortseed 20260917
run "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-cellwt-lens2/dgp.doh"

local OUT "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-cellwt-lens2"

*==================================================================================
* STEP 0 - assert dgpL's FIX arm reproduces the PUBLISHED dgpE bit for bit
*==================================================================================
dgpEchk 70001 1 1
rename y yE
keep unit t yE
tempfile ee
qui save "`ee'", replace
dgpL 70001
qui merge 1:1 unit t using "`ee'", keepusing(yE) gen(_mg)
qui count if _mg != 3
dis "MERGE_UNMATCHED=" r(N)
qui gen double dd = abs(y - yE)
qui su dd
dis "IDENT_MAXABSDIFF=" %21x r(max)
dis "IDENT_N=" r(N)
dis "TRUTH_FIX=" $TRUTH_ATTI

tempname fh
local CSV "`OUT'/L2A.csv"
file open `fh' using "`CSV'", write replace
file write `fh' "arm,m,rc,truth_real,b_att,b_atti,se_att,se_atti,nclust" _n
file close `fh'

cap prog drop doarm
prog def doarm
  syntax , arm(string) csv(string) first(int) last(int) [redraw(int 0) lo(real .4) hi(real 1.2)]
  tempname fh
  forval m = `first'/`last' {
    if `redraw' == 0 dgpL `=70000+`m''            , lo(`lo') hi(`hi')
    else             dgpL `=70000+`m'' , redraw(`=500000+`m'') lo(`lo') hi(`hi')
    local truth = $TRUTH_ATTI
    cap qui wooldid y state t ttre, att att_i cluster(state)
    local rc = _rc
    foreach v in b_att b_atti se_att se_atti nclust {
      local `v' .
    }
    if `rc' == 0 {
      tempname B V
      matrix `B' = e(b)
      matrix `V' = e(V)
      local ca = colnumb(`B', "main_att")
      local ci = colnumb(`B', "main_att_i")
      if `ca' < . {
        local b_att  = `B'[1, `ca']
        local se_att = sqrt(`V'[`ca', `ca'])
      }
      if `ci' < . {
        local b_atti  = `B'[1, `ci']
        local se_atti = sqrt(`V'[`ci', `ci'])
      }
      local nclust = e(N_clust)
    }
    file open `fh' using "`csv'", write append
    file write `fh' "`arm',`m',`rc',`truth',`b_att',`b_atti',`se_att',`se_atti',`nclust'" _n
    file close `fh'
  }
end

timer clear 1
timer on 1
doarm , arm(FIX)    csv("`CSV'") first(1) last(150)
timer off 1
qui timer list 1
dis "ARM_FIX_SECS=" r(t1)

timer clear 2
timer on 2
doarm , arm(SUP)    csv("`CSV'") first(1) last(150) redraw(1)
timer off 2
qui timer list 2
dis "ARM_SUP_SECS=" r(t2)

timer clear 3
timer on 3
doarm , arm(LOWFIX) csv("`CSV'") first(1) last(120) lo(.6) hi(1.0)
timer off 3
qui timer list 3
dis "ARM_LOW_SECS=" r(t3)

*==================================================================================
* READOUT
*==================================================================================
qui import delimited using "`CSV'", clear varnames(1) case(preserve)
qui count
dis "ROWS=" r(N)
levelsof arm, local(arms) clean
foreach a of local arms {
  qui count if arm == "`a'" & rc != 0
  dis "ARM `a' NONZERO_RC=" r(N)
  qui su b_atti if arm == "`a'"
  dis "ARM `a' n=" r(N) " mean b_atti=" %12.0g r(mean) " SD b_atti=" %12.0g r(sd)
  qui su b_att if arm == "`a'"
  dis "ARM `a'          mean b_att =" %12.0g r(mean) " SD b_att =" %12.0g r(sd)
  qui su se_atti if arm == "`a'"
  dis "ARM `a'          mean se_atti(analytic)=" %12.0g r(mean)
  qui su se_att if arm == "`a'"
  dis "ARM `a'          mean se_att (analytic)=" %12.0g r(mean)
  qui su truth_real if arm == "`a'"
  dis "ARM `a'          mean realized truth=" %12.0g r(mean) " SD realized truth=" %12.0g r(sd)
  * analytic coverage against .8
  qui gen byte cv_`a' = (b_atti - 1.96*se_atti <= .8) & (b_atti + 1.96*se_atti >= .8) if arm == "`a'"
  qui su cv_`a'
  dis "ARM `a'          analytic normal-CI coverage of .8 = " %6.4f r(mean) "  (n=" r(N) ")"
}
dis "===== L2A DONE"
