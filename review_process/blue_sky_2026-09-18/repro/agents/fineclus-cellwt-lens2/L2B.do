*! L2B - LENS2 (realism/impact) on FCW-03.  BOOTSTRAP arms.  For each dataset: the ANALYTIC
*! run and the bootstrap(standard) run at cluster(state) == i, so the two numbers a user
*! would see from the same data sit on the same row.  Rows are appended as they complete.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 20260917
set sortseed 20260917
run "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-cellwt-lens2/dgp.doh"

local OUT "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-cellwt-lens2"
local CSV "`OUT'/L2B.csv"
tempname fh0
file open `fh0' using "`CSV'", write replace
file write `fh0' "arm,m,rc_an,rc_bt,truth_real,b_atti,se_atti_an,se_atti_bt,pl_atti,pu_atti,b_att,se_att_an,se_att_bt,pl_att,pu_att,nclust,ndone,secs" _n
file close `fh0'

cap prog drop bootarm
prog def bootarm
  syntax , arm(string) csv(string) first(int) last(int) reps(int) [lo(real .4) hi(real 1.2)]
  tempname fh B V P
  forval m = `first'/`last' {
    dgpL `=70000+`m'' , lo(`lo') hi(`hi')
    local truth = $TRUTH_ATTI
    foreach v in b_att b_atti se_att_an se_atti_an se_att_bt se_atti_bt pl_att pu_att pl_atti pu_atti nclust ndone {
      local `v' .
    }
    * ---- analytic
    cap qui wooldid y state t ttre, att att_i cluster(state)
    local rc_an = _rc
    if `rc_an' == 0 {
      matrix `B' = e(b)
      matrix `V' = e(V)
      local ca = colnumb(`B', "main_att")
      local ci = colnumb(`B', "main_att_i")
      if `ca' < . local se_att_an  = sqrt(`V'[`ca',`ca'])
      if `ci' < . local se_atti_an = sqrt(`V'[`ci',`ci'])
    }
    * ---- bootstrap
    timer clear 9
    timer on 9
    cap qui wooldid y state t ttre, att att_i cluster(state) bootstrap(standard) reps(`reps') rseed(`=900000+`m'') dots(0)
    local rc_bt = _rc
    timer off 9
    qui timer list 9
    local secs = r(t9)
    if `rc_bt' == 0 {
      matrix `B' = e(b)
      matrix `V' = e(V)
      local ca = colnumb(`B', "main_att")
      local ci = colnumb(`B', "main_att_i")
      if `ca' < . {
        local b_att      = `B'[1,`ca']
        local se_att_bt  = sqrt(`V'[`ca',`ca'])
      }
      if `ci' < . {
        local b_atti     = `B'[1,`ci']
        local se_atti_bt = sqrt(`V'[`ci',`ci'])
      }
      local nclust = e(N_clust)
      local ndone  = e(boot_completed)
      cap matrix `P' = e(boot_ci_percentile)
      if _rc == 0 {
        if `ca' < . {
          local pl_att = `P'[1,`ca']
          local pu_att = `P'[2,`ca']
        }
        if `ci' < . {
          local pl_atti = `P'[1,`ci']
          local pu_atti = `P'[2,`ci']
        }
      }
    }
    file open `fh' using "`csv'", write append
    file write `fh' "`arm',`m',`rc_an',`rc_bt',`truth',`b_atti',`se_atti_an',`se_atti_bt',`pl_atti',`pu_atti',`b_att',`se_att_an',`se_att_bt',`pl_att',`pu_att',`nclust',`ndone',`secs'" _n
    file close `fh'
    dis "ROWDONE arm=`arm' m=`m' rc_an=`rc_an' rc_bt=`rc_bt' secs=`secs' se_an=`se_atti_an' se_bt=`se_atti_bt'"
  }
end

bootarm , arm(HETBOOT) csv("`CSV'") first(1) last(25) reps(99)
bootarm , arm(LOWBOOT) csv("`CSV'") first(1) last(15) reps(99) lo(.6) hi(1.0)

*==================================================================================
* READOUT
*==================================================================================
qui import delimited using "`CSV'", clear varnames(1) case(preserve)
qui count
dis "ROWS=" r(N)
levelsof arm, local(arms) clean
foreach a of local arms {
  qui count if arm == "`a'" & (rc_an != 0 | rc_bt != 0)
  dis "ARM `a' NONZERO_RC=" r(N)
  qui gen double sq_bt = se_atti_bt^2 if arm == "`a'"
  qui su sq_bt
  dis "ARM `a' n=" r(N) " RMS se_atti_boot=" %12.0g sqrt(r(mean))
  qui gen double sq_an = se_atti_an^2 if arm == "`a'"
  qui su sq_an
  dis "ARM `a'          RMS se_atti_analytic=" %12.0g sqrt(r(mean))
  qui gen double rat = se_atti_bt/se_atti_an if arm == "`a'"
  qui su rat
  dis "ARM `a'          boot/analytic se_atti: mean=" %8.4f r(mean) " min=" %8.4f r(min) " max=" %8.4f r(max)
  qui gen double ratt = se_att_bt/se_att_an if arm == "`a'"
  qui su ratt
  dis "ARM `a'          boot/analytic se_att  : mean=" %8.4f r(mean)
  qui gen byte cvp = (pl_atti <= .8) & (pu_atti >= .8) if arm == "`a'"
  qui su cvp
  dis "ARM `a'          boot percentile CI covers .8 : " %6.4f r(mean) " (n=" r(N) ")"
  qui gen byte cva = (b_atti - 1.96*se_atti_an <= .8) & (b_atti + 1.96*se_atti_an >= .8) if arm == "`a'"
  qui su cva
  dis "ARM `a'          analytic CI covers .8        : " %6.4f r(mean) " (n=" r(N) ")"
  qui drop sq_bt sq_an rat ratt cvp cva
}
dis "===== L2B DONE"
