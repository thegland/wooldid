* mc_fineclus.do - prereg F3 (P4/P5): coverage of the att_i bootstrap interval on a
* DGP with genuinely heterogeneous cohort effects, clustered STRICTLY FINER than i
* (P4) and, as the regression leg, at i (P5).  Shape follows mc_bootse2.do: one CSV
* row per MC draw, appended as it completes, so a run may be killed at any moment and
* analysed from the rows that landed, or resumed by raising <last>.
*
* usage: <stata> -e do mc_fineclus.do <ado> <cfg> <first> <last> <outdir> [reps]
*   cfg FINE = cluster(unit), strictly finer than i=state   (P4)
*   cfg NEST = cluster(state), i.e. cluster == i            (P5)
*
* M AND reps, FIXED BEFORE ANY DRAW EXISTS (2026-09-16): M = 150 per arm, reps(199),
* run as four shards FINE:1:75 FINE:76:150 NEST:1:75 NEST:76:150. Both arms carry the
* SAME M and the same reps or P5 stops being a control for P4. M was chosen off the
* measured per-draw cost (48.9 s FINE, 44.3 s NEST at reps 199 on this design) against
* the inherited band convention, which at n = 150 gives Binomial(150, .95) -> [137, 148];
* it is NOT to be raised after reading the analyzer, since a resumed shard on the same
* seed sequence is indistinguishable in the log from running until something passes.
* reps(199) rather than 99 so the .95 percentile interval is the 5th and 195th order
* statistics rather than an interpolation between the 2nd and 3rd.
*
* MEASURED 2026-09-16, M = 150 per arm, reps(199), 150/150 rows rc 0 and 199/199
* replications complete on EVERY draw of both arms:
*   P4  FINE  cluster(unit), strictly finer than i   att_i 142/150 = .9467  PASS
*   P5  NEST  cluster(state) == i                    att_i 149/150 = .9933  FAIL
* judged against Binomial(150, .95) -> [137, 148]. Reported beside them, never judged:
* att covered 147/150 = .9800 (FINE) and 149/150 = .9933 (NEST), so the DGP's cohort
* heterogeneity is live on both.
*
* P5 fails by OVER-covering, and it is not this feature's doing. cluster(state) == i is
* the nesting regime, which the two-regime branch leaves byte-for-byte unchanged, and
* that was checked rather than asserted: the first 40 draws of the P5 arm were re-run
* against the pristine bd41e01 ado on the same seeds and came back BITWISE identical on
* every field - rc, N_clust, completions, both point estimates and all four intervals.
* So the conservatism is a pre-existing property of the estimator at cluster() == i on
* this design, where each cohort is its own resampling cluster, and it is the same
* effect the help already describes at "What level to cluster at" - at or above the i
* level the covariance stops seeing the treated cells' own variation. Here it lands on
* the conservative side. The feature's own leg, P4, sits on nominal.
clear all
set more off
set linesize 200
set graphics off

local ADO `1'
qui do `"`ADO'"'
run mc_fineclus_dgp.doh

local CFG `2'
if "`CFG'" == "" local CFG FINE
local FIRST `3'
if "`FIRST'" == "" local FIRST 1
local LAST `4'
if "`LAST'" == "" local LAST 50
local OUTDIR `5'
if `"`OUTDIR'"' == "" local OUTDIR .
local REPS `6'
if "`REPS'" == "" local REPS 99

if "`CFG'" == "FINE" local CV unit
else if "`CFG'" == "NEST" local CV state
else {
  dis as error "unknown cfg `CFG'"
  exit 198
}

local CSV "`OUTDIR'/fineclus_`CFG'_`FIRST'_`LAST'.csv"
tempname fh
file open `fh' using "`CSV'", write replace
file write `fh' "cfg,cv,m,seed,rc,secs,nclust,ndone,nfail,b_att,b_atti,pl_att,pu_att,pl_atti,pu_atti,nl_att,nu_att,nl_atti,nu_atti" _n
file close `fh'

* hoisted out of the loop on purpose: a tempname issued inside it would mint a fresh
* matrix every draw and leave <LAST-FIRST+1>*3 of them standing
tempname B P N
forval m = `FIRST'/`LAST' {
  local sd = 70000 + `m'
  qui dgpE `sd' 1 1
  timer clear 3
  timer on 3
  cap qui wooldid y state t ttre, att att_i cluster(`CV') bootstrap(standard) reps(`REPS') rseed(`=900000+`m'') dots(0)
  local rc = _rc
  timer off 3
  qui timer list 3
  local secs = r(t3)

  foreach v in nclust ndone nfail b_att b_atti pl_att pu_att pl_atti pu_atti nl_att nu_att nl_atti nu_atti {
    local `v' .
  }
  if `rc' == 0 {
    local nclust = e(N_clust)
    local ndone  = e(boot_completed)
    local nfail  = e(boot_failed)
    matrix `B' = e(b)
    local ca = colnumb(`B', "main_att")
    local ci = colnumb(`B', "main_att_i")
    if `ca' < . local b_att  = `B'[1, `ca']
    if `ci' < . local b_atti = `B'[1, `ci']
    cap matrix `P' = e(boot_ci_percentile)
    if _rc == 0 {
      if `ca' < . {
        local pl_att = `P'[1, `ca']
        local pu_att = `P'[2, `ca']
      }
      if `ci' < . {
        local pl_atti = `P'[1, `ci']
        local pu_atti = `P'[2, `ci']
      }
    }
    cap matrix `N' = e(boot_ci_normal)
    if _rc == 0 {
      if `ca' < . {
        local nl_att = `N'[1, `ca']
        local nu_att = `N'[2, `ca']
      }
      if `ci' < . {
        local nl_atti = `N'[1, `ci']
        local nu_atti = `N'[2, `ci']
      }
    }
  }
  file open `fh' using "`CSV'", write append
  file write `fh' "`CFG',`CV',`m',`sd',`rc',`secs',`nclust',`ndone',`nfail',`b_att',`b_atti',`pl_att',`pu_att',`pl_atti',`pu_atti',`nl_att',`nu_att',`nl_atti',`nu_atti'" _n
  file close `fh'
}

dis as text "MC_FINECLUS_DONE `CFG' `FIRST' `LAST'"
