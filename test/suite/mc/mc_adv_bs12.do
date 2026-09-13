* mc_adv_bs12.do - paired SIZE ladder for the BS-12 adversarial null (Stage C scaffold,
* written 2026-09-13 from the Stage A audit: known_findings KF-2, mc_drivers F6/F8, critic C1).
*
* THE QUESTION. BS-12 (2026-09-10) made the adversarial null the ONE-PARAMETER line through the
* reference period (q = 7 on the 8-horizon window); every size number on record (MC_REPORT
* 2026-08-24, REPORT_ADVBOUND 2026-08-26, sthlp) was measured on the old free-intercept
* statistic, and the ruling's prerequisite simulation was never run. This driver measures the
* size of the new statistic on the ADVBOUND m-ladder - the SAME datasets under three VCE
* geometries - with the inference engines wooldid admits for it:
*   arm an  analytic p     att <vce> ... adversarial jointtests            e(adversarial_pval_att)
*   arm pm  permutation p  att cluster(i) ... permute(shuffle) permreps(B)  (sharp null of no
*           effect at any horizon on the permutation-covariance metric). permute() REQUIRES the
*           cluster() variable to nest i (wooldid_permctl, r(459): the cluster is the unit whose
*           onset is reassigned, and onset is a cohort attribute here), and the p is V-free in
*           fast mode, so there is ONE pm arm per dataset, recorded under config PM.
*   arm bt  boottstat p    + bootstrap(standard) reps(B) boottstat         (the W statistic is
*           bootstrapped; the engine HEAD's PS-64 caution recommends under containment).
*           boottstat with cluster() ALSO requires cluster() to nest i (wooldid_bootctl,
*           r(459): each replication is estimated on a physical expansion of the resampled
*           clusters), and bootstrap(standard) WITHOUT boottstat is admissible below i but does
*           not bootstrap the adversarial statistic (e(boot_adv_absent)) - so bt runs only on
*           rungs whose cluster nests i: AB12. That is exactly the containment cell the caution
*           routes to the bootstrap, so the ladder reads: an on every rung, pm at the cohort
*           level, bt at m = 1.
* Geometries (advb_grid.json, base seed 8200001 -> bit-identical datasets across VCE):
*   AB1  cluster(unitid)  m = 6 clusters per treated cell   (certified regime; an only)
*   AB13 cluster(cl2)     m = 2                             (ladder middle rung; an only)
*   AB12 cluster(i)       m = 1 = CONTAINMENT, safety(off)  (the negative control; an + bt)
* DGP: mc_advb.do's, copied verbatim into mkdgp (G=30 cohorts, U=6 units, T=12, iid errors,
* even cohorts treated at ONSET=6, dose generated but unused) plus the PLANT: SLOPE*(relyr+1) on
* treated units in every period (mc_adv.do A2's form). It is ZERO at the reference period
* relyr = -1, so the true profile is a line THROUGH the reference period = inside the BS-12
* null for every SLOPE; SLOPE = 0 is the flat null (also inside the old null). A nonzero SLOPE
* is the null the old statistic could not distinguish from the new one, so run BOTH.
* Degenerate-replication columns: ncomp/nfail (engine's completed/failed replications),
* advabsent (p missing at rc 0, or e(boot_adv_absent)), advw_n / advw_big / advw_max (the
* replication W statistics in e(boot_advW) / e(perm_advW): count, count > 1e10 (KF-2's escaped
* 3.58e25 draw), maximum), tstud_n (e(boot_tstud) row n for main_att).
*
* One CSV row per (rep, config, arm), appended as it completes (open/write/close).
* args: 1=ado 2=shard 3=nreps 4=baseseed 5=outdir 6=B (permreps/reps; production 199; the
*       engine refuses permreps() < 20) 7=SLOPE (0 = flat null; .15 = mc_adv.do A2's slope)
*       8=geometries (a +-separated list - a comma would be read by `do` as the start of
*       options; default AB1+AB13+AB12)
* Cost measured in the 2-rep smoke at B=20 (2026-09-13): an 1.5-1.9 s, pm 15 s, bt 26 s per rep;
* scaled to B=199: pm ~2.5 min, bt ~4.5 min per rep, so 400 reps are ~47 h single-stream -
* shard 8 ways (~6 h), or run bt at M = 200 with the binomial region widened accordingly.
*
* ACCEPTANCE (to be preregistered before the production run; house M2 convention of
* analyze_mc.py: central ~95% binomial region rounded outward, 10-bin chi2 co-primary):
*   - per (config, arm): rejection COUNT at .05 and .10 inside Binomial(M, alpha)'s region;
*     the pm arm is exact under the sharp null and is the reference the other two are read
*     against on the same datasets (paired);
*   - integrity: rc = 0 on every row, ncomp = B on every pm/bt row, advabsent = 0 everywhere
*     and advw_big = 0 on PM (a nonzero advw_big on AB12/bt is a FINDING about the containment
*     escape route - KF-2's escaped 3.58e25 draw - reported, not scored).
* Interpretation: an on AB12 is expected to FAIL (containment, ruled uncalibrated 2026-08-26;
* the smoke already shows p = .017 / .002 and an e(joint_horiz) of 3.34 - a p above 1 - on
* rep 1 at SLOPE .15); an on AB1 is the re-measurement KF-2 asks for; bt on AB12 is the cell of
* interest.
*
* USAGE (2026-09-13). Positional, shared with mc_bands.do and mc_ri_ext.do; the runner is
* run_stagec.sh, which honours the machine's Stata instance cap:
*     ./run_stagec.sh bs12 <outdir> <nreps> <B> <slope> <geoms|-> <shard>...
* Direct invocation:
*     <stata> -e do mc_adv_bs12.do <ado> <shard> <nreps> <seedbase> <outdir> <B> <slope> <geoms>
* geoms is AB1+AB13+AB12 (mc_advb.do's three VCE geometries on one G=30 U=6 T=12 shape),
* + or , separated. slope 0 is the NULL arm - the one the size ladder needs. B is the
* bootstrap/permutation replication count and wooldid's permreps() floor is 20, so B below 20
* makes the pm arm refuse rather than run. Writes <outdir>/mc_adv_bs12_s<shard>.csv, one row
* per (rep, arm). SMOKE-RUN 2026-09-13 at nreps 1, B 20, slope 0, geoms AB1: rc 0, 2 rows.
* NOT yet campaigned and there is no analyzer: a scaffold, not evidence.
* GAP, recorded so it is not lost: the chunk asks for the adversarial size to be measured
* "including with continuous treatments". The geometries here are all BINARY. A contreat
* geometry (AB1's shape with a unit-level dose varying WITHIN cohort, as in mc_ri_ext.do's E1)
* still has to be added before this driver answers the question as asked.
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local SHARD = `2'
local NREPS = `3'
local BASE  = `4'
local OUTDIR `5'
local B = `6'
local SLOPE = `7'
local GEOMS `8'
if "`GEOMS'" == "" local GEOMS AB1+AB13+AB12
local GEOMS = subinstr(subinstr("`GEOMS'", "+", " ", .), ",", " ", .)

* mc_advb.do's geometry constants for AB1/AB13/AB12
local G 30
local U 6
local T 12
local ONSET = ceil(`T'/2)

global MCBS12_OUT `"`OUTDIR'/mc_adv_bs12_s`SHARD'.csv"'
global MCBS12_SHARD = `SHARD'

*----------------------------------------------------------------- DGP (mc_advb.do's, copied)
cap program drop mkdgp
program define mkdgp
  args seed slope G U T onset
  clear
  set seed `seed'
  qui set obs `G'
  gen i = _n
  * even cohorts treated, odd never-treated (the comparison pool in every period)
  gen ttre = cond(mod(_n, 2) == 0, `onset', .)
  * cohort-level intercept: absorbed by the cohort FE, kept for comparability (inert)
  gen a = rnormal()
  qui expand `U'
  bysort i: gen unit = _n
  gen unitid = (i - 1)*`U' + unit
  gen cl2 = (i - 1)*2 + cond(unit <= ceil(`U'/2), 1, 2)
  * dose varies WITHIN each cohort-period cell (across units) - generated as in mc_advb.do so
  * the RNG stream (and hence the dataset) is identical to the AB campaign's; unused here
  gen dose = exp(0.8*runiform() + 0.2*mod(i, 5))
  qui expand `T'
  bysort i unitid: gen t = _n
  gen e = rnormal()
  * PLANT: a line through the reference period relyr = -1 (value slope*(relyr+1)); with
  * slope = 0 this is mc_advb.do's PWR=0 flat null exactly
  gen y = a + .05*t + `slope'*(t - `onset' + 1)*(ttre < .) + e
  drop a e
  sort i unitid t
  ereturn clear
end

*----------------------------------------------------------------- row writer
cap program drop emitrow
program define emitrow
  args config arm rep seed slope rc secs
  local FLD p pjoint n ncomp nfail advabsent advw_n advw_big advw_max tstud_n
  foreach nm of local FLD {
    local `nm' .
  }
  if `rc' == 0 {
    cap local p = e(adversarial_pval_att)
    cap local pjoint = e(joint_horiz)
    cap local n = e(N)
    local advabsent = (`p' >= .)
    local WACC
    if "`arm'" == "pm" {
      cap local ncomp = e(perm_completed)
      cap local nfail = e(perm_failed)
      local WACC perm_advW
    }
    if "`arm'" == "bt" {
      cap local ncomp = e(boot_completed)
      cap local nfail = e(boot_failed)
      if "`e(boot_adv_absent)'" != "" local advabsent 1
      local WACC boot_advW
      tempname TS
      cap matrix `TS' = e(boot_tstud)
      if _rc == 0 {
        local rn = rownumb(`TS', "n")
        local cn = colnumb(`TS', "main_att")
        if `rn' < . & `cn' < . local tstud_n = `TS'[`rn', `cn']
      }
    }
    if "`WACC'" != "" {
      tempname WA
      cap matrix `WA' = e(`WACC')
      if _rc == 0 {
        local wc = colnumb(`WA', "adversarial_pval_att")
        if `wc' >= . local wc 1
        mata: mcbs12_wstats("`WA'", `wc')
      }
    }
  }
  local row "`config',`arm',${MCBS12_SHARD},`rep',`seed',`slope',`rc',`secs'"
  foreach nm of local FLD {
    local x = strtrim(string(``nm'', "%18.12g"))
    if "`x'" == "" local x .
    local row "`row',`x'"
  }
  tempname fh
  file open `fh' using "${MCBS12_OUT}", write append
  file write `fh' `"`row'"' _n
  file close `fh'
end

mata:
void mcbs12_wstats(string scalar wn, real scalar c)
{
  real colvector v
  v = st_matrix(wn)[., c]
  v = select(v, v :< .)
  st_local("advw_n", strofreal(rows(v)))
  if (rows(v) == 0) return
  st_local("advw_big", strofreal(sum(v :> 1e10)))
  st_local("advw_max", strofreal(max(v), "%21x"))
}
end

*----------------------------------------------------------------- header (once per shard)
cap confirm file "${MCBS12_OUT}"
if _rc {
  tempname fh0
  file open `fh0' using "${MCBS12_OUT}", write replace
  file write `fh0' "config,arm,shard,rep,seed,slope,rc,secs," ///
    "p,pjoint,n,ncomp,nfail,advabsent,advw_n,advw_big,advw_max,tstud_n" _n
  file close `fh0'
}

local ESOPT esprelength(4) espostlength(4) esfixedbaseperiod adversarial jointtests

*----------------------------------------------------------------- the loop
forval rep = 1/`NREPS' {
  * advb_grid.json's seeding (base + shard*10000 + rep) so a rep reproduces an AB dataset
  local seed = `BASE' + `SHARD'*10000 + `rep'
  local rseed = `BASE' + 500000 + `SHARD'*10000 + `rep'
  * ---- pm: permutation, once per dataset at the cohort level (permute() requires
  * cluster() to nest i; the p is V-free, so no VCE ladder exists for this arm)
  timer clear 91
  timer on 91
  mkdgp `seed' `SLOPE' `G' `U' `T' `ONSET'
  cap qui wooldid y i t ttre, att cluster(i) `ESOPT' ///
      permute(shuffle) permreps(`B') rseed(`rseed') dots(0)
  local rc = _rc
  timer off 91
  qui timer list 91
  emitrow PM pm `rep' `seed' `SLOPE' `rc' `=r(t91)'
  di "MCBS12_PROGRESS rep=`rep' PM pm rc=`rc' secs=" %6.1f r(t91)
  foreach geo of local GEOMS {
    local NESTS 0
    if "`geo'" == "AB1"  local VCE cluster(unitid)
    if "`geo'" == "AB13" local VCE cluster(cl2)
    if "`geo'" == "AB12" {
      local VCE cluster(i)
      local NESTS 1
    }
    if !inlist("`geo'", "AB1", "AB13", "AB12") {
      di as error "mc_adv_bs12: unknown geometry `geo'"
      exit 198
    }
    local EXTRA
    if "`geo'" == "AB12" local EXTRA safety(off)
    * ---- an: analytic
    timer clear 91
    timer on 91
    mkdgp `seed' `SLOPE' `G' `U' `T' `ONSET'
    cap qui wooldid y i t ttre, att `VCE' `ESOPT' `EXTRA'
    local rc = _rc
    timer off 91
    qui timer list 91
    emitrow `geo' an `rep' `seed' `SLOPE' `rc' `=r(t91)'
    di "MCBS12_PROGRESS rep=`rep' `geo' an rc=`rc' secs=" %6.1f r(t91)
    * ---- bt: boottstat, only where cluster() nests i (bootctl refuses it below i)
    if !`NESTS' continue
    timer clear 91
    timer on 91
    mkdgp `seed' `SLOPE' `G' `U' `T' `ONSET'
    cap qui wooldid y i t ttre, att `VCE' `ESOPT' `EXTRA' ///
        bootstrap(standard) reps(`B') boottstat rseed(`=`rseed' + 1000000') dots(0)
    local rc = _rc
    timer off 91
    qui timer list 91
    emitrow `geo' bt `rep' `seed' `SLOPE' `rc' `=r(t91)'
    di "MCBS12_PROGRESS rep=`rep' `geo' bt rc=`rc' secs=" %6.1f r(t91)
  }
}
di "MCBS12_SHARD_DONE `SHARD'"
exit, clear
