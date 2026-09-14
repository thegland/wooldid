* mc_adv_bs12_power.do - POWER of the bootstrapped adversarial statistic on the containment
* geometry. 2026-09-14, author ruling ("measure power first, then decide").
*
* WHY. The 2026-09-14 size ladder found AB12/bt rejecting 0 of 200 at a nominal .05 AND 0 of 200
* at .10, with a 10-bin uniformity chi2 of 63.2, where the analytic arm on the same geometry
* rejects .725. The bootstrap repairs containment over-rejection and overshoots. But the ladder
* is all-null by design, so it establishes only that the test does not reject under the null -
* NOT whether it ever rejects. A conservative test and an inert test look identical there.
*
* WHY THE LADDER CANNOT ANSWER THIS. mc_adv_bs12.do's plant is SLOPE*(relyr+1), which is zero at
* the reference period relyr = -1. That is a line THROUGH the reference period, and BS-12's
* one-parameter null contains every such line, so raising SLOPE is a second null, not power.
* The driver header records this and measures it: e(adversarial_pval_att) is numerically
* unchanged by SLOPE to ~7 significant figures on all three arms. So power needs a deviation
* that is NOT a line through the reference period.
*
* THE PLANT USED HERE: a STEP. Effect 0 in every pre period and TAU in every post period. Under
* esfixedbaseperiod with esrelativeto(-1) the reported profile is then 0 across the pre horizons
* and TAU across the post ones, which no line through the reference period can fit unless TAU is
* 0. The same positive control mc_adv_bs12.do's header already cites for the continuous rung,
* where a step "rejects at p = 2e-17..5e-13".
*
* Geometry, VCE and options are AB12's exactly - cluster(i) with safety(off), the containment
* rung - so the only difference from the size campaign is the plant. TAU 0 reproduces the size
* arm and is the built-in control: run it to confirm this harness agrees with the ladder.
*
* Both arms are recorded per rep: bt (the question) and an (the analytic arm, which rejects
* .725 even under the null, so it is a ceiling rather than a control).
*
* args: 1=ado 2=shard 3=nreps 4=seedbase 5=outdir 6=B 7=TAU (0 reproduces the null arm)
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
local TAU = `7'

local G 30
local U 6
local T 12
local ONSET = ceil(`T'/2)

global BSPW_OUT `"`OUTDIR'/mc_adv_bs12_power_s`SHARD'.csv"'

* mc_adv_bs12.do's mkdgp with the plant replaced by a STEP. Everything else is byte-identical
* so the datasets differ from the size campaign's in the outcome alone.
cap program drop mkdgp
program define mkdgp
  args seed tau G U T onset
  clear
  set seed `seed'
  qui set obs `G'
  gen i = _n
  gen ttre = cond(mod(_n, 2) == 0, `onset', .)
  gen a = rnormal()
  qui expand `U'
  bysort i: gen unit = _n
  gen unitid = (i - 1)*`U' + unit
  gen cl2 = (i - 1)*2 + cond(unit <= ceil(`U'/2), 1, 2)
  gen dose = exp(0.8*runiform() + 0.2*mod(i, 5))
  qui expand `T'
  bysort i unitid: gen t = _n
  gen e = rnormal()
  * STEP: 0 before onset, TAU from onset on. NOT a line through the reference period.
  gen y = a + .05*t + `tau'*(ttre < . & t >= `onset') + e
  drop a e
  sort i unitid t
  ereturn clear
end

cap confirm file "${BSPW_OUT}"
if _rc {
  tempname fh0
  file open `fh0' using "${BSPW_OUT}", write replace
  file write `fh0' "config,arm,shard,rep,seed,tau,rc,secs,p,pjoint,n,ncomp,nfail,advabsent" _n
  file close `fh0'
}

cap program drop emitrow
program define emitrow
  args config arm rep seed tau rc secs
  local FLD p pjoint n ncomp nfail advabsent
  foreach nm of local FLD {
    local `nm' .
  }
  if `rc' == 0 {
    cap local p = e(adversarial_pval_att)
    cap local pjoint = e(joint_horiz)
    cap local n = e(N)
    local advabsent = (`p' >= .)
    if "`arm'" == "bt" {
      cap local ncomp = e(boot_completed)
      cap local nfail = e(boot_failed)
      if "`e(boot_adv_absent)'" != "" local advabsent 1
    }
  }
  local row "`config',`arm',${BSPW_SHARD},`rep',`seed',`tau',`rc',`secs'"
  foreach nm of local FLD {
    local x = strtrim(string(``nm'', "%18.12g"))
    if "`x'" == "" local x .
    local row "`row',`x'"
  }
  tempname fh
  file open `fh' using "${BSPW_OUT}", write append
  file write `fh' `"`row'"' _n
  file close `fh'
end
global BSPW_SHARD = `SHARD'

local ESOPT esprelength(4) espostlength(4) esfixedbaseperiod adversarial jointtests
local VCE cluster(i)

forval rep = 1/`NREPS' {
  local seed  = `BASE' + `SHARD'*10000 + `rep'
  local rseed = `BASE' + 500000 + `SHARD'*10000 + `rep'
  * ---- an: the analytic arm (rejects .725 even under the null - a ceiling, not a control)
  timer clear 91
  timer on 91
  mkdgp `seed' `TAU' `G' `U' `T' `ONSET'
  cap qui wooldid y i t ttre, att `VCE' `ESOPT' safety(off)
  local rc = _rc
  timer off 91
  qui timer list 91
  emitrow AB12 an `rep' `seed' `TAU' `rc' `=r(t91)'
  di "BSPW_PROGRESS rep=`rep' AB12 an rc=`rc' secs=" %6.1f r(t91)
  * ---- bt: the question
  timer clear 91
  timer on 91
  mkdgp `seed' `TAU' `G' `U' `T' `ONSET'
  cap qui wooldid y i t ttre, att `VCE' `ESOPT' safety(off) ///
      bootstrap(standard) reps(`B') boottstat rseed(`=`rseed' + 1000000') dots(0)
  local rc = _rc
  timer off 91
  qui timer list 91
  emitrow AB12 bt `rep' `seed' `TAU' `rc' `=r(t91)'
  di "BSPW_PROGRESS rep=`rep' AB12 bt rc=`rc' secs=" %6.1f r(t91)
}
di "BSPW_SHARD_DONE `SHARD'"
exit, clear
