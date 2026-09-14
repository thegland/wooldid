* mc_ame_cellsize.do - size of the CONTINUOUS-treatment adversarial line test as a function of
* observations per treated cohort-by-period cell. 2026-09-14.
*
* WHY IT EXISTS. wooldid.sthlp:745 ("Size of the adversarial test") states a condition specific
* to the continuous case:
*
*    "Under a continuous treatment the line test additionally requires cells containing many
*     observations - ON THE ORDER OF FIFTY - because each cell's dose slope is otherwise too
*     imprecisely estimated for the test's covariance to be trustworthy."
*
* That claim shipped without ever having been measured, and mc_adv_bs12.do's AB1C rung cannot
* measure it: AB1C is G=30 x U=6 x T=12, so a treated cell holds SIX observations - about an
* order of magnitude BELOW the condition. AB1C is the violating arm. This driver is the
* compliant one. It holds the spec fixed and varies U, the units per cohort, alone; since a
* treated cohort-by-period cell holds exactly U observations, U IS observations per cell.
*
* WHAT IT MEASURES. Per draw: e(adversarial_pval_ame) and e(adversarial_pval_att) from ONE fit,
* so the two families are paired and the att p is a control that isolates the ame variance path.
* The DGP plants no dose effect at all, so the true ame profile is the zero line - inside the
* adversarial null at every horizon - and any rejection is size. The att family's true profile is
* likewise a line through the reference period.
*
* RESULTS OF RECORD (2026-09-14, 100 draws at each U, this spec, slope 0):
*     U =  6   ame 41/102 = .402 at nominal .05   KS from Uniform .520, uniformity REJECTED
*              att  5/102 = .049                  inside band
*     U = 50   ame 11/100 = .110 at nominal .05   KS .109 against a .136 critical, uniform NOT
*                                                 rejected - but .110 is still about twice
*                                                 nominal and outside the Binomial(100,.05)
*                                                 central region
*              att  8/100 = .080                  inside band
* Companion mechanism run (probe shape, 200 draws at U=6): the reported ame standard errors run
* at .800 of the sampling standard deviation, every one of eight horizons below .85, while the
* att family sits at 1.036 with every horizon in [.96, 1.14]. e(V_es_ame) is POSTED on every
* draw - the documented cell-containment screen does not fire on this geometry.
* Arithmetic worth keeping: a uniform 20 percent understatement of the ame diagonal predicts at
* most ~.22 joint rejection at df 8. The observed .402 is well above that, so e(V_es_ame) is
* MIS-SHAPED and not merely scaled down. This driver records the DIAGONAL only, so it
* establishes that the ame variance is understated and NOT by how much in which direction.
*
* CONCLUSION THE HELP NOW CARRIES (commit 6ccceec): the caveat is real and the direction is
* right; fifty observations per cell is a large improvement and not exactness.
*
* args: 1=ado 2=nreps-per-shard 3=outdir 4=U (units per cohort = observations per cell)
*       5=first-shard 6=last-shard
* Writes <outdir>/probe_ame_u<U>.csv, one row per draw. Run U = 6 and U = 50 to reproduce the
* help's numbers; intermediate U values trace the curve.
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local NREPS = `2'
local OUTDIR `3'
local U  = `4'
local S0 = `5'
local S1 = `6'

local G 30
local T 12
local ONSET = ceil(`T'/2)
local BASE 900000
local SLOPE 0

global AMEU_OUT `"`OUTDIR'/probe_ame_u`U'.csv"'

* ON THE DUPLICATED mkdgp. The program below is a THIRD copy of mc_advb.do's DGP, also carried
* by mc_adv_bs12.do and mc_adv_bs12_power.do. The duplication is DELIBERATE and must not be
* "fixed" into a shared include: each driver reproduces standalone from its own file, which is
* the property that lets a result be re-run years later without the rest of the mc/ tree. If you
* change the DGP, change it in all three and say so in each header - they are meant to stay
* byte-identical in the shared part.
*
cap program drop mkdgp
program define mkdgp
  args seed slope G U T onset
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
  gen y = a + .05*t + `slope'*(t - `onset' + 1)*(ttre < .) + e
  drop a e
  sort i unitid t
  by i unitid: gen cdose = 1 + runiform() if _n == 1
  by i unitid: replace cdose = cdose[1]
  ereturn clear
end

cap confirm file "${AMEU_OUT}"
if _rc {
  tempname fh0
  file open `fh0' using "${AMEU_OUT}", write replace
  file write `fh0' "config,arm,U,obspercell,shard,rep,seed,slope,rc,secs,p_ame,p_att,pjoint,n,ameabsent" _n
  file close `fh0'
}

local ESOPT esprelength(4) espostlength(4) esfixedbaseperiod adversarial jointtests

forval S = `S0'/`S1' {
  forval rep = 1/`NREPS' {
    local seed = `BASE' + `S'*10000 + `rep'
    timer clear 91
    timer on 91
    mkdgp `seed' `SLOPE' `G' `U' `T' `ONSET'
    cap qui wooldid y i t ttre, att cluster(unitid) `ESOPT' contreat(cdose) contreatcontrols(basic)
    local rc = _rc
    timer off 91
    qui timer list 91
    local secs = r(t91)
    foreach nm in p_ame p_att pjoint n ameabsent {
      local `nm' .
    }
    if `rc' == 0 {
      cap local p_ame = e(adversarial_pval_ame)
      cap local p_att = e(adversarial_pval_att)
      cap local pjoint = e(joint_horiz)
      cap local n = e(N)
      local ameabsent = (`p_ame' >= .)
    }
    * observations per treated cohort-by-period cell is U by construction
    local row "AB1C,an,`U',`U',`S',`rep',`seed',`SLOPE',`rc',`secs'"
    foreach nm in p_ame p_att pjoint n ameabsent {
      local x = strtrim(string(``nm'', "%18.12g"))
      if "`x'" == "" local x .
      local row "`row',`x'"
    }
    tempname fh
    file open `fh' using "${AMEU_OUT}", write append
    file write `fh' `"`row'"' _n
    file close `fh'
    di "AMEU U=`U' rep=`S'.`rep' rc=`rc' p_ame=`p_ame' secs=" %6.1f `secs'
  }
}
di "AMEU_DONE U=`U'"
exit, clear
