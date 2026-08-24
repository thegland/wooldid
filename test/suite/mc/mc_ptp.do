* mc_ptp.do — pretrendpower power-validation MC driver (PREREG_NEWFEAT_2026-08-24 §M3b).
* Fixed binary design G=40, t=1..9, onset 6 (window pool-free, pretest-verified).
* Arms: pilot/Z (gamma=0), W (gamma = pilot-median mdt_wald), N (gamma = pilot-median
* mdt_nis). The planted trend is linear through the reference period: gamma*(relyr+1),
* applied to treated units in ALL periods.
* args: 1=adopath  2=arm (pilot W N Z)  3=shard  4=nreps  5=baseseed  6=outdir  7=gamma
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local ARM `2'
local SHARD = `3'
local NREPS = `4'
local BASE = `5'
local OUTDIR `6'
local GAMMA = `7'

tempname fh
file open `fh' using `"`OUTDIR'/ptp_`ARM'_s`SHARD'.csv"', write replace
file write `fh' "arm,shard,rep,seed,gamma,rc,prep,maxtpre,mdt_wald,mdt_nis,nissize,slope,slope_se" _n

forval rep = 1/`NREPS' {
  local seed = `BASE' + `SHARD'*10000 + `rep'
  set seed `seed'
  clear
  qui set obs 40
  gen i = _n
  gen ttre = cond(mod(_n, 2) == 0, 6, .)
  gen a = rnormal()
  qui expand 9
  bysort i: gen t = _n
  gen y = a + .05*t + `GAMMA'*(t - 6 + 1)*(ttre < .) + rnormal()
  drop a
  sort i t
  ereturn clear
  cap qui wooldid y i t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod pretrendpower(80)
  local rc = _rc
  foreach nm in prep maxtpre mdtw mdtn nsz slp slpse {
    local `nm' .
  }
  if `rc' == 0 {
    local prep = e(roth_prep_att)
    local nsz = e(roth_nissize_att)
    local slp = e(roth_slope_att)
    local slpse = e(roth_slope_se_att)
    tempname RP B V
    cap matrix `RP' = e(roth_power_att)
    if _rc == 0 {
      local mdtw = `RP'[1, 2]
      local mdtn = `RP'[1, 3]
    }
    * max |t| over the pre-period ES coefficients (the NIS test statistic)
    cap matrix `B' = e(b_es_att)
    cap matrix `V' = e(V_es_att)
    local np = e(numpre_att)
    if `np' < . & `np' >= 1 {
      local maxtpre 0
      forval j = 1/`np' {
        if `V'[`j', `j'] > 0 & `V'[`j', `j'] < . {
          local tj = abs(`B'[1, `j'] / sqrt(`V'[`j', `j']))
          if `tj' > `maxtpre' local maxtpre = `tj'
        }
        else local maxtpre .
      }
    }
    cap matrix drop `RP' `B' `V'
  }
  local row "`ARM',`SHARD',`rep',`seed',`GAMMA',`rc'"
  foreach nm in prep maxtpre mdtw mdtn nsz slp slpse {
    local v : di %18.0g ``nm''
    local row "`row',`=trim("`v'")'"
  }
  file write `fh' "`row'" _n
}
file close `fh'
dis "MC_PTP_SHARD_DONE `ARM' `SHARD'"
exit, clear
