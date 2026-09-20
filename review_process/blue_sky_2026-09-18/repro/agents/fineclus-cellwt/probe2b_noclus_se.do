* probe2b_noclus_se.do -- is the newly-permitted UNCLUSTERED cell-weighted bootstrap
* SE the right one?
*
* Old ado:6353 refused bootstrap(standard) with att_it/att_i/att_itime and no
* cluster(), error 198: "observation-level resampling cannot re-identify duplicated
* treatment cells".  T2 removed it.  wdb_reid is still set (ado:378), wdb_nest is 0 by
* construction (ado:8987-8990 only computes it when clus is non-empty), so the leg
* runs `bsample' on the observations and sets __wdbnewi = i.
*
* Reference: the Monte Carlo sampling SD of the same estimand on the same DGP, which
* is what a standard error is supposed to estimate.  i = unit with one row per (i,t),
* so every treatment cell holds exactly ONE observation: a cell drawn twice is two
* identical rows whose cell mean is unchanged, and the equal-cell weighting counts it
* once.  Only cell PRESENCE survives.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
set linesize 200
which wooldid

prog def mkd
  args seed
  clear
  qui set obs 40
  gen long unit = _n
  gen int ttre = .
  qui replace ttre = 3 if mod(unit,4)==1
  qui replace ttre = 5 if mod(unit,4)==2
  gen double tau = 0.3 + 0.05*unit if !missing(ttre)
  qui replace tau = 0 if missing(tau)
  qui expand 6
  bysort unit: gen int t = _n
  set seed `seed'
  gen double eps = rnormal(0,0.8)
  sort unit t
  gen double y = tau*(t>=ttre & ttre<.) + 0.03*t + eps
end

tempname fh
file open `fh' using "p2b.csv", write replace
file write `fh' "kind,m,rc,b_att,b_atti,b_attit,se_att,se_atti,se_attit" _n
file close `fh'

prog def onerun
  args fh kind m boot rs
  local rc = .
  if "`boot'" == "0" {
    cap qui wooldid y unit t ttre, att att_i att_it
    local rc = _rc
  }
  else {
    cap qui wooldid y unit t ttre, att att_i att_it bootstrap(standard) reps(299) rseed(`rs') dots(0)
    local rc = _rc
  }
  local ba = .
  local bi = .
  local bt = .
  local sa = .
  local si = .
  local st = .
  if `rc' == 0 {
    mat b = e(b)
    mat V = e(V)
    local ba = b[1,colnumb(b,"main_att")]
    local bi = b[1,colnumb(b,"main_att_i")]
    local bt = b[1,colnumb(b,"main_att_it")]
    if "`boot'" != "0" {
      local sa = sqrt(V[colnumb(b,"main_att"),colnumb(b,"main_att")])
      local si = sqrt(V[colnumb(b,"main_att_i"),colnumb(b,"main_att_i")])
      local st = sqrt(V[colnumb(b,"main_att_it"),colnumb(b,"main_att_it")])
    }
  }
  file open `fh' using "p2b.csv", write append
  file write `fh' "`kind',`m',`rc',`ba',`bi',`bt',`sa',`si',`st'" _n
  file close `fh'
end

* Monte Carlo sampling SD
forval m = 1/120 {
  qui mkd `=600000 + `m''
  onerun `fh' "MC" `m' 0 .
}
dis "MC DONE"

* wooldid's unclustered bootstrap on the first 8 of those same datasets
forval m = 1/8 {
  qui mkd `=600000 + `m''
  onerun `fh' "BOOT" `m' 1 `=9000+`m''
}
dis "BOOT DONE"

import delimited "p2b.csv", clear varnames(1) case(preserve)
foreach v in att atti attit {
  qui sum b_`v' if kind == "MC"
  dis "MC   `v'  n=" r(N) "  mean=" %10.6f r(mean) "  SD=" %10.6f r(sd)
  qui sum se_`v' if kind == "BOOT"
  dis "BOOT `v'  n=" r(N) "  mean se=" %10.6f r(mean)
}

* PRE-REGISTERED READOUT, written before the first number landed (2026-09-17):
*   ratio = mean bootstrap se(att_it) / MC sampling SD of att_it.
*     ratio in [0.85, 1.15]  -> CLEAN.  Report as a negative result: the unclustered
*                               cell-weighted bootstrap the removed r198 guard used to
*                               forbid delivers an SE of the right size.
*     ratio < 0.85           -> UNDERSTATED SE on a newly-permitted path.  High.
*     ratio > 1.15           -> the same composition mechanism as F-E, on cells rather
*                               than cohorts.  Medium.
*   Reported for att (the observation-weighted control, which should be ~1 whatever
*   happens) and att_i alongside.
