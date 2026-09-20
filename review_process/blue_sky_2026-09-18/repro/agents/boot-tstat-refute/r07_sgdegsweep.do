* r07: does the ruled rank-deficiency guard (wooldid_sgdeg, cut 1-1e-2) let through
*      subgroup-contrast standard errors that are cancellation residue?
*      Oracle = the remedy wooldid's own caution names: bootstrap() WITHOUT boottstat,
*      same cluster().  Reported per cell: rho, analytic contrast se and p, bootstrap
*      contrast se and p, and the ratio.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat-refute"
cd "`WD'"
which wooldid

cap prog drop pwget
prog def pwget
  syntax , tag(string) row(integer)
  tempname P A B
  cap mat `P' = e(pwsgs)
  if _rc != 0 {
    dis "`tag'_NOPW=1"
    exit
  }
  mat `A' = e(wdidsgresults_0)
  mat `B' = e(wdidsgresults_1)
  local afn : rowfullnames `A'
  local eqwant = cond(`row'==1, "main", "pretreatment")
  local sea .
  local seb .
  local k 0
  foreach wn of local afn {
    local ++k
    if "`wn'" == "`eqwant':att" {
      local sea = `A'[`k',2]
      local seb = `B'[`k',2]
    }
  }
  local sed = `P'[`row',2]
  local rho = .
  if `sea'<. & `seb'<. & `sed'<. local rho = (`sea'^2+`seb'^2-`sed'^2)/(2*`sea'*`seb')
  dis "`tag'_R`row' b=" %12.8f `P'[`row',1] " se=" %12.9g `sed' " p=" %12.9f `P'[`row',4] ///
      " rho=" %12.9f `rho' " flagged=" (`rho' > 0.99) " sgdegen=" e(wdid_sgdegen)
end

forval cell = 1/8 {
  local sd = 100000 + 7919*`cell'
  set seed `sd'
  set sortseed `sd'
  local G   = cond(mod(`cell',2)==0, 30, 20)
  local rag = cond(`cell'<=4, 0.08, 0.15)
  local ut  = cond(mod(`cell',4)<2, 8, 10)
  clear
  set obs `G'
  gen int state = _n
  gen int onset = .
  local c1 = round(0.25*`G')
  local c2 = round(0.48*`G')
  local c3 = round(0.70*`G')
  qui replace onset = 2002 if state <= `c1'
  qui replace onset = 2003 if state >  `c1' & state <= `c2'
  qui replace onset = 2004 if state >  `c2' & state <= `c3'
  gen double a_i = rnormal()
  gen int nun = 2 + mod(state,4)
  expand nun
  bys state: gen int unit = _n + 100*state
  gen double u_i = rnormal()*0.6
  expand `ut'
  bys unit: gen int t = 1999 + _n
  gen double tshock = 0
  qui bys t (unit): replace tshock = rnormal() if _n==1
  qui bys t (unit): replace tshock = tshock[1]
  gen byte post = (t >= onset) & !missing(onset)
  gen double tau = 1 + 0.2*(t-onset)
  qui replace tau = 0 if !post
  gen double y = a_i + u_i + tshock + tau*post + rnormal()*0.7
  gen int ttre = onset
  gen double _rr = runiform()
  qui drop if _rr < `rag'
  drop _rr
  gen byte sg = mod(state*7 + 3, 5) <= 1
  qui count
  dis "CELL`cell' G=`G' rag=`rag' T=`ut' NOBS=" r(N)
  qui tab sg
  cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg)
  dis "CELL`cell'_ARC=" _rc
  if _rc == 0 {
    pwget , tag(CELL`cell'_ANA) row(1)
    pwget , tag(CELL`cell'_ANA) row(2)
  }
  cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) subgroups(sg) ///
     bootstrap(standard) reps(299) dots(0)
  dis "CELL`cell'_BRC=" _rc
  if _rc == 0 {
    pwget , tag(CELL`cell'_BOOT) row(1)
    pwget , tag(CELL`cell'_BOOT) row(2)
  }
  cap noi wooldid y state t ttre , robust att esprelength(2) espostlength(2) subgroups(sg)
  dis "CELL`cell'_RRC=" _rc
  if _rc == 0 {
    pwget , tag(CELL`cell'_ROB) row(1)
    pwget , tag(CELL`cell'_ROB) row(2)
  }
}
dis "R07_DONE"
