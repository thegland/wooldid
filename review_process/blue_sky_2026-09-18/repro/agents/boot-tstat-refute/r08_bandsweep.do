* r08: how often, on MY data, does the simultaneous band land INSIDE the pointwise
*      interval, and how often is it STRICTLY NARROWER?  Separates the two, because
*      sthlp:528 only claims the band "can never come out narrower".
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat-refute"
cd "`WD'"
which wooldid

cap prog drop bandcmp
prog def bandcmp
  syntax , tag(string) [fam(string)]
  if "`fam'" == "" local fam att
  tempname W BD
  cap mat `BD' = e(esband_`fam')
  if _rc != 0 {
    dis "`tag'_NOBAND=1"
    exit
  }
  mat `W' = e(wdidmainresults)
  local rfn : rowfullnames `W'
  local cn : colnames `BD'
  local m = colsof(`BD')
  local nlo 0
  local nhi 0
  local nnar 0
  forval j = 1/`m' {
    local h : word `j' of `cn'
    local want "ES_`fam':`h'"
    local rr 0
    local k 0
    foreach r of local rfn {
      local ++k
      if "`r'" == "`want'" local rr = `k'
    }
    if `rr' == 0 continue
    local plo = `W'[`rr',5]
    local phi = `W'[`rr',6]
    local blo = `BD'[1,`j']
    local bhi = `BD'[2,`j']
    local vlo = (`blo' > `plo' + 1e-12)
    local vhi = (`bhi' < `phi' - 1e-12)
    local nar = ((`bhi'-`blo') < (`phi'-`plo') - 1e-12)
    local nlo = `nlo' + `vlo'
    local nhi = `nhi' + `vhi'
    local nnar = `nnar' + `nar'
    if `vlo' | `vhi' | `nar' dis "`tag'_VIOL `h' PW[" %11.7f `plo' "," %11.7f `phi' "] BAND[" %11.7f `blo' "," %11.7f `bhi' "] lo=" `vlo' " hi=" `vhi' " narrower=" `nar'
  }
  dis "`tag'_SUM m=" `m' " inside_lo=" `nlo' " inside_hi=" `nhi' " narrower=" `nnar' ///
      " crit=" %12.8f e(esband_crit_`fam') " nband=" e(esband_n_`fam') " df=" e(boot_df) " eng=`e(esband_engine)'"
end

* ---- cells: 1-6 boottstat (BT-01), 7-10 bootdisplay(wald) small G (BT-03)
forval cell = 1/10 {
  local sd = 5000 + 3571*`cell'
  set seed `sd'
  set sortseed `sd'
  local G   = 20
  local skw = 0
  local pre = 3
  local pst = 2
  local rps = 299
  local wald 0
  if `cell'==1 {
    local skw 1
  }
  if `cell'==2 {
    local skw 1
    local G 15
  }
  if `cell'==4 {
    local skw 1
    local pre 2
    local pst 1
  }
  if `cell'==5 {
    local skw 1
    local G 12
    local pre 1
    local pst 1
  }
  if `cell'==6 {
    local skw 1
    local G 30
  }
  if `cell'>=7 {
    local wald 1
    local rps 999
    local skw 0
  }
  if `cell'==7 {
    local G 8
    local pre 1
    local pst 1
  }
  if `cell'==8 {
    local G 10
    local pre 1
    local pst 1
  }
  if `cell'==9 {
    local G 8
    local pre 2
    local pst 2
  }
  if `cell'==10 {
    local G 12
    local pre 1
    local pst 1
  }
  clear
  set obs `G'
  gen int state = _n
  gen int onset = .
  local c1 = max(1,round(0.25*`G'))
  local c2 = round(0.48*`G')
  local c3 = round(0.70*`G')
  qui replace onset = 2002 if state <= `c1'
  qui replace onset = 2003 if state >  `c1' & state <= `c2'
  qui replace onset = 2004 if state >  `c2' & state <= `c3'
  gen double a_i = rnormal()
  gen int nun = 2 + mod(state,3)
  expand nun
  bys state: gen int unit = _n + 100*state
  gen double u_i = rnormal()*0.6
  expand 9
  bys unit: gen int t = 1999 + _n
  gen double tshock = 0
  qui bys t (unit): replace tshock = rnormal() if _n==1
  qui bys t (unit): replace tshock = tshock[1]
  gen byte post = (t >= onset) & !missing(onset)
  gen double tau = 1 + 0.2*(t-onset)
  qui replace tau = 0 if !post
  if `skw' gen double eps = (exp(rnormal()*1.3) - exp(1.3^2/2))*0.7
  else     gen double eps = rnormal()*0.7
  gen double y = a_i + u_i + tshock + tau*post + eps
  gen int ttre = onset
  qui count
  dis "BCELL`cell' G=`G' skew=`skw' pre=`pre' post=`pst' wald=`wald' reps=`rps' NOBS=" r(N)
  if `wald' {
    cap noi wooldid y state t ttre , cluster(state) att esprelength(`pre') espostlength(`pst') ///
       bootstrap(standard) reps(`rps') bootdisplay(wald) simultaneousbands dots(0)
  }
  else {
    cap noi wooldid y state t ttre , cluster(state) att esprelength(`pre') espostlength(`pst') ///
       bootstrap(standard) reps(`rps') boottstat simultaneousbands dots(0)
  }
  dis "BCELL`cell'_RC=" _rc
  if _rc == 0 {
    if `wald' dis "BCELL`cell'_INVT=" %12.8f invt(e(boot_df),0.975)
    bandcmp , tag(BCELL`cell')
  }
}
dis "R08_DONE"
