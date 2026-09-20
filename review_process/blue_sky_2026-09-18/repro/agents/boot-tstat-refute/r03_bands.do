clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777333
set sortseed 777333
local WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat-refute"
cd "`WD'"
which wooldid

* ---- comparator: band vs pointwise, per horizon, for family `fam'
cap prog drop bandcmp
prog def bandcmp
  syntax , tag(string) [fam(string)]
  if "`fam'" == "" local fam att
  tempname W BD
  cap mat `BD' = e(esband_`fam')
  if _rc != 0 {
    dis "`tag'_NOBAND=1  why=`e(esband_why_`fam')'"
    exit
  }
  mat `W' = e(wdidmainresults)
  local rfn : rowfullnames `W'
  local cn : colnames `BD'
  local m = colsof(`BD')
  local nviol_lo 0
  local nviol_hi 0
  local nviol_any 0
  local nnarrow 0
  dis "`tag'_ENGINE=`e(esband_engine)'  CRIT=" %14.9f e(esband_crit_`fam') "  NBAND=" e(esband_n_`fam') "  DF=" e(boot_df)
  forval j = 1/`m' {
    local h : word `j' of `cn'
    local want "ES_`fam':`h'"
    local rr 0
    local k 0
    foreach r of local rfn {
      local ++k
      if "`r'" == "`want'" local rr = `k'
    }
    if `rr' == 0 {
      dis "`tag'_H`j' `h' NOROW"
      continue
    }
    local b   = `W'[`rr',1]
    local se  = `W'[`rr',2]
    local plo = `W'[`rr',5]
    local phi = `W'[`rr',6]
    local blo = `BD'[1,`j']
    local bhi = `BD'[2,`j']
    local pw  = `phi' - `plo'
    local bw  = `bhi' - `blo'
    local vlo = (`blo' > `plo' + 1e-12)
    local vhi = (`bhi' < `phi' - 1e-12)
    local nar = (`bw' < `pw' - 1e-12)
    local nviol_lo = `nviol_lo' + `vlo'
    local nviol_hi = `nviol_hi' + `vhi'
    local nviol_any = `nviol_any' + (`vlo' | `vhi')
    local nnarrow = `nnarrow' + `nar'
    dis "`tag'_H `h'  b=" %12.8f `b' " se=" %12.8f `se' ///
        "  PW[" %12.8f `plo' "," %12.8f `phi' "] w=" %11.7f `pw' ///
        "  BAND[" %12.8f `blo' "," %12.8f `bhi' "] w=" %11.7f `bw' ///
        "  inside_lo=" `vlo' " inside_hi=" `vhi' " narrower=" `nar'
  }
  dis "`tag'_SUMMARY nhoriz=" `m' " inside_lo=" `nviol_lo' " inside_hi=" `nviol_hi' ///
      " any_inside=" `nviol_any' " STRICTLY_NARROWER=" `nnarrow'
end

*==========================================================
dis "### S1: boottstat + simultaneousbands, SKEWED errors, G=25"
local uG 25
local uT 9
local uSkew 1
local uRagged 1
local uSGmode bal
include "`WD'/gen_unbal.do"
qui count
dis "S1_NOBS=" r(N)
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(2) ///
   bootstrap(standard) reps(399) boottstat simultaneousbands dots(0)
dis "S1_RC=" _rc
bandcmp , tag(S1)

*==========================================================
dis "### S2: SAME data, percentile display (no boottstat) -- calibration"
cap noi wooldid y state t ttre , cluster(state) att esprelength(3) espostlength(2) ///
   bootstrap(standard) reps(399) simultaneousbands dots(0)
dis "S2_RC=" _rc
bandcmp , tag(S2)

*==========================================================
dis "### S3: bootdisplay(wald) + simultaneousbands, G=10 (BT-03)"
local uG 10
local uT 9
local uSkew 0
local uRagged 1
local uSGmode bal
include "`WD'/gen_unbal.do"
qui count
dis "S3_NOBS=" r(N)
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) ///
   bootstrap(standard) reps(999) bootdisplay(wald) simultaneousbands dots(0)
dis "S3_RC=" _rc
dis "S3_bootdf=" e(boot_df) "  invt=" %14.9f invt(e(boot_df),0.975) "  bandcrit=" %14.9f e(esband_crit_att)
bandcmp , tag(S3)

*==========================================================
dis "### S4: same data, bootdisplay(percentile) default -- calibration"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) ///
   bootstrap(standard) reps(999) simultaneousbands dots(0)
dis "S4_RC=" _rc
bandcmp , tag(S4)

*==========================================================
dis "### S5: same data, bootdisplay(wald) bootdf(normal) -- does the crit match z?"
cap noi wooldid y state t ttre , cluster(state) att esprelength(2) espostlength(2) ///
   bootstrap(standard) reps(999) bootdisplay(wald) bootdf(normal) simultaneousbands dots(0)
dis "S5_RC=" _rc
dis "S5_bootdf=" e(boot_df) "  z=" %14.9f invnormal(0.975) "  bandcrit=" %14.9f e(esband_crit_att)
bandcmp , tag(S5)
dis "R03_DONE"
