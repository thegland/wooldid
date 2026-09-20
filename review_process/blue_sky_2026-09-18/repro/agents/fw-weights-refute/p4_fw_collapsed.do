* SELF-CONTAINED.  FW-03: does [aw=n]/[pw=n] on collapsed cells substitute for [fw=n]?
* And WOULD [fw=n] have reproduced the microdata answer in the first place?
* Plus: does esttab actually read e(vcetype)?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* ---- microdata: 30 states x 12 years x ~25 people ----
set obs 30
gen long state = _n
gen ttre = .
replace ttre = 5  if mod(_n,4)==1
replace ttre = 8  if mod(_n,4)==2
replace ttre = 10 if mod(_n,4)==3
gen double as = rnormal()
expand 12
bys state: gen int year = _n
gen int nppl = 15 + int(20*runiform())
expand nppl
gen double y = as + 0.2*year + rnormal() + 0.4*(year>=ttre & ttre<.)
gen byte post = (year>=ttre & ttre<.)
count
local NMICRO = r(N)
dis "NMICRO=" `NMICRO'
save micro.dta, replace

* ---- oracle legs on plain regress ----
qui regress y post i.state i.year, cluster(state)
local b_m_clu = _b[post]
local s_m_clu = _se[post]
local n_m     = e(N)
qui regress y post i.state i.year
local b_m_iid = _b[post]
local s_m_iid = _se[post]

use micro.dta, clear
collapse (mean) ybar=y (count) n=y, by(state year ttre post)
count
dis "NCELL=" r(N)
save cells.dta, replace

qui regress ybar post i.state i.year [fw=n], cluster(state)
local b_fw_clu = _b[post]
local s_fw_clu = _se[post]
local n_fw     = e(N)
qui regress ybar post i.state i.year [fw=n]
local b_fw_iid = _b[post]
local s_fw_iid = _se[post]
qui regress ybar post i.state i.year [aw=n], cluster(state)
local b_aw_clu = _b[post]
local s_aw_clu = _se[post]
local n_aw     = e(N)
qui regress ybar post i.state i.year [aw=n]
local b_aw_iid = _b[post]
local s_aw_iid = _se[post]
qui regress ybar post i.state i.year [pw=n], cluster(state)
local b_pw_clu = _b[post]
local s_pw_clu = _se[post]

dis "ORACLE_micro_clu  b=" %21x `b_m_clu'  " se=" %21x `s_m_clu'  " N=" `n_m'
dis "ORACLE_fwcell_clu b=" %21x `b_fw_clu' " se=" %21x `s_fw_clu' " N=" `n_fw'
dis "ORACLE_awcell_clu b=" %21x `b_aw_clu' " se=" %21x `s_aw_clu' " N=" `n_aw'
dis "ORACLE_pwcell_clu b=" %21x `b_pw_clu' " se=" %21x `s_pw_clu'
dis "ORACLE_micro_iid  b=" %21x `b_m_iid'  " se=" %21x `s_m_iid'
dis "ORACLE_fwcell_iid b=" %21x `b_fw_iid' " se=" %21x `s_fw_iid'
dis "ORACLE_awcell_iid b=" %21x `b_aw_iid' " se=" %21x `s_aw_iid'
dis "RATIO_fw_over_micro_CLU = " %12.9f (`s_fw_clu'/`s_m_clu')
dis "RATIO_aw_over_micro_CLU = " %12.9f (`s_aw_clu'/`s_m_clu')
dis "RATIO_fw_over_micro_IID = " %12.9f (`s_fw_iid'/`s_m_iid')
dis "RATIO_aw_over_micro_IID = " %12.9f (`s_aw_iid'/`s_m_iid')
dis "BRELDIF_fw_clu = " %12.4e (reldif(`b_fw_clu',`b_m_clu'))
dis "BRELDIF_aw_clu = " %12.4e (reldif(`b_aw_clu',`b_m_clu'))

* ---- esttab / e(vcetype) check ----
use micro.dta, clear
qui regress y post i.state i.year, cluster(state)
est store rclu
qui regress y post i.state i.year
est store riid
est restore rclu
dis "RCLU_EVCETYPE=[" e(vcetype) "]"
capture noisily esttab rclu riid, se nostar keep(post)
dis "ESTTAB_RC=" _rc

* ---- wooldid legs ----
prog def grab
  args tag
  tempname B V
  matrix `B' = e(b)
  matrix `V' = e(V)
  foreach c in att att_it att_i att_itime {
    local cn = colnumb(`B', "main_`c'")
    if `cn' < . dis "W_`tag'_`c' b=" %21x (`B'[1,`cn']) " se=" %21x (sqrt(`V'[`cn',`cn'])) " N=" e(N)
  }
end

use micro.dta, clear
dis _n "=== wooldid on MICRODATA, cluster(state) ==="
capture noisily wooldid y state year ttre, att att_it cluster(state)
dis "WM_RC=" _rc
grab MICRO

use cells.dta, clear
dis _n "=== wooldid on CELLS [fw=n] ==="
capture noisily wooldid ybar state year ttre [fw=n], att att_it cluster(state)
dis "WFW_RC=" _rc
dis _n "=== wooldid on CELLS [aw=n], cluster(state) ==="
capture noisily wooldid ybar state year ttre [aw=n], att att_it cluster(state)
dis "WAW_RC=" _rc
grab AWCELL
dis _n "=== wooldid on CELLS [pw=n], cluster(state) ==="
capture noisily wooldid ybar state year ttre [pw=n], att att_it cluster(state)
dis "WPW_RC=" _rc
grab PWCELL
dis _n "=== wooldid on CELLS [aw=n], iid ==="
capture noisily wooldid ybar state year ttre [aw=n], att att_it
dis "WAWI_RC=" _rc
grab AWCELLIID

use micro.dta, clear
dis _n "=== wooldid on MICRODATA, iid ==="
capture noisily wooldid y state year ttre, att att_it
dis "WMI_RC=" _rc
grab MICROIID
dis "P4_DONE"
