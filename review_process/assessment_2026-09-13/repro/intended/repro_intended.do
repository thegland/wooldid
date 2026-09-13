* repro_intended.do - old-vs-new reproduction for the intended-change / storage-noise /
* rc-by-design proof sections (2026-09-12). ONE build per Stata process:
*   cd <workdir>/old ; stata-mp -e do <repro>/intended/repro_intended.do old <repo> <workdir>
*   cd <workdir>/new ; stata-mp -e do <repro>/intended/repro_intended.do new <repo> <workdir>
* ARGUMENTS (positional, all three required):
*   `1' build   - "old" or "new"; adopath ++ <workdir>/<build>
*   `2' repo    - the wooldid_dev checkout root
*   `3' workdir - scratch directory holding old/ and new/
* Mirrors test/suite/mkrun.py: same seeds/scrub per config, same dumpres format.
clear all
set more off
set maxvar 20000
set graphics off
local build "`1'"
local REPO "`2'"
local W "`3'"
if "`build'" == "" | "`REPO'" == "" | "`W'" == "" {
  di as error "usage: do repro_intended.do <old|new> <repo-root> <workdir>"
  exit 198
}
adopath ++ "`W'/`build'"
cap mkdir "`W'/`build'/dumps"
global DUMPDIR "`W'/`build'/dumps/"
global BUILDLABEL "`build'"
if "`build'" == "old" global BUILDSHA "301b706"
else global BUILDSHA "08ec0b4"
which wooldid
di "REPRO_BUILD=`build' sha=$BUILDSHA"

cap program drop dumpres
program define dumpres
  gettoken cfg 0 : 0
  gettoken rc  0 : 0
  tempname fh
  file open `fh' using `"${DUMPDIR}`cfg'.csv"', write replace
  file write `fh' "rc|`rc'" _n
  file write `fh' "prov|${BUILDLABEL}|${BUILDSHA}" _n
  qui ds
  file write `fh' `"varlist|`r(varlist)'"' _n
  cap qui count if e(sample)
  if _rc == 0 {
    local v : di %21x r(N)
    file write `fh' "esamp|__N|`v'" _n
    foreach var of varlist _all {
      cap qui sum `var' if e(sample), meanonly
      if _rc == 0 {
        local v : di %21x r(sum)
        file write `fh' "esamp|`var'|`v'" _n
      }
    }
  }
  local escalars : e(scalars)
  foreach s of local escalars {
    local v : di %21x e(`s')
    file write `fh' "scalar|`s'|`v'" _n
  }
  local emacros : e(macros)
  foreach m of local emacros {
    file write `fh' `"macro|`m'|`e(`m')'"' _n
  }
  local emats : e(matrices)
  foreach m of local emats {
    tempname M
    cap matrix `M' = e(`m')
    if _rc != 0 continue
    local R = rowsof(`M')
    local C = colsof(`M')
    local rn : rowfullnames `M'
    local cn : colfullnames `M'
    file write `fh' `"matnames|`m'|`R'x`C'|`rn'|`cn'"' _n
    forval r = 1/`R' {
      forval c = 1/`C' {
        local v : di %21x `M'[`r',`c']
        file write `fh' "mat|`m'|`r'|`c'|`v'" _n
      }
    }
    cap matrix drop `M'
  }
  file write `fh' "end|`cfg'" _n
  file close `fh'
end

cap program drop scrub
program define scrub
  set seed 98765
  set sortseed 98765
  ereturn clear
  cap macro drop WOOLDID_LASTRC WOOLDID_CLEANUP
  cap matrix drop _all
  cap scalar drop _all
  cap estimates drop _all
end

* ======================= basic : PS62 + CORE10 =======================
scrub
use "`REPO'/test/data/dta_binary.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) summarystats jointtests safety(off)
local rc = _rc
dumpres basic `rc'
di "REPRO_CORE10 build=`build' rc=`rc' e(cmdline)=[`e(cmdline)']"
di "REPRO_CORE10_DONE"
tempname V
matrix `V' = e(V)
local R = rowsof(`V')
local nm : colfullnames `V'
di "REPRO_PS62 build=`build' e(V) is `R'x`R'; names: `nm'"
di "REPRO_PS62 build=`build' e(wdid_vcvstate)=[`e(wdid_vcvstate)'] e(wdid_vblocks)=[`e(wdid_vblocks)']"
forval r = 1/`R' {
  local v : di %21x `V'[`r',`r']
  local vd : di %18.0g `V'[`r',`r']
  di "REPRO_PS62 build=`build' diag[`r',`r'] = `v'  (`vd')"
}
mata:
V = st_matrix("`V'")
R = rows(V)
off_es = 0; nz_es = 0; off_other = 0; nz_other = 0
for (a = 1; a <= R; a++) {
  for (b = 1; b <= R; b++) {
    if (a == b) continue
    if (a >= 3 & b >= 3) {
      off_es = off_es + abs(V[a, b]); nz_es = nz_es + (V[a, b] != 0)
    }
    else {
      off_other = off_other + abs(V[a, b]); nz_other = nz_other + (V[a, b] != 0)
    }
  }
}
printf("REPRO_PS62 build=%s ES block rows/cols 3-%g: sum|offdiag| = %21.15g, nonzero offdiag cells = %g\n", "`build'", R, off_es, nz_es)
printf("REPRO_PS62 build=%s outside the ES block (main_att / pretreatment_att rows+cols): sum|offdiag| = %21.15g, nonzero = %g\n", "`build'", off_other, nz_other)
end
matrix list `V', format(%9.2e) title(REPRO_PS62 e(V) `build')
di "REPRO_PS62_DONE"

* ======================= cov_icbad / cov_icbad2 : B21b-ic8 =======================
scrub
use "`REPO'/test/suite/data/dta_s_covars.dta", clear
qui gen double q = icbad - float(icbad)
qui sum q
di "REPRO_B21b_ic8 build=`build' cov_icbad  level~1e6: max|icbad-float(icbad)| = " max(abs(r(min)), abs(r(max))) "  (icbad min/max: " %12.0f r(N) " obs)"
drop q
cap noi wooldid y i t ttre, att cluster(i) interactivecontrols(icbad)
local rc = _rc
dumpres cov_icbad `rc'
tempname B VV
cap matrix `B' = e(b)
cap matrix `VV' = e(V)
local b21 : di %21x `B'[1,1]
local bd : di %18.0g `B'[1,1]
cap local se = sqrt(`VV'[1,1])
di "REPRO_B21b_ic8 build=`build' cov_icbad  rc=`rc' att b[1,1]=`b21' = `bd'  se=`se'  properties=[`e(properties)']"

scrub
use "`REPO'/test/suite/data/dta_s_covars.dta", clear
qui replace icbad = icbad*100
qui gen double q = icbad - float(icbad)
qui sum q
di "REPRO_B21b_ic8 build=`build' cov_icbad2 level~1e8: max|icbad-float(icbad)| = " max(abs(r(min)), abs(r(max)))
drop q
cap noi wooldid y i t ttre, att cluster(i) interactivecontrols(icbad)
local rc = _rc
dumpres cov_icbad2 `rc'
tempname B2 VV2
cap matrix `B2' = e(b)
cap matrix `VV2' = e(V)
local b21 : di %21x `B2'[1,1]
local bd : di %18.0g `B2'[1,1]
local se .
cap local se = sqrt(`VV2'[1,1])
di "REPRO_B21b_ic8 build=`build' cov_icbad2 rc=`rc' att b[1,1]=`b21' = `bd'  se=`se'  properties=[`e(properties)']"
di "REPRO_B21b_ic8_DONE"

* ======================= pois : E-B12 =======================
scrub
use "`REPO'/test/data/dta_binary.dta", clear
cap noi wooldid ycount i t ttre, att cluster(i) poisson
local rc = _rc
dumpres pois `rc'
cap confirm variable _ppmlhdfe_d
di "REPRO_E_B12 build=`build' rc=`rc' confirm-variable _ppmlhdfe_d after the run: rc=" _rc "   (0 = LEFT in the user's data; 111 = absent)"
di "REPRO_E_B12_DONE"

* ======================= allest_aw : A8-B21 =======================
scrub
use "`REPO'/test/data/dta_binary.dta", clear
cap noi wooldid y i t ttre [aw=w], att att_it att_i att_itime cluster(i) esprelength(2) espostlength(2)
local rc = _rc
dumpres allest_aw `rc'
tempname B3
cap matrix `B3' = e(b)
local nm : colfullnames `B3'
di "REPRO_A8_B21 build=`build' rc=`rc' b names: `nm'"
foreach c in 1 3 4 7 {
  local v : di %21x `B3'[1,`c']
  local vd : di %18.0g `B3'[1,`c']
  di "REPRO_A8_B21 build=`build' b[1,`c'] = `v'  (`vd')"
}
di "REPRO_A8_B21_DONE"

* ======================= impv : A9-svmat =======================
scrub
use "`REPO'/test/data/dta_binary.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) impvalues
local rc = _rc
dumpres impv `rc'
tempname B4 V4
cap matrix `B4' = e(b)
cap matrix `V4' = e(V)
local v : di %21x `B4'[1,1]
local vd : di %18.0g `B4'[1,1]
local w : di %21x `V4'[1,1]
local wd : di %18.0g `V4'[1,1]
di "REPRO_A9_svmat build=`build' rc=`rc' b[1,1] = `v' (`vd')   V[1,1] = `w' (`wd')   e(ll)=" e(ll)
di "REPRO_A9_svmat_DONE"

* ======================= trig_cfx : A6-cfxstripe =======================
scrub
use "`REPO'/test/data/dta_cont.dta", clear
cap noi wooldid y i t ttre, att att_it contreat(intensity) contreatcontrols(basic) lattice(all 3) makecfxplots cfxplottypes(post 2) safety(off) cluster(i) esprelength(2) espostlength(2)
local rc = _rc
dumpres trig_cfx `rc'
tempname C1 C2
cap matrix `C1' = e(cfxatt_cfxpost)
if _rc == 0 {
  local rn : rowfullnames `C1'
  local cn : colfullnames `C1'
  di "REPRO_A6 build=`build' rc=`rc' e(cfxatt_cfxpost) " rowsof(`C1') "x" colsof(`C1') " rownames=[`rn'] colnames=[`cn']"
}
else di "REPRO_A6 build=`build' rc=`rc' e(cfxatt_cfxpost) MISSING"
cap matrix `C2' = e(cfxatt_it_cfxpost)
if _rc == 0 {
  local rn : rowfullnames `C2'
  local cn : colfullnames `C2'
  di "REPRO_A6 build=`build' e(cfxatt_it_cfxpost) " rowsof(`C2') "x" colsof(`C2') " rownames=[`rn'] colnames=[`cn']"
}
di "REPRO_A6_DONE"

* ======================= rc-by-design probes =======================
scrub
use "`REPO'/test/data/dta_binary.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) fastnose
di "REPRO_RC build=`build' scr_fastnose (new-only option) rc=" _rc
scrub
use "`REPO'/test/data/dta_cont.dta", clear
cap noi wooldid y i t ttre, contreat(intensity_uz) contreatcontrols(untreatedZero) safety(off) cluster(i) esprelength(2) espostlength(2) summarystats
di "REPRO_RC build=`build' cont_uz_ss (crash cure B5) rc=" _rc
scrub
use "`REPO'/test/data/dta_cont.dta", clear
cap noi wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) lattice(all 0) cluster(i)
di "REPRO_RC build=`build' scr_lat0 (informative refusal d1f7542) rc=" _rc
scrub
use "`REPO'/test/data/dta_binary.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) impvalues subgroups(sg)
di "REPRO_RC build=`build' scr_impsg (both-error parity) rc=" _rc
di "REPRO_RC_DONE"

di "REPRO_INTENDED_ALL_DONE build=`build'"
exit, clear
