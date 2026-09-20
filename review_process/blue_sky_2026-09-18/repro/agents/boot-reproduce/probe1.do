* probe1.do -- DIFFERENTIAL TEST of wooldid's two bootstrap resampling legs.
* count (fweight) leg vs physical re-identification leg, SAME rseed.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local GEN "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/recon/gen"

which wooldid

* ---------------- data ----------------
local gseed 111
local gN 60
local gT 8
local gt0 2000
local gonsets 2003 2005
local gysd 0.30
local gnever 0.40
include "`GEN'/gen_core.do"
qui count
dis "NOBS=" r(N)
qui levelsof unit, local(uu)
dis "NUNITS=" wordcount("`uu'")
save base.dta, replace

* ---------------- runs ----------------
local R  30
local SD 4242

cap erase bA.dta
cap erase bB.dta
cap erase bC.dta
cap erase bD.dta
cap erase bE.dta
cap erase bG.dta
cap erase bH.dta

* A: att only, cluster(unit)  -> COUNT leg
use base.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(bA) replace
dis "A_RC=" _rc
mat bA = e(b)
mat VA = e(V)
dis "A_b_att = " %21x bA[1,colnumb(bA,"main_att")]
dis "A_se_att= " %21x sqrt(VA[colnumb(bA,"main_att"),colnumb(bA,"main_att")])
dis "A_colnames: " "`: colfullnames bA'"

* B: att + att_i, cluster(unit) -> REID leg (nesting regime)
use base.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(unit) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(bB) replace
dis "B_RC=" _rc
mat bB = e(b)
mat VB = e(V)
dis "B_b_att = " %21x bB[1,colnumb(bB,"main_att")]
dis "B_se_att= " %21x sqrt(VB[colnumb(bB,"main_att"),colnumb(bB,"main_att")])

* C: att + boottstat, cluster(unit) -> REID leg (boottstat route)
use base.dta, clear
cap noi wooldid ync unit t ttre, att cluster(unit) bootstrap(standard) boottstat reps(`R') rseed(`SD') dots(0) bootsave(bC) replace
dis "C_RC=" _rc

* D: att only, cluster(region) -> COUNT leg, coarser cluster
use base.dta, clear
cap noi wooldid ync unit t ttre, att cluster(region) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(bD) replace
dis "D_RC=" _rc
mat bD = e(b)
mat VD = e(V)
dis "D_se_att= " %21x sqrt(VD[colnumb(bD,"main_att"),colnumb(bD,"main_att")])

* E: att + att_i, cluster(region) -> REID leg, coarser cluster (nesting)
use base.dta, clear
cap noi wooldid ync unit t ttre, att att_i cluster(region) bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(bE) replace
dis "E_RC=" _rc
mat bE = e(b)
mat VE = e(V)
dis "E_se_att= " %21x sqrt(VE[colnumb(bE,"main_att"),colnumb(bE,"main_att")])

* G: att only, NO cluster -> COUNT leg, observation resampling
use base.dta, clear
cap noi wooldid ync unit t ttre, att bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(bG) replace
dis "G_RC=" _rc
mat bG = e(b)
mat VG = e(V)
dis "G_se_att= " %21x sqrt(VG[colnumb(bG,"main_att"),colnumb(bG,"main_att")])

* H: att + att_i, NO cluster -> REID leg, plain bsample on rows
use base.dta, clear
cap noi wooldid ync unit t ttre, att att_i bootstrap(standard) reps(`R') rseed(`SD') dots(0) bootsave(bH) replace
dis "H_RC=" _rc
mat bH = e(b)
mat VH = e(V)
dis "H_se_att= " %21x sqrt(VH[colnumb(bH,"main_att"),colnumb(bH,"main_att")])

* ---------------- comparisons ----------------
cap prog drop cmp2
prog def cmp2
  syntax , f1(string) f2(string) tag(string)
  qui use "`f1'.dta", clear
  qui ds
  local v1 `r(varlist)'
  qui rename * X1_*
  qui rename X1__rep _rep
  tempfile t1
  qui save `t1'
  qui use "`f2'.dta", clear
  qui ds
  local v2 `r(varlist)'
  qui merge 1:1 _rep using `t1', nogen
  local common
  foreach v of local v2 {
    if "`v'" == "_rep" continue
    cap confirm variable X1_`v'
    if _rc == 0 local common `common' `v'
  }
  dis "`tag': common columns = " wordcount("`common'")
  local worst = 0
  local worstv
  foreach v of local common {
    qui gen double __d = reldif(`v', X1_`v')
    qui replace __d = 0 if missing(`v') & missing(X1_`v')
    qui sum __d, meanonly
    if r(max) > `worst' {
      local worst = r(max)
      local worstv `v'
    }
    qui drop __d
  }
  dis "`tag': MAXRELDIF = " %12.4e `worst' "   on " "`worstv'"
  * report the w_ columns separately
  local wc
  foreach v of local common {
    if substr("`v'",1,2) == "w_" local wc `wc' `v'
  }
  local wworst = 0
  foreach v of local wc {
    qui gen double __d = abs(`v' - X1_`v')
    qui sum __d, meanonly
    if r(max) > `wworst' local wworst = r(max)
    qui drop __d
  }
  dis "`tag': W-DRAW columns = " wordcount("`wc'") "   MAXABSDIFF = " %12.4e `wworst'
  * main_att column explicitly
  cap confirm variable main_att
  if _rc == 0 {
    cap confirm variable X1_main_att
    if _rc == 0 {
      qui gen double __d = reldif(main_att, X1_main_att)
      qui sum __d, meanonly
      dis "`tag': main_att draws MAXRELDIF = " %12.4e r(max)
      list _rep main_att X1_main_att in 1/6, noobs
      qui drop __d
    }
  }
end

cap noi cmp2 , f1(bA) f2(bB) tag("A(count,clu=unit) vs B(reid,clu=unit)")
dis "CMPAB_RC=" _rc
cap noi cmp2 , f1(bA) f2(bC) tag("A(count) vs C(reid via boottstat)")
dis "CMPAC_RC=" _rc
cap noi cmp2 , f1(bD) f2(bE) tag("D(count,clu=region) vs E(reid,clu=region)")
dis "CMPDE_RC=" _rc
cap noi cmp2 , f1(bG) f2(bH) tag("G(count,nocluster) vs H(reid,nocluster)")
dis "CMPGH_RC=" _rc

dis "PROBE1_DONE"
