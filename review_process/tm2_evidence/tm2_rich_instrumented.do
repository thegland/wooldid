clear all
set more off
set maxvar 20000
set seed 98765
set sortseed 98765
do "/Users/hegland/.claude/jobs/640b9294/tmp/tm2/wooldid_instr.ado"

cap program drop dumpres
program define dumpres
  args cfg rc
  tempname fh
  file open `fh' using "`cfg'.csv", write replace
  file write `fh' "rc|`rc'" _n
  qui ds
  file write `fh' `"varlist|`r(varlist)'"' _n
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
  file close `fh'
end

clear
set seed 20260824
set obs 15
gen i = _n
gen ttre = cond(i<=12, i+3, .)
expand 8
bysort i: gen unit = _n
gen dose = exp(0.6*runiform() + 0.15*mod(i+unit,6))
expand 16
bysort i unit: gen t = _n
gen sg = mod(i,3)+1
gen post = (t >= ttre & ttre < .)
gen pret = cond(ttre < . & t < ttre, (0.08 + 0.06*sg)*dose*(t-ttre), 0)
gen y = 0.5*dose + 0.6*dose*post + pret + 0.1*i + 0.05*t + rnormal()*0.35
drop post pret
ereturn clear
cap noi wooldid y i t ttre, att att_it contreat(dose) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests
local rc = _rc
dumpres rich_instr `rc'

* record what e() stored for the joint pre-trend p-values (expect "." on stock behavior)
foreach nm in joint_Choriz joint_Choriz_i joint_Choriz_s0 joint_Choriz_s1 joint_Choriz_s2 joint_Choriz_s3 joint_Choriz_i_s0 joint_Choriz_i_s1 joint_Choriz_i_s2 joint_Choriz_i_s3 {
  dis `"TM2ERET `nm' = [`e(`nm')']"'
}
dis "TM2RUN COMPLETE rich_instr rc=`rc'"
