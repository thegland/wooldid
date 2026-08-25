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

use "/Users/hegland/Documents/Research/wooldid_dev/test/data/dta_cont", clear
ereturn clear
qui gen sg = mod(i,2)+1
cap noi wooldid y i t ttre, att att_it contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) subgroups(sg) jointtests
local rc = _rc
dumpres cont_sg_attit_instr `rc'

* record what e() stored for the joint pre-trend p-values (expect "." on stock behavior)
foreach nm in joint_Choriz joint_Choriz_i joint_Choriz_s0 joint_Choriz_s1 joint_Choriz_s2 joint_Choriz_s3 joint_Choriz_i_s0 joint_Choriz_i_s1 joint_Choriz_i_s2 joint_Choriz_i_s3 {
  dis `"TM2ERET `nm' = [`e(`nm')']"'
}
dis "TM2RUN COMPLETE cont_sg_attit_instr rc=`rc'"
