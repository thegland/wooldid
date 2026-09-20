clear all
local LV "1 3 7 99"
foreach v in 1 3 7 99 {
  tempname SG`v'
  matrix `SG`v'' = (1,2\3,4)
}
local r 1
foreach nm in att att_it {
  dis ""
  dis "############ estimand: `nm' ############"
  * level estimates and SEs
  local k 0
  foreach v of local LV {
    local k = `k' + 1
    local b`k' = `SG`v''[`r',1]
    local s`k' = `SG`v''[`r',2]
    dis "  level `v': b=" %18.14f `b`k'' " se=" %18.14f `s`k''
  }
  * pairwise: locate each row of e(pwsgs) whose equation is exactly `nm'
  tempname V RR BB
  matrix `V' = J(4,4,.)
  forval a = 1/4 {
    matrix `V'[`a',`a'] = (`s`a'')^2
  }
  dis "MARK OK `nm'"
  local r = `r' + 1
}
