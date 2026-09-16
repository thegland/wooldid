clear all
do "./w_fire.ado"
clear
set obs 240
gen long unit = _n
gen coh = ceil(unit/40)
gen st  = ceil(unit/8)
gen ttre = cond(coh==6, ., 4 + coh)
expand 16
bysort unit: gen t = _n
set seed 11
gen y = 2*(t>=ttre & ttre<.) + 0.1*unit + 0.05*t + rnormal()
* BREAK IT: give one unit inside cohort 2 a different onset, so ttre is no
* longer constant within i=coh. This is exactly the F13 condition.
replace ttre = 9 if coh==2 & unit==50
tempfile bad
save `bad'

prog def rpt
  args tag
  di "RESULT `tag' F13n=${F13N} F13fired=[${F13FIRE}]"
end
prog def zero
  global F13N 0
  global SI8N 0
  global F13FIRE
  global SI8FIRE
end

* CALIBRATION 1: on broken data the check must FIRE on the observed run
zero
use `bad', clear
cap noi qui wooldid y coh t ttre, att cluster(st)
rpt BADDATA_ANALYTIC_rc`=_rc'

* CALIBRATION 2: and safety(off) must silence it, proving the gate is the gate
zero
use `bad', clear
cap noi qui wooldid y coh t ttre, att cluster(st) safety(off)
rpt BADDATA_SAFETYOFF_rc`=_rc'
