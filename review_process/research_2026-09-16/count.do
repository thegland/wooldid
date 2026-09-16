clear all
do "./wooldid_counted.ado"
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
tempfile p
save `p'

prog def rpt
  args tag
  di "RESULT `tag' worker=${WRKN} F13=${F13N} SI8=${SI8N}"
end

global F13N 0
global SI8N 0
global WRKN 0
use `p', clear
cap noi qui wooldid y coh t ttre, att cluster(st)
rpt ANALYTIC_rc`=_rc'

global F13N 0
global SI8N 0
global WRKN 0
use `p', clear
cap noi qui wooldid y coh t ttre, att_it cluster(st) bootstrap(standard) reps(20) rseed(4242)
rpt BOOT20_rc`=_rc'

global F13N 0
global SI8N 0
global WRKN 0
use `p', clear
cap noi qui wooldid y coh t ttre, att cluster(st) permute(sharp) permreps(20) rseed(4242)
rpt PERM20_rc`=_rc'
