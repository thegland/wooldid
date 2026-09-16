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
tempfile p
save `p'

prog def rpt
  args tag
  di "RESULT `tag' F13n=${F13N} SI8n=${SI8N} F13fired=[${F13FIRE}] SI8fired=[${SI8FIRE}]"
end
prog def zero
  global F13N 0
  global SI8N 0
  global F13FIRE
  global SI8FIRE
end

zero
use `p', clear
cap noi qui wooldid y coh t ttre, att_it cluster(st) bootstrap(standard) reps(40) rseed(4242)
rpt BOOTSTD40_rc`=_rc'

zero
use `p', clear
cap noi qui wooldid y coh t ttre, att cluster(coh) permute(shuffle) permreps(40) rseed(4242)
rpt PERMSHUF40_rc`=_rc'

zero
use `p', clear
cap noi qui wooldid y coh t ttre, att cluster(st) permute(shuffle) permreps(40) permonsetcohorts rseed(4242)
rpt PERMONSET40_rc`=_rc'

zero
use `p', clear
cap noi qui wooldid y coh t ttre, att_it cluster(st) boottstat bootstrap(standard) reps(40) rseed(4242)
rpt BOOTTSTAT40_rc`=_rc'
