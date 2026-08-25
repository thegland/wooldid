* follow-up: what does wooldid actually leave in e() on the basic+cluster(i) config?
clear all
set more off
adopath + "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall"
set seed 987654321
set sortseed 987654321
set obs 180
gen i = ceil(_n/6)
gen unitid = 1000 + _n
gen dose = rnormal(1,1)
gen cohfe = rnormal(0,1)
bysort i (unitid): replace cohfe = cohfe[1]
expand 12
bysort unitid: gen t = _n
gen ttre = 6 if mod(i,2)==0
gen y = cohfe + .05*t + rnormal(0,1)

wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) esprelength(4) espostlength(4) esfixedbaseperiod adversarial
ereturn list
