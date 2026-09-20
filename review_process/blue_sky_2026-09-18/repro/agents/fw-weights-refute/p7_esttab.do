* SELF-CONTAINED. Does esttab's footer key on e(vcetype)? Single-model tables, so a
* mixed-vcetype table cannot confound the answer.
clear all
set more off
set processors 2
set seed 20260917
set sortseed 20260917
set obs 30
gen long state = _n
expand 12
bys state: gen int year = _n
gen double y = rnormal() + 0.2*year
gen byte post = (year>=6)
qui regress y post i.state i.year, cluster(state)
dis "CLU_EVCETYPE=[" e(vcetype) "] CLU_EVCE=[" e(vce) "]"
est store rclu
qui regress y post i.state i.year, robust
dis "ROB_EVCETYPE=[" e(vcetype) "] ROB_EVCE=[" e(vce) "]"
est store rrob
qui regress y post i.state i.year
dis "IID_EVCETYPE=[" e(vcetype) "] IID_EVCE=[" e(vce) "]"
est store riid
dis _n "### esttab rclu ALONE ###"
capture noisily esttab rclu, se nostar keep(post)
dis "A_RC=" _rc
dis _n "### esttab rrob ALONE ###"
capture noisily esttab rrob, se nostar keep(post)
dis "B_RC=" _rc
dis _n "### esttab riid ALONE ###"
capture noisily esttab riid, se nostar keep(post)
dis "C_RC=" _rc
dis _n "### estout rrob, stats(vcetype) ###"
capture noisily estout rrob, cells(b se) stats(vcetype N)
dis "D_RC=" _rc
dis "P7_DONE"
