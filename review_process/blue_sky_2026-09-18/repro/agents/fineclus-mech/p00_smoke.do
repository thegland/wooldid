clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid
dis "WHICHDONE"

* ---- tiny firm/industry panel: i = firm, cluster = industry (NESTING) ----
* 40 industries x 5 firms = 200 firms, T = 8
clear
set obs 200
gen long firm = _n
gen long ind  = ceil(_n/5)
gen int onset = .
replace onset = 2004 if mod(firm,4)==1
replace onset = 2006 if mod(firm,4)==2
replace onset = 2008 if mod(firm,4)==3
gen double ai = rnormal()
gen double indshk = 0
bysort ind: replace indshk = rnormal() if _n==1
bysort ind: replace indshk = indshk[1]
expand 8
bysort firm: gen int year = 2001 + _n - 1
gen double yr = rnormal()*0.3
gen byte D = (year >= onset) & !missing(onset)
gen double y = ai + indshk + 0.05*(year-2001) + 1.0*D + yr
compress
save "fpanel.dta", replace

cap noisily wooldid y firm year onset, att_i cluster(ind)
dis "BASE_RC=" _rc
dis "SE_ANALYTIC=" _se[att_i]

cap noisily wooldid y firm year onset, att_i cluster(ind) bootstrap(standard) reps(100) rseed(11) dots(0)
dis "BOOT_RC=" _rc
mat b = e(b)
mat V = e(V)
dis "BOOT_SE=" sqrt(V[1,1])
dis "NCLUST=" e(N_clust)
dis "ALLDONE"
