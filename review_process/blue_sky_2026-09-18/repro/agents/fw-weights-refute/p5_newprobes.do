* SELF-CONTAINED. Three fresh weight-surface probes nobody suggested:
*  (A) is the aw-vs-fw SE gap on collapsed data ENTIRELY Stata's (N-1)/(N-k) multiplier?
*  (B) sthlp:710 says "Poisson estimation accepts pweights but not aweights (refused with
*      an error)". Is it?
*  (C) zero survey weights: ado:683 drops rows with `w' <= 0. Are they counted/announced?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

set obs 30
gen long state = _n
gen ttre = .
replace ttre = 5  if mod(_n,4)==1
replace ttre = 8  if mod(_n,4)==2
replace ttre = 10 if mod(_n,4)==3
gen double as = rnormal()
expand 12
bys state: gen int year = _n
gen int nppl = 15 + int(20*runiform())
expand nppl
gen double y = as + 0.2*year + rnormal() + 0.4*(year>=ttre & ttre<.)
gen double ycount = rpoisson(exp(0.5 + 0.1*year + 0.3*(year>=ttre & ttre<.)))
gen double svywt = 0.5 + runiform()
save m5.dta, replace

****** (A) df multiplier ******
dis _n "=== A1: wooldid microdata, cluster(state) ==="
capture noisily wooldid y state year ttre, att cluster(state)
dis "A1_RC=" _rc
local seA1 = sqrt(e(V)[1,1])
dis "A1 N=" e(N) " df_r=" e(df_r) " Nclust=" e(N_clust) " se=" %21x `seA1'
scalar sA1 = `seA1'
scalar nA1 = e(N)
scalar dA1 = e(df_r)

use m5.dta, clear
collapse (mean) ybar=y (count) n=y, by(state year ttre)
save c5.dta, replace
dis _n "=== A2: wooldid cells [aw=n], cluster(state) ==="
capture noisily wooldid ybar state year ttre [aw=n], att cluster(state)
dis "A2_RC=" _rc
local seA2 = sqrt(e(V)[1,1])
dis "A2 N=" e(N) " df_r=" e(df_r) " Nclust=" e(N_clust) " se=" %21x `seA2'
scalar sA2 = `seA2'
scalar nA2 = e(N)
scalar dA2 = e(df_r)
dis "SE_RATIO_obs   = " %14.10f (sA2/sA1)
dis "DF_PRED_RATIO  = " %14.10f sqrt( ((nA2-1)/dA2) / ((nA1-1)/dA1) )
dis "RESIDUAL       = " %12.4e (reldif(sA2/sA1, sqrt( ((nA2-1)/dA2) / ((nA1-1)/dA1) )))

****** (B) poisson weight types vs sthlp:710 ******
use m5.dta, clear
dis _n "=== B1: poisson + [pw=svywt] ==="
capture noisily wooldid ycount state year ttre [pw=svywt], att poisson cluster(state)
dis "B1_RC=" _rc
dis _n "=== B2: poisson + [aw=svywt]  (help says REFUSED) ==="
capture noisily wooldid ycount state year ttre [aw=svywt], att poisson cluster(state)
dis "B2_RC=" _rc
dis _n "=== B3: poisson + [aw=svywt] + safety(off) ==="
capture noisily wooldid ycount state year ttre [aw=svywt], att poisson cluster(state) safety(off)
dis "B3_RC=" _rc

****** (C) zero weights ******
use m5.dta, clear
set seed 4242
gen double w0 = svywt
replace w0 = 0 if runiform() < 0.20
count if w0 == 0
local nzero = r(N)
dis "NZERO=" `nzero'
count
dis "NTOT=" r(N)
dis _n "=== C1: [aw=w0] on the full data ==="
capture noisily wooldid y state year ttre [aw=w0], att cluster(state)
dis "C1_RC=" _rc
dis "C1 N=" e(N) " Ndropped=" e(N_dropped) " b=" %21x (e(b)[1,1]) " se=" %21x (sqrt(e(V)[1,1]))
dis _n "=== C2: [aw=w0] if w0>0  (hand-restricted) ==="
capture noisily wooldid y state year ttre if w0>0 [aw=w0], att cluster(state)
dis "C2_RC=" _rc
dis "C2 N=" e(N) " Ndropped=" e(N_dropped) " b=" %21x (e(b)[1,1]) " se=" %21x (sqrt(e(V)[1,1]))
dis _n "=== C3: plain regress comparison ==="
gen byte post = (year>=ttre & ttre<.)
qui regress y post i.state i.year [aw=w0], cluster(state)
dis "C3 N=" e(N)
dis "P5_DONE"
