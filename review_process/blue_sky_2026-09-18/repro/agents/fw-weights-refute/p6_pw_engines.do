* SELF-CONTAINED. Does a pweight survive the engines that treat weights specially?
*  L1/L2: sthlp:268 says pweights satisfy unconditionalse's robust requirement (ado:4676
*         exempts them). Does the unconditionalse leg actually RUN under a bare [pw=]?
*  L3/L4: sthlp:532 says bootstrap() falls back to physical duplication when a user weight
*         is present. Does it run, and does it differ from the unweighted leg sensibly?
*  L5:    permute() with a pweight.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

set obs 24
gen long state = _n
gen ttre = .
replace ttre = 5  if mod(_n,4)==1
replace ttre = 8  if mod(_n,4)==2
replace ttre = 10 if mod(_n,4)==3
gen double as = rnormal()
expand 12
bys state: gen int year = _n
gen double y = as + 0.2*year + rnormal() + 0.4*(year>=ttre & ttre<.)
gen double svywt = 0.5 + runiform()

prog def rep
  args tag
  tempname B V
  matrix `B' = e(b)
  matrix `V' = e(V)
  dis "R_`tag' b=" %21x (`B'[1,1]) " se=" %21x (sqrt(`V'[1,1])) " N=" e(N) " vce=[" e(vce) "] vcetype=[" e(vcetype) "]"
end

dis _n "=== L1: [pw=svywt], att, unconditionalse, robust ==="
capture noisily wooldid y state year ttre [pw=svywt], att unconditionalse robust
dis "L1_RC=" _rc
if _rc==0 rep L1

dis _n "=== L2: [pw=svywt], att, unconditionalse, NO robust/cluster (relies on ado:4676) ==="
capture noisily wooldid y state year ttre [pw=svywt], att unconditionalse
dis "L2_RC=" _rc
if _rc==0 rep L2

dis _n "=== L2b: [aw=svywt], att, unconditionalse, NO robust/cluster (should be REFUSED) ==="
capture noisily wooldid y state year ttre [aw=svywt], att unconditionalse
dis "L2b_RC=" _rc

dis _n "=== L3: [pw=svywt], att att_it, bootstrap(standard) reps(25) cluster(state) ==="
capture noisily wooldid y state year ttre [pw=svywt], att att_it bootstrap(standard) reps(25) rseed(77) cluster(state)
dis "L3_RC=" _rc
if _rc==0 rep L3

dis _n "=== L4: unweighted, same bootstrap ==="
capture noisily wooldid y state year ttre, att att_it bootstrap(standard) reps(25) rseed(77) cluster(state)
dis "L4_RC=" _rc
if _rc==0 rep L4

dis _n "=== L5: [pw=svywt], permute(shuffle) permreps(25) cluster(state) ==="
capture noisily wooldid y state year ttre [pw=svywt], att permute(shuffle) permreps(25) rseed(77) cluster(state)
dis "L5_RC=" _rc
if _rc==0 rep L5
dis "P6_DONE"
