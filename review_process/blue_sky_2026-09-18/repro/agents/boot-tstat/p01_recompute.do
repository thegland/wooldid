clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/boot-tstat"

which wooldid

*---- DGP: 24 states, 5 units per state, 10 periods, staggered onsets by state
clear
set obs 24
gen int state = _n
gen int onset = .
replace onset = 2004 if mod(state,4)==1
replace onset = 2006 if mod(state,4)==2
replace onset = 2008 if mod(state,4)==3
* mod==0 -> never treated
expand 5
bys state: gen int unit = _n + 100*state
gen double a_i = rnormal()
gen double s_e = rnormal()*0.4
expand 10
bys unit: gen int t = 1999 + _n
gen double tshock = 0
bys t (unit): replace tshock = rnormal() if _n==1
bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1.0 + 0.2*(t - onset)
replace tau = 0 if !post
gen double y = a_i + s_e + tshock + tau*post + rnormal()*0.7
gen int ttre = onset

dis "NOBS=" _N

*---- RUN 1: boottstat + bootsave
cap noi wooldid y unit t ttre , cluster(state) bootstrap(standard) reps(199) boottstat rseed(4242) bootsave(d1) replace att att_i att_it esprelength(3) espostlength(3) dots(0)
dis "R1_RC=" _rc
mat MR = e(wdidmainresults)
mat list MR, format(%21x)
dis "LEVEL=" e(level)
dis "BOOTDF=" "`e(bootdf)'"
dis "NCLUS=" e(boot_N_clust)
ereturn list
mat B = e(b)
mat list B, format(%21x)
estimates store R1
