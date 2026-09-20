* probe14: pweights force a robust covariance (sthlp:268) - measured bitwise in probe7 -
* but what does wooldid SAY about it?  Header line, e(vce), e(vcetype), e(wdid_vcespec).
* Also attribute: does reghdfe itself label the pweight run robust?
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

use "p2_A.dta", clear
dis "===== [pw=k], no vce option ====="
wooldid y i t ttre [pw=k], att
dis "PW: e(vce)=[" e(vce) "] e(vcetype)=[" e(vcetype) "] wdid_vcespec=[" e(wdid_vcespec) "]"
dis "PW: se=" %14.10f sqrt(e(V)[1,1])

dis "===== [aw=k], no vce option ====="
wooldid y i t ttre [aw=k], att
dis "AW: e(vce)=[" e(vce) "] e(vcetype)=[" e(vcetype) "] wdid_vcespec=[" e(wdid_vcespec) "]"
dis "AW: se=" %14.10f sqrt(e(V)[1,1])

dis "===== [aw=k], robust ====="
wooldid y i t ttre [aw=k], att robust
dis "AWR: e(vce)=[" e(vce) "] e(vcetype)=[" e(vcetype) "] wdid_vcespec=[" e(wdid_vcespec) "]"
dis "AWR: se=" %14.10f sqrt(e(V)[1,1])

dis "===== [pw=k], robust explicitly ====="
wooldid y i t ttre [pw=k], att robust
dis "PWR: e(vce)=[" e(vce) "] e(vcetype)=[" e(vcetype) "] wdid_vcespec=[" e(wdid_vcespec) "]"
dis "PWR: se=" %14.10f sqrt(e(V)[1,1])

dis "===== attribution: what does reghdfe itself say? ====="
qui reghdfe y post [pw=k], absorb(i t)
dis "REGHDFE_pw: e(vce)=[" e(vce) "] e(vcetype)=[" e(vcetype) "] se=" %14.10f _se[post]
qui reghdfe y post [aw=k], absorb(i t) vce(robust)
dis "REGHDFE_awrob: e(vce)=[" e(vce) "] e(vcetype)=[" e(vcetype) "] se=" %14.10f _se[post]
qui reghdfe y post [aw=k], absorb(i t)
dis "REGHDFE_aw: e(vce)=[" e(vce) "] e(vcetype)=[" e(vcetype) "] se=" %14.10f _se[post]
dis "PROBE14_DONE"
