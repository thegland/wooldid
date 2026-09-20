* p21: the Mata / scalar namespace leftovers the dirs-diff found on a PLAIN run
* (Mata globals __00000L, __00000M and `k'; Stata scalars __gtools_gfile_*).
* Who creates them, and does a user's own Mata `k' survive?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"
which wooldid

cap prog drop wdns
prog def wdns
  args tag
  dis "NS|`tag'|--MATA--"
  capture noisily mata: mata describe
  dis "NS|`tag'|--SCALARS--"
  capture noisily scalar dir
  dis "NS|`tag'|--END--"
end

wdgen
wdns 00_fresh

* ---- reghdfe alone
qui reghdfe y post, absorb(unit year)
wdns 01_after_reghdfe

* ---- gtools alone, in a session that has NOT run reghdfe
clear all
run "common.do"
cap prog drop wdns
prog def wdns
  args tag
  dis "NS|`tag'|--MATA--"
  capture noisily mata: mata describe
  dis "NS|`tag'|--SCALARS--"
  capture noisily scalar dir
  dis "NS|`tag'|--END--"
end
wdgen
wdns 02_fresh2
capture noisily gegen double gm = mean(y), by(unit)
wdns 03_after_gtools

* ---- wooldid, in a session that has run neither
clear all
run "common.do"
cap prog drop wdns
prog def wdns
  args tag
  dis "NS|`tag'|--MATA--"
  capture noisily mata: mata describe
  dis "NS|`tag'|--SCALARS--"
  capture noisily scalar dir
  dis "NS|`tag'|--END--"
end
wdgen
wdns 04_fresh3
capture noisily wooldid y unit year ttre, cluster(state) att
dis "W1_RC=" _rc
wdns 05_after_wooldid

* ---- does a USER's Mata k survive a wooldid run?
clear all
run "common.do"
wdgen
mata: k = 12345
mata: myvec = (1, 2, 3)
capture noisily mata: k
capture noisily wooldid y unit year ttre, cluster(state) att
dis "W2_RC=" _rc
dis "USER_K_AFTER:"
capture noisily mata: k
dis "USER_MYVEC_AFTER:"
capture noisily mata: myvec

* ---- does wooldid_Mtmp (ado:6689/8121/8791/8795) survive a gteffects engine run?
clear all
run "common.do"
wdgen
capture noisily wooldid y unit year ttre, cluster(state) att gteffects bootstrap(standard) reps(20) dots(0) rseed(7)
dis "W3_RC=" _rc
dis "MTMP_CHECK:"
capture noisily mata: mata describe wooldid_Mtmp
dis "MTMP_RC=" _rc
capture noisily mata: mata describe

dis "P21_DONE"
