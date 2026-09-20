* p20: ATTRIBUTION for the two session-setting changes the full snapshot diff found
* on a plain wooldid run:  c(matastrict) off -> on, and c(sortrngstate) advancing.
* Are they wooldid's, or its dependencies' / Stata's own?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/state-leak"
run "common.do"

cap prog drop wdms
prog def wdms
  args tag
  dis "MS|`tag'|matastrict=" "`c(matastrict)'"
  dis "MS|`tag'|matalnum=" "`c(matalnum)'"
  dis "MS|`tag'|sortrng=" substr("`c(sortrngstate)'",1,12)
end

wdms 00_fresh_session

* ---- a plain `sort' with ties
qui gen byte tievar = mod(_n,3)
wdms 01_before_sort
qui sort tievar
wdms 02_after_plain_sort

* ---- reghdfe alone
wdgen
wdms 03_before_reghdfe
qui reghdfe y post, absorb(unit year)
wdms 04_after_reghdfe

* ---- gtools alone
clear all
run "common.do"
cap prog drop wdms
prog def wdms
  args tag
  dis "MS|`tag'|matastrict=" "`c(matastrict)'"
  dis "MS|`tag'|matalnum=" "`c(matalnum)'"
  dis "MS|`tag'|sortrng=" substr("`c(sortrngstate)'",1,12)
end
wdgen
wdms 05_before_gtools
capture noisily gegen double gm = mean(y), by(unit)
wdms 06_after_gtools

* ---- margins alone
clear all
run "common.do"
cap prog drop wdms
prog def wdms
  args tag
  dis "MS|`tag'|matastrict=" "`c(matastrict)'"
  dis "MS|`tag'|matalnum=" "`c(matalnum)'"
  dis "MS|`tag'|sortrng=" substr("`c(sortrngstate)'",1,12)
end
wdgen
wdms 07_before_margins
qui regress y i.post##i.sg
qui margins, dydx(post)
wdms 08_after_margins

* ---- and wooldid, in a session where NONE of the above ran
clear all
run "common.do"
cap prog drop wdms
prog def wdms
  args tag
  dis "MS|`tag'|matastrict=" "`c(matastrict)'"
  dis "MS|`tag'|matalnum=" "`c(matalnum)'"
  dis "MS|`tag'|sortrng=" substr("`c(sortrngstate)'",1,12)
end
wdgen
wdms 09_before_wooldid
capture noisily wooldid y unit year ttre, cluster(state) att
dis "W_RC=" _rc
wdms 10_after_wooldid

* ---- what a user's own Mata code does afterwards
capture noisily mata: zz = 1 + 1
dis "MATA_UNDECLARED_RC=" _rc

dis "P20_DONE"
