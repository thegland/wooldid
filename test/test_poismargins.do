* test_poismargins.do - the poisson margins reference-scale gate.
* Run from test/: stata-mp -e do test_poismargins.do
*
* Since the 2026-09-20 removals the gate is UNCONDITIONAL: every margins call wooldid
* makes under poisson is a plain dydx or dyex of the linear index, so both injectors name
* predict(xb) whenever poisson is on and are empty off it. predict() reaches margins'
* ANALYTIC Jacobian, which is exact for a contrast that is linear in b; the
* expression(predict(xb)) branch that used to stand beside it forced a NUMERICAL Jacobian
* and existed only for the ey* operators, which logged the sign-unconstrained linear index.
* Those operators are gone (poisexpresults and semielasticity removed outright,
* contreatelasticitytype narrowed to dydx/dyex), so nothing reaches margins on a non-xb
* scale any more and the branch went with them.
*
* That collapse is numerically inert by construction, not merely by measurement: on the
* surviving path semielasticity defaulted to dydx and BS-6 remapped eydx -> dydx, so the
* expression() branch was already unreachable under poisson-without-poisexpresults. It was
* also checked empirically - pinned c1d8979 against the patched build, in separate Stata
* processes, over four poisson paths (att, att+contreat, att+contreat+dyex, poisson event
* study), b and diagonal V dumped at %21x: bit-identical on every cell.
*
* P1 the xb reporting scale is unchanged (the documented purpose of the gate)
* P2 plain dydx poisson event study runs and produces finite inference
* P3 the surviving contreat operators dydx and dyex both run under poisson on the xb scale
* P4 non-poisson runs are untouched by the gate entirely
* P5 one cluster stops under poisson (rc 459) as it does under reghdfe: the too-few-clusters
*    guard reads e(N_clust) as well as e(df_r), because ppmlhdfe never posts e(df_r)
* P6 under poisson the dose effects are titled for what they are on the linear index: dyex
*    "Elasticities (dyex)", dydx "Semi-elasticities (dydx)" (author ruling 2026-09-23)
*
* The retired option NAMES are pinned in test_si89.do (T1/T2), which is their cert host:
* test/suite/configs_ext.csv cannot carry them, because coverage_check.py treats an option
* token it cannot resolve against the syntax line as FATAL (author ruling 2026-09-11).
* Prints POISMARGINS_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0

use "`REPO'/test/data/dta_binary.dta", clear
tempfile b
save `b'

* ---- P1: the xb scale. The att on the log scale must stay the log-scale value:
* exp() of it is the count-scale ratio, so a scale regression would move it hugely.
qui wooldid ycount i t ttre, att cluster(i) poisson
scalar p1b = _b[main_att]
di "P1_xb_scale_att=" p1b
if abs(p1b) > 5 {
  di as error "P1 FAIL: att is implausible for a log-scale report (scale regression?)"
  local ++fails
}

* ---- P2: plain dydx with an event study - the gated path - runs with finite inference
use `b', clear
cap qui wooldid ycount i t ttre, att cluster(i) poisson esprelength(2) espostlength(2)
local p2rc = _rc
matrix P2 = e(wdidmainresults)
di "P2_dydx_es_rc=" `p2rc' " se=" P2[1,2]
if `p2rc' != 0 | P2[1,2] >= . | P2[1,2] <= 0 {
  di as error "P2 FAIL: gated dydx path rc " `p2rc' " or non-finite SE"
  local ++fails
}

* ---- P3: both surviving contreat operators run under poisson. dyex logs the REGRESSOR,
* not the prediction, so it takes the predict(xb) form too (F032, author ruling 2026-09-06)
* and is not affected by the ey* removal. This leg is what stops a future narrowing of the
* gate from quietly taking dyex with it.
foreach v in dydx dyex {
  use "`REPO'/test/data/dta_cont.dta", clear
  cap qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) ///
      poisson contreatelasticitytype(`v') safety(off)
  local p3rc = _rc
  di "P3_contreat_`v'_rc=" `p3rc'
  if `p3rc' != 0 {
    di as error "P3 FAIL: contreatelasticitytype(`v') under poisson aborted (rc `p3rc')"
    local ++fails
  }
}

* ---- P4: non-poisson is untouched - the gate never fires
use `b', clear
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
scalar p4 = _b[main_att]
di "P4_nonpois_att=" p4
if p4 >= . {
  di as error "P4 FAIL: non-poisson run broken"
  local ++fails
}

* ---- P5: onecluster-guard-poisson. the worker's too-few-clusters guard keyed on e(df_r)
* alone, which ppmlhdfe never posts (missing, and `. < 1' is false), so the same one-cluster
* input stopped the reghdfe path with rc 459 and completed under poisson with every SE
* missing. the guard now also reads e(N_clust), which ppmlhdfe does post. P5a is the leg
* that flips (rc 0 on the unpatched build); P5b is the reghdfe parity anchor; P5c is an
* unclustered poisson run, where e(N_clust) is unposted and the added clause must stay false.
use `b', clear
gen byte one = 1
cap qui wooldid ycount i t ttre, att cluster(one) poisson esprelength(1) espostlength(1)
local p5arc = _rc
di "P5a_pois_onecluster_rc=" `p5arc'
if `p5arc' != 459 {
  di as error "P5a FAIL: poisson + cluster(one) rc " `p5arc' " (expected 459, as on the reghdfe path)"
  local ++fails
}
use `b', clear
gen byte one = 1
cap qui wooldid y i t ttre, att cluster(one) esprelength(1) espostlength(1)
local p5brc = _rc
di "P5b_reghdfe_onecluster_rc=" `p5brc'
if `p5brc' != 459 {
  di as error "P5b FAIL: reghdfe + cluster(one) rc " `p5brc' " (expected 459)"
  local ++fails
}
use `b', clear
cap qui wooldid ycount i t ttre, att robust poisson esprelength(1) espostlength(1)
local p5crc = _rc
local p5cse .
cap matrix P5 = e(wdidmainresults)
if _rc == 0 local p5cse = P5[1,2]
di "P5c_pois_robust_rc=" `p5crc' " se=" `p5cse'
if `p5crc' != 0 | `p5cse' >= . | `p5cse' <= 0 {
  di as error "P5c FAIL: unclustered poisson rc " `p5crc' " or non-finite SE (the cluster-count clause must not fire when e(N_clust) is unposted)"
  local ++fails
}

* ---- P6: poisson-dyex-title-semielasticity (change audit round 2, U13/S3; author ruling
* 2026-09-23, text approved). The contreat labels keyed on contreatelasticitytype alone, so under
* poisson a dyex figure was titled "Semi-elasticities (dyex)" and a dydx one "Effects", the
* reverse of the help: on the linear index dydx is the semi-elasticity of E[y] and dyex the
* elasticity. Under poisson with contreat() dyex is now "Elasticities (dyex)" and dydx
* "Semi-elasticities (dydx)" in the cfx titles (e(wdid_cfxttl1)), the ES plot's AME title, the
* ame histogram's axis title and the table header; linear runs and binary poisson runs keep
* their labels. The plots are read from the saved .gph files, the header from a log (none of
* the needles is in the echoed command). Each check is taken as the file is read, so no file's
* contents ever pass through a macro
local p6ls = c(linesize)
set linesize 255
local p6n 0
* case: poisson flag, operator, then the needles that must be present in the cfx title, the
* header log, the ES ame plot and the ame_i histogram ("-" = no check), then a needle that must
* be ABSENT from the header log
local p6c1 pdyex poisson dyex "Post-treatment Elasticities (dyex) Across" "Presented as Elasticities (dyex)" "AME Elasticities (dyex)" "Treatment Elasticities (dyex)" "Semi-elasticities"
local p6c2 pdydx poisson dydx "Post-treatment Semi-elasticities (dydx) Across" "Presented as Semi-elasticities (dydx)" "AME Semi-elasticities (dydx)" "Treatment Semi-elasticities (dydx)" "Presented as Elasticities"
local p6c3 ldyex "" dyex "Post-treatment Semi-elasticities (dyex) Across" "Presented as Semi-elasticities (dyex)" "AME Semi-elasticities (dyex)" "Treatment Semi-elasticities (dyex)" "Presented as Elasticities"
local p6c4 ldydx "" dydx "Post-treatment Effects Across" "-" "AME Effects" "Treatment Effects" "Presented as"
local p6c5 pbin poisson "" "-" "-" "-" "-" "Presented as"
forval k = 1/5 {
  local c : word 1 of `p6c`k''
  local p6opt : word 2 of `p6c`k''
  local p6cet : word 3 of `p6c`k''
  local nt : word 4 of `p6c`k''
  local nh : word 5 of `p6c`k''
  local ne : word 6 of `p6c`k''
  local ni : word 7 of `p6c`k''
  local nx : word 8 of `p6c`k''
  use "`REPO'/test/data/dta_cont.dta", clear
  qui xtile p6bin = intensity, nq(3)
  local p6ct contreat(intensity) contreatcontrols(basic) contreatelasticitytype(`p6cet') makecfxplots cfxlattice(p6bin)
  if "`p6cet'" == "" local p6ct
  cap erase p6_`c'.txt
  qui log using p6_`c'.txt, text replace name(p6log)
  cap noi wooldid y i t ttre, att att_i `p6ct' cluster(i) `p6opt' safety(off) esprelength(1) espostlength(1) makeplots histogramcohorteffects saveplots(p6`c') replace
  local p6rc = _rc
  qui log close p6log
  local p6t `"`e(wdid_cfxttl1)'"'
  local miss
  if `p6rc' != 0 local miss `miss' rc`p6rc'
  if "`nt'" != "-" & strpos(`"`p6t'"', "`nt'") == 0 local miss `miss' cfxtitle
  if "`nh'" != "-" & strpos(fileread("p6_`c'.txt"), "`nh'") == 0 local miss `miss' header
  if strpos(fileread("p6_`c'.txt"), "`nx'") > 0 local miss `miss' header-absent
  if "`ne'" != "-" {
    cap confirm file p6`c'_ES_ame.gph
    if _rc != 0 local miss `miss' esfile
    else if strpos(fileread("p6`c'_ES_ame.gph"), "`ne'") == 0 local miss `miss' esplot
  }
  if "`ni'" != "-" {
    cap confirm file p6`c'_hist_ame_i.gph
    if _rc != 0 local miss `miss' histfile
    else if strpos(fileread("p6`c'_hist_ame_i.gph"), "`ni'") == 0 local miss `miss' hist
  }
  di `"P6_`c' rc=`p6rc' cfxttl1=[`p6t'] misses=[`miss']"'
  if "`miss'" != "" local ++p6n
  local p6f : dir . files "p6`c'*"
  foreach f of local p6f {
    cap erase "`f'"
  }
  cap erase p6_`c'.txt
}
set linesize `p6ls'
if `p6n' > 0 {
  di as error "P6 FAIL: `p6n' of 5 poisson/linear dose-effect label cases missed a title check"
  local ++fails
}

if `fails' == 0 di as result "POISMARGINS_ALL_PASS"
else di as error "POISMARGINS_FAILS=`fails'"
