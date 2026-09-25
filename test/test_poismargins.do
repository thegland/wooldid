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
* P7 poisson + unconditionalse + cluster() posts e(ll_0) and e(r2pseud) and prints the
*    pseudo-R2: the stage-1 command spells vce(cluster ...) (author ruling 2026-09-23)
* P8 poisson runs whose absorbed terms carry two or more slopes complete instead of stopping
*    r(3201): the dose polynomial enters explicitly, ccc_absorb() is split one slope per term
*    (author ruling 2026-09-23)
* P10 plain basic poisson fits first time with the explicit pooled dose control: no r(3200)
*    caution, and the results pinned from the pre-7647081 re-fit; a one-level fe() slope of the
*    user's still cautions and stops r(3200) (author ruling 2026-09-23)
* P11 a dose with a large level next to its spread within the fixed effects makes ppmlhdfe drop
*    the pooled dose control; the run now cautions, and stays silent on an ordinary run and where
*    the fixed effects really do carry the pooled slope (author ruling 2026-09-24)
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

* ---- P7: poisson-uncondse-pseudor2-dot (change audit round 2, S3; author ruling 2026-09-23).
* On the unconditionalse path the stage-1 poisson got the legacy cluster() spelling, under which
* Stata's poisson posts e(ll_0) and e(r2_p) as missing, so the header printed "pseudo-R2 = ."
* while robust and customvce(cluster ...) runs printed a value. The command now spells
* vce(cluster ...). e(b) and e(V) do not move (measured bitwise when the change was made, on
* regress and poisson, weighted and unweighted); what moves is that e(ll_0) and e(r2pseud) are
* posted and the header prints the number. e(wdid_vcespec), the header's "Standard Errors"
* label, keeps the user's spelling
use `b', clear
cap erase p7.txt
qui log using p7.txt, text replace name(p7log)
cap noi wooldid ycount i t ttre, att cluster(i) poisson unconditionalse esprelength(1) espostlength(1)
local p7rc = _rc
qui log close p7log
scalar p7r2 = e(r2pseud)
scalar p7ll0 = e(ll_0)
local p7ok 0
if `p7rc' == 0 & p7r2 < . & p7ll0 < . {
  local p7lbl = strpos(fileread("p7.txt"), "pseudo-R2 = " + string(p7r2, "%-6.4f")) > 0
  if reldif(p7r2, 1 - e(ll)/p7ll0) < 1e-12 & `p7lbl' & "`e(wdid_vcespec)'" == "cluster(i)" local p7ok 1
}
cap erase p7.txt
di "P7_uncond_cluster_rc=" `p7rc' " r2pseud=" p7r2 " ll_0=" p7ll0 " vcespec=`e(wdid_vcespec)'"
if !`p7ok' {
  di as error "P7 FAIL: poisson + unconditionalse + cluster() must post e(ll_0) and e(r2pseud) (= 1 - ll/ll_0) and print the pseudo-R2, with e(wdid_vcespec) still cluster(i)"
  local ++fails
}

* ---- P8: ppml-multislope-3201 (author ruling 2026-09-23). ppmlhdfe 2.3.0's slope
* preconditioner stops r(3201) on any absorbed factor carrying two or more slopes, and wooldid
* built two such terms itself: the dose polynomial under contreatpoly(2+) with a lattice or a
* byX control, and ccc_absorb() of two or more variables. Both runs stopped r(3201) and posted
* nothing. The dose polynomial now enters the regression explicitly with the ccnoabsorb spelling
* (splitting it across copies of its factor converged far too slowly), so P8a must equal its
* ccnoabsorb twin bit for bit. ccc_absorb() is split into one absorbed term per variable, so
* P8b's point estimates must equal the explicit coarsecohortcontrols(c.x1 c.x2) twin to 1e-8,
* and its standard error, which the explicit twin loses here (all-zero VCE), must match Stata's
* own poisson on the explicit terms (customvce(robust)) to that command's convergence floor
local CCD "`REPO'/test/suite/data/continuous_control_dense.dta"
local p8n 0
use "`CCD'", clear
cap qui wooldid y_count cohort t ttre, contreat(dose) contreatcontrols(byTreated) contreatpoly(2) poisson robust
local p8arc = _rc
tempname p8ab p8bb p8bv p8cb p8ob p8ov
if `p8arc' == 0 matrix `p8ab' = e(b)
use "`CCD'", clear
cap qui wooldid y_count cohort t ttre, contreat(dose) contreatcontrols(byTreated) contreatpoly(2) poisson robust ccnoabsorb
local p8anrc = _rc
local p8ad .
if `p8arc' == 0 & `p8anrc' == 0 local p8ad = mreldif(`p8ab', e(b))
di "P8a_poly_bytreated rc=`p8arc' twin_rc=`p8anrc' mreldif_b_vs_ccnoabsorb=`p8ad'"
if `p8arc' != 0 | `p8ad' != 0 {
  di as error "P8a FAIL: poisson + byTreated + contreatpoly(2) must run and equal its ccnoabsorb twin bit for bit"
  local ++p8n
}
use "`CCD'", clear
cap qui wooldid y_count cohort t ttre, poisson robust ccc_absorb(x1 x2)
local p8brc = _rc
local p8bse .
if `p8brc' == 0 {
  matrix `p8bb' = e(b)
  cap matrix `p8bv' = e(V)
  if _rc == 0 local p8bse = sqrt(`p8bv'[1,1])
}
use "`CCD'", clear
cap qui wooldid y_count cohort t ttre, poisson robust coarsecohortcontrols(c.x1 c.x2)
local p8crc = _rc
local p8bd .
if `p8brc' == 0 & `p8crc' == 0 local p8bd = mreldif(`p8bb', e(b))
use "`CCD'", clear
cap qui wooldid y_count cohort t ttre, poisson customvce(robust) ccc_absorb(x1 x2)
local p8orc = _rc
local p8sd .
if `p8orc' == 0 & `p8bse' < . {
  matrix `p8ov' = e(V)
  local p8sd = reldif(`p8bse', sqrt(`p8ov'[1,1]))
}
di "P8b_ccc_two rc=`p8brc' se=`p8bse' mreldif_b_vs_explicit=`p8bd' reldif_se_vs_poisson=`p8sd'"
if `p8brc' != 0 | `p8bd' > 1e-8 | `p8sd' > 1e-4 {
  di as error "P8b FAIL: poisson + ccc_absorb(x1 x2) must run, match the explicit spelling to 1e-8 and Stata's poisson SE to 1e-4"
  local ++p8n
}
if `p8n' > 0 local ++fails

* ---- P10: ppml-plainbasic-first-call (commit 7647081, author ruling 2026-09-23). On plain basic
* poisson the pooled dose control i.null##c.(ringer) is a slope on a one-level factor, which
* ppmlhdfe 2.3.0's slope preconditioner cannot absorb: every such run failed its first fit
* r(3200), printed the four-line fallback caution and reported the re-fit with the control as
* an explicit regressor. That spelling is now built up front. P10a: no caution prints and the
* results are the old re-fit's, e(b) and the diagonal of e(V) pinned at %21x from a run of
* 7647081^ (bfbaa81) in its own Stata session. P10b: a one-level slope of the user's own in fe()
* still prints the caution and stops r(3200), as before. The needle is read from logs closed
* before it is used, and appears in no command echoed while they are open
local p10ls = c(linesize)
set linesize 255
local p10needle "initial estimation of ppmlhdfe failed"
local p10pin main_att +1.6029ba44662f0X-003 +1.c1b42ebbbbe6aX-013 pretreatment_att -1.85401f6bde400X-009 +1.349c47b0b6b2aX-013 ///
  main_ame +1.3df4d64decd86X-005 +1.65ba10d9a9715X-00e pretreatment_ame -1.48c44fcd7c90eX-005 +1.06bcda5ae0e79X-00e ///
  ES_att_pre1 -1.85401f6bde400X-009 +1.349c47b0b6b2aX-013 ES_att_post0 +1.c1c9585e38ae0X-004 +1.1412fc44a2574X-012 ///
  ES_ame_pre1 -1.48c44fcd7c90eX-005 +1.06bcda5ae0e79X-00e ES_ame_post0 +1.095c60902dbb9X-004 +1.44e79eef59917X-00d
use "`CCD'", clear
cap erase p10a.txt
qui log using p10a.txt, text replace name(p10log)
cap noi wooldid y_count cohort t ttre, contreat(dose) contreatcontrols(basic) poisson robust esprelength(1) espostlength(1)
local p10arc = _rc
qui log close p10log
local p10acau = strpos(fileread("p10a.txt"), "`p10needle'") > 0
cap erase p10a.txt
local p10miss
tempname p10b p10v
if `p10arc' == 0 {
  matrix `p10b' = e(b)
  cap matrix `p10v' = e(V)
  if _rc != 0 local p10miss noV
  else {
    tokenize `p10pin'
    local j 1
    while "``j''" != "" {
      local jb = `j' + 1
      local jv = `j' + 2
      local c = colnumb(`p10b', "``j''")
      local hb .
      local hv .
      if `c' < . {
        local hb : di %21x `p10b'[1, `c']
        local hv : di %21x `p10v'[`c', `c']
      }
      if "`hb'" != "``jb''" | "`hv'" != "``jv''" local p10miss `p10miss' ``j''
      local j = `j' + 3
    }
  }
}
use "`CCD'", clear
gen byte one = 1
cap erase p10u.txt
qui log using p10u.txt, text replace name(p10log)
cap noi wooldid y_count cohort t ttre, contreat(dose) contreatcontrols(basic) poisson robust fe(i.one#c.x1)
local p10urc = _rc
qui log close p10log
local p10ucau = strpos(fileread("p10u.txt"), "`p10needle'") > 0
* P10c (author ruling 2026-09-25, card m034): the fallback's second line names only the
* ccc_absorb() and continuous treatment controls; time-trend controls never enter absorb()
local p10unew = strpos(fileread("p10u.txt"), "Falling back to estimation with ppmlhdfe, with any ccc_absorb() controls and continuous treatment controls moved out of the absorb() statement and into the main regression.") > 0
local p10uold = strpos(fileread("p10u.txt"), "time-trend controls") > 0
cap erase p10u.txt
set linesize `p10ls'
di "P10a_plainbasic rc=`p10arc' caution=`p10acau' pin_misses=[`p10miss'] | P10b_user_onelevel rc=`p10urc' caution=`p10ucau'"
di "P10c fallback second line: new=`p10unew' old=`p10uold'"
if !`p10unew' | `p10uold' {
  di as error "P10c FAIL: the r(3200) fallback's second line must carry the 2026-09-25 text"
  local ++fails
}
if `p10arc' != 0 | `p10acau' | "`p10miss'" != "" | `p10urc' != 3200 | !`p10ucau' {
  di as error "P10 FAIL: plain basic poisson must fit first time with the old re-fit's results, and a one-level fe() slope must still caution and stop r(3200)"
  local ++fails
}

* ---- P11: ppml-dropped-pooled-slope (author ruling 2026-09-24, text approved). Wherever the
* pooled dose control i.null##c.(ringer) is an explicit regressor, ppmlhdfe zeroes it when its sum
* of squares after the fixed effects is at most 1e-6 of its raw, uncentred sum of squares, and the
* c.dose twin goes with it, so a dose whose level is large next to its spread within the fixed
* effects silently loses the pooled slope. On P11a's panel (the 2026-09-24 probe's sweep panel,
* seed 505, dose + 300) ame comes out .5196 against .2373 unshifted, 6.5 standard errors off, with
* no caution under robust. The run now cautions. P11a: that panel prints the caution. P11b: P10's
* ordinary plain basic design does not. P11d: the shifted panel under contreatpoly(2), where only the
* linear pooled term is dropped, prints it too. P11c: the same panel unshifted under fe(pid), where the
* dose is constant within pid, so ppmlhdfe drops the control for the right reason (the fixed
* effects carry the pooled slope and ame is correct): no caution. P11a is the leg that fails on
* the unpatched build. The needle is read from logs closed before it is used, and appears in no
* command echoed while they are open
local p11ls = c(linesize)
set linesize 255
local p11needle "dropped the pooled"
clear
set seed 505
set obs 200
gen int i = ceil(_n/20)
gen int pid = _n
gen double base = 0.5 + runiform()
bys i: egen double cb = mean(base)
gen double dose = cb + (base - cb)
drop base cb
expand 8
bys pid: gen int t = _n
gen ttre = cond(i <= 2, 4, cond(i <= 4, 5, cond(i <= 6, 6, .)))
gen byte D = t >= ttre & ttre < .
gen double mu = exp(1.5 + 0.05*i + 0.04*t + 0.3*dose + D*(0.2 + 0.25*dose))
gen long y = rpoisson(mu)
tempfile p11d
save `p11d'
replace dose = dose + 300
cap erase p11a.txt
qui log using p11a.txt, text replace name(p11log)
cap noi wooldid y i t ttre, contreat(dose) contreatcontrols(basic) poisson robust
local p11arc = _rc
qui log close p11log
local p11acau = strpos(fileread("p11a.txt"), "`p11needle'") > 0
cap erase p11a.txt
use "`CCD'", clear
cap erase p11b.txt
qui log using p11b.txt, text replace name(p11log)
cap noi wooldid y_count cohort t ttre, contreat(dose) contreatcontrols(basic) poisson robust esprelength(1) espostlength(1)
local p11brc = _rc
qui log close p11log
local p11bcau = strpos(fileread("p11b.txt"), "`p11needle'") > 0
cap erase p11b.txt
use `p11d', clear
cap erase p11c.txt
qui log using p11c.txt, text replace name(p11log)
cap noi wooldid y i t ttre, contreat(dose) contreatcontrols(basic) poisson robust fe(pid)
local p11crc = _rc
qui log close p11log
local p11ccau = strpos(fileread("p11c.txt"), "`p11needle'") > 0
cap erase p11c.txt
* P11d: under contreatpoly(2) ppmlhdfe can drop the linear pooled term and keep the square;
* the pooled slope is gone all the same, so the caution must print
use `p11d', clear
replace dose = dose + 300
cap drop _ppmlhdfe_d
cap erase p11dd.txt
qui log using p11dd.txt, text replace name(p11log)
cap noi wooldid y i t ttre, contreat(dose) contreatcontrols(basic) contreatpoly(2) poisson robust
local p11drc = _rc
qui log close p11log
local p11dcau = strpos(fileread("p11dd.txt"), "`p11needle'") > 0
cap erase p11dd.txt
cap drop _ppmlhdfe_d
set linesize `p11ls'
di "P11a_shifted rc=`p11arc' caution=`p11acau' | P11b_ordinary rc=`p11brc' caution=`p11bcau' | P11c_fe_unit rc=`p11crc' caution=`p11ccau' | P11d_poly2_shifted rc=`p11drc' caution=`p11dcau'"
if `p11arc' != 0 | !`p11acau' | `p11brc' != 0 | `p11bcau' | `p11crc' != 0 | `p11ccau' | `p11drc' != 0 | !`p11dcau' {
  di as error "P11 FAIL: the shifted-dose runs (linear and contreatpoly(2)) must print the dropped-slope caution, and the ordinary and fe(pid) runs must not"
  local ++fails
}

if `fails' == 0 di as result "POISMARGINS_ALL_PASS"
else di as error "POISMARGINS_FAILS=`fails'"
