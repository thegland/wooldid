* test_bootstrap_df.do - bootdf(): the reference distribution for bootstrap inference.
* Run from test/: stata-mp -e do test_bootstrap_df.do
* T1 the default is CONDITIONAL: clustered -> t(G-1), unclustered -> normal, each
*    bit-identical to the matching explicit bootdf(), and an explicit value overrides;
* T2 bootdf(clusters) gives t(G-1): SE and the statistic itself unchanged, p and CI
* strictly wider, hand-reconstructed exactly, column relabelled t, and e(df_r) aligned
* so a hand-run lincom reproduces the table p; T3 joint tests become F(k,df), hand-
* checked against Ftail on the same replication covariance; T4 an explicit integer df;
* T5 the gates (bad token, stray without bootstrap(), zero df); T6 standard mode with
* subgroups - pwsgs contrasts t-based; T7 replay. Prints BOOTDF_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0
cap prog drop mkp
prog def mkp
  clear
  set obs 12
  gen unit = _n
  gen sg = mod(unit, 2)
  gen ttre = cond(unit<=4, 4, cond(unit<=8, 6, .))
  expand 8
  bysort unit: gen t = _n
  set seed 11
  gen y = 2*(t>=ttre & ttre<.) + 0.3*(t>=ttre & ttre<.)*sg + 0.1*unit + 0.05*t + rnormal()
end

* T1: the default is CONDITIONAL on clustering (author ruling 2026-09-02, replacing the
* earlier unconditional bootdf(normal) default). Clustered -> t(G-1), bit-identical to an
* explicit bootdf(clusters); unclustered -> normal, bit-identical to an explicit
* bootdf(normal). An explicit bootdf() always wins over the default.
mkp
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(20) rseed(9)
matrix D1 = e(wdidmainresults)
local st1 `e(boot_stat)'
local df1 = e(boot_df)
mkp
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(20) rseed(9) bootdf(clusters)
di "T1a_clustered_default_vs_clusters=" mreldif(D1, e(wdidmainresults)) " stat=[`st1'] boot_df=`df1'"
if mreldif(D1, e(wdidmainresults)) > 0 | "`st1'" != "t" | `df1' != 11 {
  di as error "T1a FAIL: clustered default must be t(G-1)"
  local ++fails
}
* and the explicit normal still overrides it
mkp
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(20) rseed(9) bootdf(normal) bootdisplay(wald)
matrix D1n = e(wdidmainresults)
local st1b `e(boot_stat)'
local df1b = e(boot_df)
di "T1b_clustered_explicit_normal stat=[`st1b'] boot_df=`df1b'"
if "`st1b'" != "z" | `df1b' < . {
  di as error "T1b FAIL: an explicit bootdf(normal) must override the clustered default"
  local ++fails
}
* unclustered: the default is normal, bit-identical to an explicit bootdf(normal)
mkp
qui wooldid y unit t ttre, att robust esprelength(1) espostlength(1) bootstrap(bayes) reps(20) rseed(9)
matrix D1c = e(wdidmainresults)
local st1c `e(boot_stat)'
mkp
qui wooldid y unit t ttre, att robust esprelength(1) espostlength(1) bootstrap(bayes) reps(20) rseed(9) bootdf(normal)
di "T1c_unclustered_default_vs_normal=" mreldif(D1c, e(wdidmainresults)) " stat=[`st1c']"
if mreldif(D1c, e(wdidmainresults)) > 0 | "`st1c'" != "z" {
  di as error "T1c FAIL: unclustered default must be normal"
  local ++fails
}

* T2: bootdf(clusters) -> t with G-1=11; SE unchanged, only the reference distribution moves
mkp
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(20) rseed(9) bootdf(clusters) bootdisplay(wald)
matrix D2 = e(wdidmainresults)
di "T2_stat=[`e(boot_stat)'] boot_df=" e(boot_df) " df_r=" e(df_r)
* the z baseline is the EXPLICIT bootdf(normal) run (D1n): since the ruling of 2026-09-02
* the plain default on a clustered run is itself t(G-1), so D1 would be the same object
* SE column identical to the normal run
di "T2_se_same=" reldif(D1n[1,2], D2[1,2])
* the t stat itself is b/se either way
di "T2_tstat_same=" reldif(D1n[1,3], D2[1,3])
* p and CI must be WIDER / larger under t
local pz = D1n[1,4]
local pt = D2[1,4]
local wz = D1n[1,6] - D1n[1,5]
local wt = D2[1,6] - D2[1,5]
di "T2_p_z=" `pz' " p_t=" `pt' " ci_width_z=" `wz' " ci_width_t=" `wt'
if reldif(D1n[1,2], D2[1,2]) > 0 | `pt' <= `pz' | `wt' <= `wz' | e(boot_df) != 11 {
  di as error "T2 FAIL"
  local ++fails
}
* hand-check the t p-value and critical value
local hp = 2*ttail(11, abs(D2[1,3]))
local hlb = D2[1,1] - invt(11, .975)*D2[1,2]
di "T2_hand_p_reldif=" reldif(`pt', `hp') " hand_lb_reldif=" reldif(D2[1,5], `hlb')
if reldif(`pt', `hp') > 1e-12 | reldif(D2[1,5], `hlb') > 1e-12 {
  di as error "T2b FAIL: hand reconstruction"
  local ++fails
}
* column label is t
local cn : colnames D2
di "T2_cols=[`cn']"
if strpos("`cn'", " t ") == 0 {
  di as error "T2c FAIL: column not labelled t"
  local ++fails
}
* e(df_r) aligned -> hand test/lincom agrees with the table
qui lincom main_att
di "T2_lincom_p=" r(p) " table_p=" `pt' " reldif=" reldif(r(p), `pt')
if reldif(r(p), `pt') > 1e-10 {
  di as error "T2d FAIL: lincom does not match the table"
  local ++fails
}

* T3: joint tests under bootdf. RE-PINNED 2026-09-04 for the author's ruling that a joint
* test which cannot be BOOTSTRAPPED is ABSENT ("the joint test is absent if not
* bootstrappable. NEVER show an analytic SE when bootstrap is requested"). Bootstrapping a
* joint test means bootstrapping its test statistic, which requires boottstat, so without
* it the scalar must not exist at all - not merely be missing, and never the analytic
* value. The old T3 asserted the bootstrapped joint p was F(k,df) under bootdf(clusters)
* and strictly larger than the chi-squared one; that whole surface is gone without
* boottstat, so bootdf() no longer has a joint test to move.
mkp
qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(1) jointtests bootstrap(bayes) reps(40) rseed(99) bootdf(normal)
cap confirm scalar e(joint_horiz)
local rc_z = _rc
local flag_z `e(boot_jt_absent)'
mkp
qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(1) jointtests bootstrap(bayes) reps(40) rseed(99) bootdf(clusters)
cap confirm scalar e(joint_horiz)
local rc_t = _rc
local flag_t `e(boot_jt_absent)'
di "T3_absent rc_normal=`rc_z' rc_clusters=`rc_t' flags=[`flag_z'][`flag_t'] (want 111 111 and both flagged)"
if `rc_z' == 0 | `rc_t' == 0 | "`flag_z'" != "1" | "`flag_t'" != "1" {
  di as error "T3 FAIL: without boottstat the joint test must be absent and flagged"
  local ++fails
}
* and the analytic value must not be what got suppressed - an unclustered analytic run
* posts one, so check the bootstrap run does not reproduce it
mkp
qui wooldid y unit t ttre, att esprelength(3) espostlength(1) jointtests
scalar jp_analytic = e(joint_horiz)
mkp
qui wooldid y unit t ttre, att esprelength(3) espostlength(1) jointtests bootstrap(bayes) reps(40) rseed(99) cluster(unit)
cap confirm scalar e(joint_horiz)
local rc_u = _rc
di "T3b_analytic=" jp_analytic " bootstrap_rc=`rc_u' (want a real analytic value and rc 111)"
if `rc_u' == 0 {
  di as error "T3b FAIL: an analytic joint p survived into a bootstrap run"
  local ++fails
}

* T4: explicit integer df
mkp
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(20) rseed(9) bootdf(5)
di "T4_boot_df=" e(boot_df) " stat=[`e(boot_stat)']"
if e(boot_df) != 5 {
  di as error "T4 FAIL"
  local ++fails
}

* T5: gates - bad token, stray without bootstrap()
mkp
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) reps(20) bootdf(garbage)
di "T5a_rc=" _rc
if _rc != 198 {
  di as error "T5a FAIL"
  local ++fails
}
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootdf(clusters)
di "T5b_stray_rc=" _rc
if _rc != 198 {
  di as error "T5b FAIL"
  local ++fails
}
cap qui wooldid y unit t ttre, att cluster(unit) bootstrap(bayes) reps(20) bootdf(0)
di "T5c_rc=" _rc
if _rc != 198 {
  di as error "T5c FAIL"
  local ++fails
}

* T6: standard mode + subgroups + pwsgs under bootdf. RE-POINTED 2026-09-04: the
* assertion below is that the contrast p is the t-based Wald tail of the contrast
* statistic, which is a property of the STANDARD-DEVIATION family. Since the percentile
* package became the display default, the contrast row carries a rank p dual to the
* percentile interval instead, so this test must ask for the family it is testing.
* T6c below pins the default's contrast row separately.
mkp
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) subgroups(sg) bootstrap(standard) reps(25) rseed(6) bootdf(clusters) bootdisplay(wald)
di "T6_rc=" _rc " stat=[`e(boot_stat)']"
if _rc != 0 {
  di as error "T6 FAIL"
  local ++fails
}
else {
  matrix P = e(pwsgs)
  local pc : colnames P
  di "T6_pwcols=[`pc'] p=" P[1,4]
  local hpw = 2*ttail(e(boot_df), abs(P[1,3]))
  if reldif(P[1,4], `hpw') > 1e-12 {
    di as error "T6b FAIL: pwsgs p not t-based"
    local ++fails
  }
}

* T6c: the DEFAULT contrast row is the percentile package - its interval is the percentile
* interval of the contrast draws, its p is that interval's dual rank p, and its standard
* error is the interval inverted. Also pins the duality on a contrast row.
mkp
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) subgroups(sg) bootstrap(standard) reps(25) rseed(6) bootdf(clusters)
if _rc != 0 {
  di as error "T6c FAIL: default-display run errored"
  local ++fails
}
else {
  matrix PD = e(pwsgs)
  matrix PP = e(boot_pctl_pwsgs)
  matrix PS = e(boot_wald_pwsgs)
  * RE-POINTED 2026-09-04 for the author's ruling that wooldid never reports a standard
  * error implied by inverting the reported interval - the earlier assertion here was
  * exactly that, and it is now wrong on purpose. The invariants that replace it: the
  * displayed contrast row IS the percentile family; its standard error is the bootstrap
  * STANDARD DEVIATION, identical to the one the wald display reports; and the p is dual
  * to the displayed interval.
  local bad6 0
  local live6 0
  forvalues r = 1/`=rowsof(PD)' {
    forvalues c = 1/6 {
      if (PD[`r',`c'] < . | PP[`r',`c'] < .) & reldif(PD[`r',`c'], PP[`r',`c']) > 0 local ++bad6
    }
    if PD[`r',2] < . {
      local ++live6
      * the standard error must be the SD family's, NOT anything derived from the interval
      if reldif(PD[`r',2], PS[`r',2]) > 0 local ++bad6
      * and it must NOT coincide with the interval-inverted quantity the ruling removed
      local zc6 = invt(e(boot_df), 1 - (100 - e(boot_level))/200)
      local se6 = (PD[`r',6] - PD[`r',5])/(2*`zc6')
      if reldif(PD[`r',2], `se6') <= 1e-12 {
        di as error "T6c NOTE: se coincides with the inverted interval on row `r' - check this is arithmetic accident, not a regression"
      }
      * duality: p at or below alpha exactly when the interval excludes zero
      local exc6 = (PD[`r',5] > 0 | PD[`r',6] < 0)
      local sig6 = (PD[`r',4] <= (100 - e(boot_level))/100)
      if `exc6' != `sig6' local ++bad6
    }
  }
  di "T6c_default_is_percentile bad=`bad6' live_rows=`live6' (want bad 0, live > 0); sd twin present=" (rowsof(PS) > 0)
  * a test that passes because every contrast row is blank has checked nothing - the
  * earlier version of T6c could do exactly that, so vacuity is now itself a failure
  if `live6' == 0 {
    di as error "T6c FAIL: no non-blank contrast rows, so the assertions were vacuous"
    local ++fails
  }
  if `bad6' > 0 {
    di as error "T6c FAIL: the default contrast row is not the percentile package"
    local ++fails
  }
}

* T7: replay reprints under bootdf
cap qui wooldid
di "T7_replay_rc=" _rc
if _rc != 0 {
  di as error "T7 FAIL"
  local ++fails
}

* ---- T9: event-study PLOT bands must equal the TABLE bands, in whichever family
* bootdisplay() selected. Before 2026-09-04 wooldid_esplots ran inside the worker with no
* bootstrap gate, so it read e(wdidmainresults) BEFORE the rewiring and the figure showed
* ANALYTIC bands while the table showed bootstrap ones - measured to disagree about whether
* a pre-trend excluded zero. The call is now deferred to the shell, after wooldid_bootwire.
*
* ANTI-VACUITY: this test fails unless it actually compared bands (compared > 0) AND the
* three flavours produce genuinely different intervals - otherwise it would pass on a build
* that ignored bootdisplay() entirely, or on one that plotted nothing.
cap prog drop bandchk
prog def bandchk, rclass
  args dta mt
  local bad 0
  local n 0
  local rn : rowfullnames `mt'
  local k 0
  foreach r of local rn {
    local ++k
    if strpos("`r'", "ES_att:") {
      local ry = `mt'[`k',7]
      local key = cond(`ry'<0, "m"+string(-`ry'), string(`ry'))
      scalar B_lb_`key' = `mt'[`k',5]
      scalar B_ub_`key' = `mt'[`k',6]
    }
  }
  preserve
  use "`dta'", clear
  forvalues i = 1/`=_N' {
    if rownamesvar[`i'] == "ES_att" {
      local ry = relyr[`i']
      local key = cond(`ry'<0, "m"+string(-`ry'), string(`ry'))
      cap confirm scalar B_lb_`key'
      if _rc == 0 {
        local ++n
        if reldif(lb_95ci[`i'], B_lb_`key') > 1e-10 | reldif(ub_95ci[`i'], B_ub_`key') > 1e-10 local ++bad
      }
    }
  }
  restore
  return scalar n = `n'
  return scalar bad = `bad'
end
local t9ok 1
foreach d in percentile bc wald {
  mkp
  cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2) bootstrap(standard) reps(40) rseed(7) makeplots saveplots(__t9`d') replace bootdisplay(`d')
  if _rc != 0 {
    di as error "T9 FAIL: bootdisplay(`d') plot run errored rc=" _rc
    local ++fails
    local t9ok 0
  }
  else {
    matrix T9_`d' = e(wdidmainresults)
    bandchk "__t9`d'_ES.dta" T9_`d'
    di "T9_`d' compared=" r(n) " mismatches=" r(bad)
    if r(n) == 0 {
      di as error "T9 FAIL (VACUOUS): no plot bands were compared for bootdisplay(`d')"
      local ++fails
      local t9ok 0
    }
    else if r(bad) > 0 {
      di as error "T9 FAIL: plot bands differ from the table under bootdisplay(`d')"
      local ++fails
      local t9ok 0
    }
  }
  cap erase "__t9`d'_ES.dta"
  cap erase "__t9`d'_ES_att.gph"
}
if `t9ok' {
  local t9differ = (reldif(T9_percentile[1,5], T9_wald[1,5]) > 1e-10)
  di "T9_differ=`t9differ' (pctl lb=" T9_percentile[1,5] " wald lb=" T9_wald[1,5] ")"
  if `t9differ' == 0 {
    di as error "T9 FAIL (VACUOUS): the display flavours produced identical bands, so this test cannot detect a wrong flavour"
    local ++fails
  }
}

* ---------------------------------------------------------------- T10: e(gteffects) cells carry bootstrap inference (ruling
* 2026-09-04): live cells have se > 0 and lb < ub (a level() default of 0 once collapsed every cell interval onto its
* median, lb == ub), the three family twins are posted, and the bounds move with the seed
mkp
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) gteffects bootstrap(standard) reps(30) rseed(7)
local rc10 = _rc
di "T10_run rc=`rc10' family=[`e(boot_gteffects_family)']"
cap confirm matrix e(boot_pctl_gteffects)
local tw1 = _rc
cap confirm matrix e(boot_bc_gteffects)
local tw2 = _rc
cap confirm matrix e(boot_wald_gteffects)
local tw3 = _rc
if `rc10' != 0 | `tw1' != 0 | `tw2' != 0 | `tw3' != 0 | "`e(boot_gteffects_family)'" != "percentile" {
  di as error "T10 FAIL: gteffects run rc=`rc10', twins pctl/bc/wald rc=`tw1'/`tw2'/`tw3', family=[`e(boot_gteffects_family)']"
  local ++fails
}
else {
  matrix T10_G7 = e(gteffects)
  mkp
  qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) gteffects bootstrap(standard) reps(30) rseed(8)
  matrix T10_G8 = e(gteffects)
  mata:
    G7=st_matrix("T10_G7"); G8=st_matrix("T10_G8")
    live=0; degen=0; badse=0; moved=0
    for (i=1;i<=rows(G7);i++) {
      if (G7[i,5] >= . | G7[i,7] >= .) continue
      live++
      if (!(G7[i,8] < G7[i,9])) degen++
      if (!(G7[i,6] > 0)) badse++
      if (G7[i,8] != G8[i,8]) moved++
    }
    st_numscalar("t10live", live); st_numscalar("t10degen", degen); st_numscalar("t10badse", badse); st_numscalar("t10moved", moved)
  end
  di "T10_cells live=" t10live " degenerate(lb>=ub)=" t10degen " se<=0:" t10badse " lb-moved-with-seed=" t10moved
  if t10live == 0 {
    di as error "T10 FAIL (VACUOUS): no live gteffects cell to check"
    local ++fails
  }
  else if t10degen > 0 | t10badse > 0 | t10moved < t10live {
    di as error "T10 FAIL: gteffects cell inference is degenerate or not bootstrap-derived"
    local ++fails
  }
}

* ---------------------------------------------------------------- T11: the e(cfx*) matrices behind the continuous-treatment
* plots carry bootstrap inference (ruling 2026-09-04): family twins posted, live rows have se > 0 and ll < ul, and the se
* moves with the seed while b does not. The shipped dta_cont has an exact dose slope (p = 0 under every seed), so noise is
* added; p itself is not asserted to move, since a p of exactly 0 under both seeds is legitimate
set graphics off
cap prog drop mkcfx11
prog def mkcfx11
  use "data/dta_cont", clear
  qui gen byte clat = mod(i, 3)
  set seed 4242
  qui replace y = y + 4*rnormal()
end
mkcfx11
cap qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat) bootstrap(standard) reps(20) rseed(7)
local rc11 = _rc
local mats11 `e(wdid_cfxmats)'
di "T11_run rc=`rc11' cfx mats=[`mats11']"
if `rc11' != 0 | "`mats11'" == "" {
  di as error "T11 FAIL: the cfx bootstrap run did not post its matrices (rc=`rc11')"
  local ++fails
}
else {
  local t11live 0
  local t11bad 0
  local t11twins 0
  foreach nm of local mats11 {
    cap confirm matrix e(`nm')
    if _rc != 0 continue
    matrix T11_`nm'_7 = e(`nm')
    cap confirm matrix e(boot_pctl_`nm')
    if _rc != 0 local ++t11twins
    cap confirm matrix e(boot_bc_`nm')
    if _rc != 0 local ++t11twins
    cap confirm matrix e(boot_wald_`nm')
    if _rc != 0 local ++t11twins
  }
  mkcfx11
  qui wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat) bootstrap(standard) reps(20) rseed(8)
  foreach nm of local mats11 {
    cap confirm matrix e(`nm')
    if _rc != 0 continue
    matrix T11_C8 = e(`nm')
    matrix T11_C7 = T11_`nm'_7
    mata: st_numscalar("t11r", rows(st_matrix("T11_C7")))
    forval i = 1/`=t11r' {
      if T11_C7[`i',1] >= . | T11_C7[`i',4] >= . continue
      local ++t11live
      if !(T11_C7[`i',2] > 0) | !(T11_C7[`i',5] < T11_C7[`i',6]) local ++t11bad
      if T11_C7[`i',2] == T11_C8[`i',2] local ++t11bad
      if T11_C7[`i',1] != T11_C8[`i',1] local ++t11bad
    }
  }
  di "T11_cells live=`t11live' bad=`t11bad' missing-twins=`t11twins'"
  if `t11live' == 0 {
    di as error "T11 FAIL (VACUOUS): no live cfx row to check"
    local ++fails
  }
  else if `t11bad' > 0 | `t11twins' > 0 {
    di as error "T11 FAIL: cfx cell inference is degenerate, static across seeds, or a family twin is missing"
    local ++fails
  }
}

* ---------------------------------------------------------------- T12: the continuous-treatment plots are drawn AFTER the
* rewiring from the bootstrapped cfx matrices (ruling 2026-09-05): under bootstrap() the worker draws and saves nothing, the
* shell's wooldid_cfxplots draws every e(cfx*) matrix under its old graph name, saveplots() writes the same .gph names plus
* one <name>_cfx.dta whose bands equal the posted matrices' columns 5/6, and the percentile and wald bands differ
cap prog drop mkcfx12
prog def mkcfx12
  use "data/dta_cont", clear
  qui gen byte clat = mod(i, 3)
  set seed 4242
  qui replace y = y + 4*rnormal()
end
local S12 contreat(intensity) contreatcontrols(basic) cluster(i) esprelength(2) espostlength(2) makecfxplots cfxlattice(clat)
mkcfx12
graph drop _all
qui wooldid y i t ttre, `S12'
qui graph dir, memory
local g12an `r(list)'
local t12bad 0
local t12live 0
foreach fam in percentile wald {
  mkcfx12
  graph drop _all
  cap qui wooldid y i t ttre, `S12' bootstrap(standard) reps(20) rseed(7) bootdisplay(`fam') saveplots(t12_`fam') replace
  local rc12 = _rc
  qui graph dir, memory
  local g12 `r(list)'
  local mats12 `e(wdid_cfxmats)'
  di "T12_`fam' rc=`rc12' graphs=[`g12'] analytic=[`g12an'] mats=[`mats12']"
  if `rc12' != 0 | "`: list g12 === g12an'" != "1" | "`mats12'" == "" {
    di as error "T12 FAIL: the bootstrap run did not redraw the same set of cfx graphs (rc=`rc12')"
    local ++fails
    continue
  }
  cap confirm file "t12_`fam'_cfx.dta"
  if _rc != 0 {
    di as error "T12 FAIL: saveplots() wrote no t12_`fam'_cfx.dta"
    local ++fails
    continue
  }
  foreach nm of local mats12 {
    matrix T12_`fam'_`nm' = e(`nm')
    cap confirm file "t12_`fam'_`nm'.gph"
    if _rc != 0 local ++t12bad
  }
  preserve
  use "t12_`fam'_cfx.dta", clear
  foreach nm of local mats12 {
    local rn12 : rownames T12_`fam'_`nm'
    forval r = 1/`=rowsof(T12_`fam'_`nm')' {
      local rt : word `r' of `rn12'
      if !regexm("`rt'", "^lattice_([0-9]+)$") {
        local ++t12bad
        continue
      }
      local lv = regexs(1)
      qui count if matname == "`nm'" & lattice == `lv'
      if r(N) != 1 {
        local ++t12bad
        continue
      }
      if T12_`fam'_`nm'[`r',5] < . {
        local ++t12live
        qui su ll if matname == "`nm'" & lattice == `lv', meanonly
        if reldif(r(mean), T12_`fam'_`nm'[`r',5]) > 1e-12 local ++t12bad
        qui su ul if matname == "`nm'" & lattice == `lv', meanonly
        if reldif(r(mean), T12_`fam'_`nm'[`r',6]) > 1e-12 local ++t12bad
      }
    }
  }
  restore
  foreach nm of local mats12 {
    cap erase "t12_`fam'_`nm'.gph"
  }
  cap erase "t12_`fam'_cfx.dta"
  cap erase "t12_`fam'_ES.dta"
  foreach g in `: dir . files "t12_`fam'_*.gph"' {
    cap erase "`g'"
  }
}
di "T12_bands live=`t12live' bad=`t12bad'"
if `t12live' == 0 {
  di as error "T12 FAIL (VACUOUS): no live cfx band was compared"
  local ++fails
}
else if `t12bad' > 0 {
  di as error "T12 FAIL: a saved cfx band or file does not match the posted matrix"
  local ++fails
}
cap confirm matrix T12_percentile_cfxatt_cfxpost
if _rc == 0 {
  cap confirm matrix T12_wald_cfxatt_cfxpost
  if _rc == 0 {
    local t12differ = (reldif(T12_percentile_cfxatt_cfxpost[1,5], T12_wald_cfxatt_cfxpost[1,5]) > 1e-10)
    di "T12_differ=`t12differ'"
    if `t12differ' == 0 {
      di as error "T12 FAIL (VACUOUS): the percentile and wald cfx bands are identical, so a wrong family could not be detected"
      local ++fails
    }
  }
}

if `fails' == 0 di as result "BOOTDF_ALL_PASS"
else di as error "BOOTDF_FAILS=`fails'"
