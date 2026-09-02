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
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(20) rseed(9) bootdf(normal)
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
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(20) rseed(9) bootdf(clusters)
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

* T3: joint tests become F(k, df) under bootdf
mkp
qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(1) jointtests bootstrap(bayes) reps(40) rseed(99) bootdf(normal)
scalar jpz = e(joint_horiz)
mkp
qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(1) jointtests bootstrap(bayes) reps(40) rseed(99) bootdf(clusters)
scalar jpt = e(joint_horiz)
matrix BD = e(bootdist)
matrix OB = e(b)
local cn : colfullnames e(b)
local pcols
local k 0
foreach nm of local cn {
  local ++k
  if strpos("`nm'", "ES_att_pre") local pcols `pcols' `k'
}
mata: A = st_matrix("BD")[., strtoreal(tokens(st_local("pcols")))']
mata: b = st_matrix("OB")[1, strtoreal(tokens(st_local("pcols")))']'
mata: CC = select(A, rowmissing(A) :== 0)
mata: dev = CC :- mean(CC)
mata: S = (dev'*dev)/(rows(CC)-1)
mata: st_numscalar("W", b' * invsym(S) * b)
mata: st_numscalar("kk", cols(A))
local hF = Ftail(kk, 11, W/kk)
di "T3_jp_z=" jpz " jp_t=" jpt " hand_F=" `hF' " reldif=" reldif(jpt, `hF')
if reldif(jpt, `hF') > 1e-12 | jpt <= jpz {
  di as error "T3 FAIL"
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

* T6: standard mode + subgroups + pwsgs under bootdf
mkp
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) subgroups(sg) bootstrap(standard) reps(25) rseed(6) bootdf(clusters)
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

* T7: replay reprints under bootdf
cap qui wooldid
di "T7_replay_rc=" _rc
if _rc != 0 {
  di as error "T7 FAIL"
  local ++fails
}

if `fails' == 0 di as result "BOOTDF_ALL_PASS"
else di as error "BOOTDF_FAILS=`fails'"
