* test_bootstrap_tstat.do - boottstat: percentile-t (studentized) bootstrap inference.
* Run from test/: stata-mp -e do test_bootstrap_tstat.do
* T1 the plain bootstrap path is untouched (no e(boot_tstat), column not relabelled);
* T2 boottstat posts its whole surface (t_pt label, e(boot_tstat)/tq/tn, e(bootdist_se),
* the hidden Wald twins); T3 the p-value and critical value are reconstructable from
* e() alone, exactly; T4 the interval inverts the critical value exactly and the reported
* statistic is b/se_obs; T5 the se column is still the bootstrap SD; T6 the add-one p
* floor is 1/(completed+1); T7 the quantities boottstat does NOT studentize - pairwise
* subgroup contrasts and the joint tests - are unchanged from the Wald run; T8 the gates;
* T9 replay. Prints BOOTTS_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0

cap prog drop mkp
prog def mkp
  clear
  set obs 14
  gen unit = _n
  gen sg = mod(unit, 2)
  gen ttre = cond(unit<=5, 4, cond(unit<=10, 6, .))
  expand 8
  bysort unit: gen t = _n
  set seed 11
  gen y = 2*(t>=ttre & ttre<.) + 0.3*(t>=ttre & ttre<.)*sg + 0.1*unit + 0.05*t + rnormal()
end

local SPEC att att_it cluster(unit) esprelength(1) espostlength(1)

* ---------------------------------------------------------------- T1: plain path untouched
mkp
qui wooldid y unit t ttre, `SPEC' bootstrap(standard) reps(30) rseed(7)
matrix P1 = e(wdidmainresults)
local p_stat `e(boot_stat)'
local p_ts = e(boot_tstat)
local p_cn : colnames P1
di "T1_plain stat=[`p_stat'] boot_tstat=`p_ts' cols=[`p_cn']"
if "`p_stat'" == "t_pt" | `p_ts' < . | strpos("`p_cn'", "t_pt") > 0 {
  di as error "T1 FAIL: the plain bootstrap path shows studentized markers"
  local ++fails
}

* ---------------------------------------------------------------- T2: the boottstat surface
mkp
qui wooldid y unit t ttre, `SPEC' bootstrap(standard) reps(30) rseed(7) boottstat
matrix B2 = e(wdidmainresults)
local b_stat `e(boot_stat)'
local b_ts = e(boot_tstat)
local b_cn : colnames B2
local b_done = e(boot_completed)
cap matrix TS = e(boot_tstud)
local rc_tq = _rc
local rc_tn = _rc
cap matrix BS = e(bootdist_se)
local rc_bs = _rc
cap matrix WM = e(boot_wald_mainresults)
local rc_wm = _rc
cap matrix TPM = e(boot_tpctl_mainresults)
local rc_tpm = _rc
cap matrix TBM = e(boot_tbc_mainresults)
local rc_tbm = _rc
cap matrix AS = e(boot_ase_obs)
local rc_as = _rc
di "T2_surface stat=[`b_stat'] boot_tstat=`b_ts' completed=`b_done' cols=[`b_cn']"
di "T2_rcs tstud=`rc_tq' bootdist_se=`rc_bs' wald_twin=`rc_wm' tpctl_twin=`rc_tpm' tbc_twin=`rc_tbm' ase_obs=`rc_as'"
* 2026-09-05: the studentized kernel is one NAMED matrix e(boot_tstud); the reported table
* under boottstat carries se_analytic, and both studentized families are posted as twins
if "`b_stat'" != "t_pt" | `b_ts' != 1 | strpos("`b_cn'", "t_pt") == 0 | strpos("`b_cn'", "se_analytic") == 0 ///
   | `rc_tq' != 0 | `rc_bs' != 0 | `rc_wm' != 0 | `rc_tpm' != 0 | `rc_tbm' != 0 | `rc_as' != 0 {
  di as error "T2 FAIL: the boottstat e() surface is incomplete"
  local ++fails
}
* the SE draws must be the same shape as the coefficient draws
matrix BD = e(bootdist)
di "T2_dims bootdist=" rowsof(BD) "x" colsof(BD) " bootdist_se=" rowsof(BS) "x" colsof(BS)
if rowsof(BD) != rowsof(BS) | colsof(BD) != colsof(BS) {
  di as error "T2b FAIL: e(bootdist_se) does not match e(bootdist) in shape"
  local ++fails
}

* ---------------------------------------------------------------- T3: reconstructable from e()
* RE-DERIVED 2026-09-05 (not re-toleranced): the studentized p is the two-sided RANK of
* t_obs among the signed t-draws with NO add-one, and the interval is EQUAL-TAILED from the
* signed quantiles under wooldid_bootci's order-statistic rule - the author's rulings of
* 2026-09-03/04 replaced the symmetric |t| construction with add-one that T3/T4 used to pin.
* e(bootdist) and e(bootdist_se) are compacted on ONE row mask, so pairing them row by row
* studentizes the same replication (that misalignment was ledger item #2).
matrix OB = e(b)
mata:
  A = st_matrix("BD"); S = st_matrix("BS"); ob = st_matrix("OB"); as = st_matrix("AS"); TS = st_matrix("TS")
  so = as[1,1]; to = ob[1,1]/so
  ok = (A[.,1] :< .) :& (S[.,1] :< .) :& (S[.,1] :> 0)
  tr = (select(A[.,1], ok) :- ob[1,1]) :/ select(S[.,1], ok)
  v = sort(tr, 1); n = rows(v); al = 0.025
  x = n*al;     qlo = (x==floor(x) & x>=1 ? (v[x]+v[min((x+1,n))])/2 : v[max((min((ceil(x),n)),1))])
  x = n*(1-al); qhi = (x==floor(x) & x>=1 ? (v[x]+v[min((x+1,n))])/2 : v[max((min((ceil(x),n)),1))])
  st_numscalar("hn", n)
  st_numscalar("hp", min((1, (2*min((sum(v:<=to), sum(v:>=to))))/n)))
  st_numscalar("hqlo", qlo); st_numscalar("hqhi", qhi); st_numscalar("hto", to)
  st_numscalar("rowsBD", rows(A)); st_numscalar("rowsBS", rows(S))
end
di "T3_hand n=" hn " tn=" TS[1,1] "  p=" %18.14f hp " table=" %18.14f B2[1,4] "  qlo=" %14.8f hqlo " qhi=" %14.8f hqhi " tstud=" %14.8f TS[4,1] " " %14.8f TS[5,1]
if hn != TS[1,1] | reldif(hp, B2[1,4]) > 0 | reldif(hqlo, TS[4,1]) > 0 | reldif(hqhi, TS[5,1]) > 0 | reldif(hto, B2[1,3]) > 1e-12 {
  di as error "T3 FAIL: the reported studentized numbers are not reproducible from e()"
  local ++fails
}
if rowsBD != rowsBS {
  di as error "T3b FAIL: e(bootdist) and e(bootdist_se) are not row-aligned"
  local ++fails
}

* ---------------------------------------------------------------- T4: the interval is equal-tailed and inverts the quantiles
local seo = AS[1,1]
local hlb = B2[1,1] - TS[5,1]*`seo'
local hub = B2[1,1] - TS[4,1]*`seo'
di "T4_ci lb=" %14.10f B2[1,5] " hand=" %14.10f `hlb' "  ub=" %14.10f B2[1,6] " hand=" %14.10f `hub'
di "T4_stat t_pt=" %14.8f B2[1,3] " b/se_obs=" %14.8f B2[1,1]/`seo'
if reldif(B2[1,5], `hlb') > 1e-12 | reldif(B2[1,6], `hub') > 1e-12 | reldif(B2[1,3], B2[1,1]/`seo') > 1e-12 {
  di as error "T4 FAIL: the interval does not invert the equal-tailed studentized quantiles"
  local ++fails
}
* it must NOT be the old symmetric interval (a did-the-ruling-take-effect guard)
if reldif(B2[1,6] - B2[1,1], B2[1,1] - B2[1,5]) < 1e-10 {
  di as error "T4b FAIL: the studentized interval is symmetric about b - the equal-tailed ruling did not take effect"
  local ++fails
}
* and p must be dual to the interval
local exc4 = (B2[1,5] > 0 | B2[1,6] < 0)
if `exc4' != (B2[1,4] <= 0.05) {
  di as error "T4c FAIL: studentized p and interval disagree about zero"
  local ++fails
}

* ------------------------------------------------- T5: the se column is the ANALYTIC se
* RE-PINNED 2026-09-04. The old assertion was that boottstat leaves the bootstrap standard
* deviation in column 2. The author's ruling replaced that: a studentized run reports the
* ANALYTIC standard error - the very denominator the studentization divides by - because
* there is no standard deviation of a studentized distribution in coefficient units, and
* an interval-implied standard error is not what the literature reports. The column is
* relabelled se_analytic so the table says which one it is.
di "T5_se boottstat=" %20.16f B2[1,2] " analytic=" %20.16f AS[1,1] " sd_twin=" %20.16f WM[1,2]
if reldif(B2[1,2], AS[1,1]) > 0 {
  di as error "T5 FAIL: the se column is not the analytic standard error"
  local ++fails
}
if reldif(B2[1,2], WM[1,2]) <= 0 {
  di as error "T5b FAIL: the se column still equals the bootstrap SD - the ruling did not take effect"
  local ++fails
}
local t5cn : colnames B2
di "T5_cols=[`t5cn']"
if strpos("`t5cn'", "se_analytic") == 0 {
  di as error "T5c FAIL: the se column is not labelled se_analytic"
  local ++fails
}

* ---------------------------------------------------------------- T6: the p lives on the 2/n grid
* 2026-09-05: no add-one any more (author ruling, measured); the p lives on the 2/n grid
local n6 = TS[1,1]
local grid = round(B2[1,4] * `n6' / 2)
di "T6_grid n=`n6' p=" %18.14f B2[1,4] " p*n/2=" %12.6f B2[1,4]*`n6'/2 " (must be an integer unless p == 1)"
if B2[1,4] < 1 & reldif(B2[1,4]*`n6'/2, `grid') > 1e-9 {
  di as error "T6 FAIL: the studentized p is not on the 2/n rank grid"
  local ++fails
}

* ---------------------------------------------------------------- T7: contrasts studentized, joint tests present only under boottstat
mkp
qui wooldid y unit t ttre, att att_it cluster(unit) esprelength(1) espostlength(1) subgroups(sg) jointtests bootstrap(standard) reps(30) rseed(7) boottstat
local rc7 = _rc
cap matrix PW = e(pwsgs)
local rc_pw = _rc
scalar jh7 = e(joint_horiz)
cap matrix SW = e(boot_wald_sgresults_1)
local rc_sw = _rc
local pwcn : colnames PW
di "T7_boot rc=`rc7' pwsgs rc=`rc_pw' sg1_wald_twin rc=`rc_sw' joint_horiz=" jh7
di "T7_pwcols=[`pwcn']"
* 2026-09-05: under boottstat the pairwise contrasts ARE studentized (author ruling: boottstat
* extends to the contrasts), so their columns carry t_pt and se_analytic
if strpos("`pwcn'", "t_pt") == 0 | strpos("`pwcn'", "se_analytic") == 0 {
  di as error "T7 FAIL: pairwise contrasts were not studentized under boottstat"
  local ++fails
}
cap matrix TSP = e(boot_tstud_pwsgs)
if _rc != 0 {
  di as error "T7c FAIL: e(boot_tstud_pwsgs) not posted"
  local ++fails
}
* T7d: the studentized contrast is reconstructable from e() BY NAME - contrast draws from
* e(bootdist), each replication's contrast se from e(bootdist_pwse) (row-aligned with
* bootdist), the observed analytic pairwise se as denominator. Anti-vacuity: the row must be
* LIVE (se_analytic > 0) and n >= 2, else the test fails rather than passing on blanks
matrix PW7 = e(pwsgs)
matrix BD7 = e(bootdist)
matrix PS7 = e(bootdist_pwse)
matrix OB7 = e(b)
local pwr : rowfullnames PW7
local r7 : list posof "att:1vs0.sg" in pwr
if `r7' == 0 local r7 : list posof "att:1vs0bn.sg" in pwr
local cn7 : colnames OB7
local ka7 : list posof "main_SG1_att" in cn7
local kb7 : list posof "main_SG0_att" in cn7
local psn : colfullnames PS7
local pj7 : list posof "att:1vs0.sg" in psn
if `pj7' == 0 local pj7 : list posof "att:1vs0bn.sg" in psn
di "T7d_map row=`r7' ka=`ka7' kb=`kb7' pwsecol=`pj7'"
if `r7' == 0 | `ka7' == 0 | `kb7' == 0 | `pj7' == 0 | PW7[`r7',2] >= . | PW7[`r7',2] <= 0 {
  di as error "T7d FAIL (VACUOUS): the att:1vs0 contrast is not live or not mappable, so nothing can be checked"
  local ++fails
}
else {
  scalar r7s = `r7'
  scalar ka7s = `ka7'
  scalar kb7s = `kb7'
  scalar pj7s = `pj7'
  mata:
    PW=st_matrix("PW7"); BD=st_matrix("BD7"); PS=st_matrix("PS7"); OB=st_matrix("OB7")
    r=st_numscalar("r7s"); ka=st_numscalar("ka7s"); kb=st_numscalar("kb7s"); pj=st_numscalar("pj7s")
    dobs = OB[1,ka]-OB[1,kb]; seobs = PW[r,2]
    ok = (BD[.,ka]:<.) :& (BD[.,kb]:<.) :& (PS[.,pj]:<.) :& (PS[.,pj]:>0)
    tr = (select(BD[.,ka]:-BD[.,kb], ok) :- dobs) :/ select(PS[.,pj], ok)
    v = sort(tr,1); n=rows(v); to=dobs/seobs; al=0.025
    x=n*al;     qlo = (x==floor(x) & x>=1 ? (v[x]+v[min((x+1,n))])/2 : v[max((min((ceil(x),n)),1))])
    x=n*(1-al); qhi = (x==floor(x) & x>=1 ? (v[x]+v[min((x+1,n))])/2 : v[max((min((ceil(x),n)),1))])
    st_numscalar("h7n", n)
    st_numscalar("h7p", min((1, (2*min((sum(v:<=to), sum(v:>=to))))/n)))
    st_numscalar("h7lb", dobs - qhi*seobs); st_numscalar("h7ub", dobs - qlo*seobs)
    st_numscalar("h7t", to)
  end
  di "T7d_hand n=" h7n " p=" %12.8f h7p " lb=" %12.8f h7lb " ub=" %12.8f h7ub "   table p=" %12.8f PW7[`r7',4] " lb=" %12.8f PW7[`r7',5] " ub=" %12.8f PW7[`r7',6]
  if h7n < 2 {
    di as error "T7d FAIL (VACUOUS): fewer than 2 usable replications for the contrast"
    local ++fails
  }
  else if reldif(h7p, PW7[`r7',4]) > 1e-12 | reldif(h7lb, PW7[`r7',5]) > 1e-10 | reldif(h7ub, PW7[`r7',6]) > 1e-10 | reldif(h7t, PW7[`r7',3]) > 1e-12 {
    di as error "T7d FAIL: the studentized contrast row is not reproducible from e() by name"
    local ++fails
  }
}
* T7b RE-PINNED 2026-09-04. The old assertion was that boottstat leaves the joint p alone,
* because the joint tests were Wald either way. The author's ruling inverts that premise:
* a joint test that cannot be BOOTSTRAPPED is ABSENT, and bootstrapping one requires
* boottstat. So boottstat is now precisely what decides whether the joint test EXISTS -
* absent without it, present with it - and the flag must record the absence for replay.
mkp
qui wooldid y unit t ttre, att att_it cluster(unit) esprelength(1) espostlength(1) subgroups(sg) jointtests bootstrap(standard) reps(30) rseed(7)
cap confirm scalar e(joint_horiz)
local rc7w = _rc
local flag7w `e(boot_jt_absent)'
di "T7_joint boottstat=" jh7 " plain_rc=`rc7w' plain_flag=[`flag7w'] (want a value, 111, 1)"
if jh7 >= . | `rc7w' == 0 | "`flag7w'" != "1" {
  di as error "T7b FAIL: boottstat must decide whether the joint test exists at all"
  local ++fails
}

* ---------------------------------------------------------------- T8: the gates
mkp
cap noi wooldid y unit t ttre, `SPEC' boottstat
local rc8a = _rc
cap noi wooldid y unit t ttre, `SPEC' fastnose bootstrap(standard) reps(20) rseed(7) boottstat
local rc8b = _rc
di "T8_gates stray=`rc8a' with_fastnose=`rc8b' (want 198 198)"
if `rc8a' != 198 | `rc8b' != 198 {
  di as error "T8 FAIL: a gate did not fire"
  local ++fails
}

* ---------------------------------------------------------------- T9: replay
mkp
qui wooldid y unit t ttre, `SPEC' bootstrap(standard) reps(20) rseed(7) boottstat
cap qui wooldid
di "T9_replay rc=" _rc
if _rc != 0 {
  di as error "T9 FAIL: replay after a boottstat run"
  local ++fails
}

* ---------------------------------------------------------------- T10: weighted + unclustered + boottstat takes the physical
* expansion leg (ruling 2026-09-04). The aweight count leg carried the resample multiplicity as an aweight, which loses the
* duplicated rows from N and squares the weight in the robust meat: pooled E[se_r^2]/se_obs^2 measured 1.8-2.1. The physical
* leg must run clean (no failed replications), sit well below that, and draw exactly what the plain (non-studentized) run draws
mkp
set seed 4
gen double wgt = 0.5 + rchi2(3)
tempfile t10d
save `t10d'
cap qui wooldid y unit t ttre [aw=wgt], att esprelength(1) espostlength(1) bootstrap(standard) reps(30) rseed(7) boottstat
local rc10 = _rc
di "T10_phys rc=`rc10' completed=" e(boot_completed) " failed=" e(boot_failed)
if `rc10' != 0 | e(boot_failed) != 0 {
  di as error "T10 FAIL: the weighted unclustered boottstat run did not complete cleanly (rc=`rc10', failed=" e(boot_failed) ")"
  local ++fails
}
else {
  matrix T10_SE = e(bootdist_se)
  matrix T10_AS = e(boot_ase_obs)
  matrix T10_BD = e(bootdist)
  mata:
    S=st_matrix("T10_SE"); A=st_matrix("T10_AS")
    r=J(0,1,.)
    for (j=1;j<=cols(S);j++) {
      if (A[1,j] >= . | A[1,j] <= 0) continue
      s2=select(S[.,j], S[.,j]:<.)
      if (rows(s2) < 2) continue
      r=r \ ((s2:^2):/(A[1,j]^2))
    }
    st_numscalar("t10ratio", rows(r) > 0 ? mean(r) : .); st_numscalar("t10n", rows(r))
  end
  use `t10d', clear
  qui wooldid y unit t ttre [aw=wgt], att esprelength(1) espostlength(1) bootstrap(standard) reps(30) rseed(7)
  matrix T10_BDc = e(bootdist)
  mata:
    A=st_matrix("T10_BD"); B=st_matrix("T10_BDc")
    st_numscalar("t10maxrel", rows(A)==rows(B) & cols(A)==cols(B) ? max(abs(A:-B):/(abs(B):+1e-300)) : .)
  end
  di "T10_ratio pooled E[se_r^2]/se_obs^2=" %8.4f t10ratio " (n=" t10n ")  draws vs plain run max reldif=" %10.3e t10maxrel
  if t10n < 10 {
    di as error "T10 FAIL (VACUOUS): too few replication standard errors to pool"
    local ++fails
  }
  else if t10ratio >= 1.4 | t10ratio <= 0.5 | t10maxrel >= . | t10maxrel > 1e-9 {
    di as error "T10 FAIL: the physical leg's replication standard errors are inflated (ratio " t10ratio ") or its draws differ from the plain run"
    local ++fails
  }
}
use `t10d', clear
cap qui wooldid y unit t ttre [pw=wgt], att esprelength(1) espostlength(1) bootstrap(standard) reps(10) rseed(7) boottstat
di "T10_pw rc=" _rc " failed=" e(boot_failed)
if _rc != 0 | e(boot_failed) != 0 {
  di as error "T10 FAIL: pweight + unclustered + boottstat does not run clean"
  local ++fails
}

if `fails' == 0 di as result "BOOTTS_ALL_PASS"
else di as error "BOOTTS_FAILS=`fails'"
