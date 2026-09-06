* test_bootstrap_joint.do - commit-E validation: bootstrapped joint tests, subgroup
* contrasts, and CI flavors. Run from test/: stata-mp -e do test_bootstrap_joint.do
* ET1 joint_horiz == hand Wald from e(bootdist) (exact); ET2 SI-A1 carve-out (analytic
* refuses r198, bootstrap runs with a real p); ET3 pwsgs z-column fill + allsgssame
* boot Wald; ET4 percentile AND bias-corrected CIs EXACT vs the bstat oracle fed our
* own replication matrix + rep-0 comparator (the odd-rep-count case that separates the
* classic percentile rule from interpolation); ET5 no joint scalars without jointtests;
* ET6 analytic joint path unchanged off-state; ET7 ties z0 convention == bstat.
* Prints BOOTJOINT_ALL_PASS iff all pass.
clear all
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0
cap prog drop mkpan
prog def mkpan
  clear
  set obs 12
  gen unit = _n
  gen s = ceil(unit/4)
  gen sg = mod(unit, 2)
  gen ttre = cond(unit<=4, 4, cond(unit<=8, 6, .))
  expand 8
  bysort unit: gen t = _n
  set seed 11
  gen y = 2*(t>=ttre & ttre<.) + 0.3*(t>=ttre & ttre<.)*sg + 0.1*unit + 0.05*t + rnormal()
end

* ---- E-T1: joint_horiz under bootstrap == hand Wald from bootdist (exact)
* the reference distribution is requested explicitly: this leg checks the Wald construction
* against chi2, and since the ruling of 2026-09-02 the default on a CLUSTERED run is
* bootdf(clusters), which reports F(k, G-1) instead. bootdf() itself is covered in
* test_bootstrap_df.do (T1 both branches, T3 the F reference).
* RE-PINNED 2026-09-04. Author ruling: a joint test that cannot be BOOTSTRAPPED is ABSENT,
* and bootstrapping one means bootstrapping its test statistic, which requires boottstat.
* So the no-boottstat run must post no joint scalar at all, and the hand-Wald identity is
* checked on the boottstat run. NOTE FOR LAYER 2: under boottstat the statistic is slated
* to become the recentered studentized W_r rather than the draw-covariance Wald checked
* here, at which point the ET1 identity below must be re-derived, not merely re-toleranced.
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(1) jointtests bootstrap(bayes) reps(40) rseed(99) bootdf(normal)
* RE-PINNED 2026-09-06 (F034): without boottstat the joint test is PRESENT as the
* fixed-metric rank test (e(boot_jt_family) = fixedmetric) and reconstructable by name from
* e(bootdist) and e(b): S = listwise covariance of the pre-horizon draws, W_obs = b' S^-1 b,
* W_r = (b_r - b)' S^-1 (b_r - b), p = (1 + #{W_r >= W_obs}) / (n + 1); e(boot_jointWobs)
* holds W_obs and e(boot_jointW) the W_r
cap confirm scalar e(joint_horiz)
local et1_present_rc = _rc
local et1_fam `e(boot_jt_family)'
di "ET1_FIXEDMETRIC_WITHOUT_BOOTTSTAT rc=`et1_present_rc' family=[`et1_fam'] p=" e(joint_horiz)
if `et1_present_rc' != 0 | "`et1_fam'" != "fixedmetric" | "`e(boot_jt_absent)'" != "" {
  di as error "ET1z FAIL: the fixed-metric joint test must be present and labelled without boottstat"
  local ++fails
}
else {
  matrix ET1F_B = e(b)
  matrix ET1F_A = e(bootdist)
  matrix ET1F_WO = e(boot_jointWobs)
  * the accumulator columns follow e(wdidmainresults)' rows (the ado's positional walk):
  * the ES_att pre-horizon columns are the rows with equation ES_att and relyr (col 7) < 0
  matrix ET1F_M = e(wdidmainresults)
  local et1frn : rowfullnames ET1F_M
  local et1fpre
  local et1fk 0
  foreach rn of local et1frn {
    local ++et1fk
    if strpos("`rn'", "ES_att:") == 1 & ET1F_M[`et1fk', 7] < 0 local et1fpre `et1fpre' `et1fk'
  }
  mata: idx = strtoreal(tokens(st_local("et1fpre"))); A = st_matrix("ET1F_A")[., idx']; b0 = st_matrix("ET1F_B")[1, idx']; ok = rowmissing(A) :== 0; D = select(A, ok); n = rows(D); S = variance(D); Si = invsym(S); Wo = b0 * Si * b0'; Wr = J(n, 1, .); for (r = 1; r <= n; r++) Wr[r] = (D[r, .] - b0) * Si * (D[r, .] - b0)'; st_numscalar("et1f_p", (1 + sum(Wr :>= Wo)) / (n + 1)); st_numscalar("et1f_wo", Wo); st_numscalar("et1f_n", n)
  di "ET1_FIXEDMETRIC hand p=" %18.14f et1f_p " table=" %18.14f e(joint_horiz) " W_obs hand=" %14.8f et1f_wo " posted=" %14.8f ET1F_WO[1,1] " n=" et1f_n " cols=[`et1fpre']"
  if abs(et1f_p - e(joint_horiz)) > 1e-12 | reldif(et1f_wo, ET1F_WO[1,1]) > 1e-10 | et1f_n < 2 {
    di as error "ET1y FAIL: the fixed-metric joint test is not reconstructable by name from e(bootdist)/e(b)"
    local ++fails
  }
}
* RE-DERIVED 2026-09-05 (Layer 2): under boottstat the joint statistic is the bootstrap of the
* TEST STATISTIC - each replication's recentred Wald (b_r - b)' V_r^-1 (b_r - b) on its own
* e(V_es_att), stored by name in e(boot_jointW) - and the p is the add-one share at or beyond
* the observed Wald on the analytic e(V_es_analytic_att). boottstat is refused with
* bootstrap(bayes) (task 2), so the identity runs on bootstrap(standard)
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(1) jointtests bootstrap(standard) reps(40) rseed(99) bootdf(normal) boottstat
scalar jp = e(joint_horiz)
matrix JW = e(boot_jointW)
matrix VA = e(V_es_analytic_att)
matrix BE = e(b_es_att)
local cn : colnames BE
local pcols
local k 0
foreach nm of local cn {
  local ++k
  if substr("`nm'",1,3) == "pre" local pcols `pcols' `k'
}
local jwn : colnames JW
local jcol : list posof "att" in jwn
scalar jc = `jcol'
di "ET1_map pcols=[`pcols'] jointW col=`jcol' of [`jwn']"
if `jcol' == 0 | "`pcols'" == "" {
  di as error "ET1 FAIL (VACUOUS): e(boot_jointW) has no att column or e(b_es_att) has no pre horizons"
  local ++fails
}
else {
  mata: JW = st_matrix("JW"); VA = st_matrix("VA"); BE = st_matrix("BE"); idx = strtoreal(tokens(st_local("pcols")))
  mata: b = BE[1, idx']'; V = VA[idx', idx']; V = (V + V')/2
  mata: sc = sqrt(diagonal(V))'; R = V :/ (sc' * sc); Si = invsym(R) :/ (sc' * sc)
  mata: Wobs = b' * Si * b; w = JW[., st_numscalar("jc")]; w = select(w, w :< .)
  mata: st_numscalar("jphand", (1 + sum(w :>= Wobs))/(rows(w) + 1)); st_numscalar("jpn", rows(w))
  di "ET1_JOINT_RELDIF=" reldif(jp, jphand) " usable=" jpn
  if jpn < 2 {
    di as error "ET1 FAIL (VACUOUS): fewer than 2 usable replication statistics"
    local ++fails
  }
  else if reldif(jp, jphand) > 1e-12 {
    di as error "ET1 FAIL"
    local ++fails
  }
}

* ---- E-T2: SI-A1 carve-out
mkpan
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(1) jointtests
local t2a = _rc
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(1) jointtests bootstrap(standard) reps(20) rseed(4) boottstat
local t2b = _rc
di "ET2_ANALYTIC_RC=" `t2a' " BOOT_RC=" `t2b' " jp=" e(joint_horiz)
* the carve-out is unchanged - analytic refuses r(198), the bootstrap route does not - but
* since the 2026-09-04 ruling the bootstrap route only PRODUCES a joint test under boottstat
if `t2a' != 198 | `t2b' != 0 | e(joint_horiz) >= . {
  di as error "ET2 FAIL"
  local ++fails
}

* ---- E-T3: pwsgs filled + allsgssame boot Wald + per-sg macro split
mkpan
* RE-POINTED 2026-09-04: "z == est/se" is a property of the STANDARD-DEVIATION family, so
* this run asks for it explicitly; and the across-subgroup equality scalar now requires
* boottstat, as the joint tests do
* SPLIT 2026-09-04: this test was checking two different things at once, and the author's
* ruling made their preconditions mutually exclusive - "z == est/se" is a property of the
* WALD display, while the across-subgroup equality scalar now requires boottstat, and
* bootdisplay(wald) + boottstat is refused because a studentized run has no bootstrap
* standard deviation to build a normal/t approximation on. So: run A for the contrast
* algebra, run B for the equality scalar.
* run B RE-POINTED 2026-09-05: boottstat is refused with bootstrap(bayes), so it runs on
* bootstrap(standard); unclustered and att-only, because att_it needs cluster() under
* bootstrap(standard) and the analytic pairwise se the studentized equality test needs is
* degenerate under clustering at or above i on some designs (pre-existing at HEAD)
qui wooldid y unit t ttre, att esprelength(1) espostlength(1) subgroups(sg) bootstrap(standard) reps(30) rseed(7) safety(off) boottstat
scalar allsg_B = e(joint_att_allsgssame)
qui wooldid y unit t ttre, att att_it cluster(unit) esprelength(1) espostlength(1) subgroups(sg) bootstrap(bayes) reps(30) rseed(7) safety(off) bootdisplay(wald)
matrix P = e(pwsgs)
local pc : colnames P
di "ET3_PWCOLS=[`pc'] se11=" P[1,2] " allsg(boottstat run)=" allsg_B
if P[1,2] >= . | P[1,2] <= 0 | allsg_B >= . {
  di as error "ET3 FAIL"
  local ++fails
}
if reldif(P[1,3], P[1,1]/P[1,2]) > 1e-12 {
  di as error "ET3b FAIL: z != est/se"
  local ++fails
}

* ---- E-T4: CI trio vs bstat oracle on OUR bootdist + comparator
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) bootstrap(bayes) reps(35) rseed(12) level(90)
matrix BD = e(bootdist)
matrix CMP = e(boot_comparator)
matrix CIP = e(boot_ci_percentile)
matrix CIB = e(boot_ci_bc)
preserve
clear
qui svmat double BD, names(bsv)
local nb = colsof(BD)
forval c = 1/`nb' {
  rename bsv`c' _bs_`c'
}
local cn : colfullnames BD
local k 0
foreach nm of local cn {
  local ++k
  label var _bs_`k' "[x]`nm'"
}
matrix OBSST = CMP
matrix colnames OBSST = `cn'
qui bstat _bs_*, stat(OBSST) n(96) level(90)
matrix OCIP = e(ci_percentile)
matrix OCIB = e(ci_bc)
restore
di "ET4_CIP_DIFF=" mreldif(CIP, OCIP) " CIB_DIFF=" mreldif(CIB, OCIB)
if mreldif(CIP, OCIP) > 1e-12 | mreldif(CIB, OCIB) > 1e-12 {
  di as error "ET4 FAIL"
  local ++fails
}

* ---- E-T5: no-jointtests bootstrap run posts NO joint scalars
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(1) bootstrap(bayes) reps(15) rseed(5)
if "`e(joint_horiz)'" != "" {
  di as error "ET5 FAIL: joint posted without jointtests"
  local ++fails
}
di "ET5_OK"

* ---- E-T6: off-state analytic joint tests unchanged
mkpan
qui wooldid y unit t ttre, att robust esprelength(3) espostlength(1) jointtests
di "ET6_ANALYTIC_JP=" e(joint_horiz)
if e(joint_horiz) >= . {
  di as error "ET6 FAIL"
  local ++fails
}

* ---- E-T8 (verifier F1 regression): with THREE OR MORE subgroups Stata marks the
* base level in the pairwise rownames (1vs0bn.sg, not 1vs0.sg). A mapping that does
* not tolerate the marker leaves those rows carrying ANALYTIC se/t/p/CI inside a
* table the display attributes to the bootstrap covariance - measured at 46x too
* small, z=-67, p=2e-33. No row may retain the analytic value, and any row that
* cannot be mapped must be blank rather than analytic.
use "`REPO'/test/data/dta_path_relyrgap.dta", clear
cap drop sg3
gen sg3 = mod(i, 3)
tempfile d3sg
save `d3sg'
qui wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(2) jointtests subgroups(sg3) safety(off)
matrix PA = e(pwsgs)
use `d3sg', clear
qui wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(2) jointtests subgroups(sg3) safety(off) bootstrap(bayes) reps(30) rseed(21)
matrix PB = e(pwsgs)
local nr = rowsof(PB)
local et8leak 0
local et8halfblank 0
forval r = 1/`nr' {
  if PA[`r',2] < . & PB[`r',2] < . {
    if reldif(PA[`r',2], PB[`r',2]) == 0 local ++et8leak
  }
  if PB[`r',2] >= . & PB[`r',4] < . local ++et8halfblank
}
di "ET8_rows=" `nr' " analytic_leaks=" `et8leak' " half_blank=" `et8halfblank'
if `et8leak' > 0 | `et8halfblank' > 0 {
  di as error "ET8 FAIL: analytic inference survived into e(pwsgs)"
  local ++fails
}

* ---- E-T7: ties z0 convention pinned vs bstat (half-weighted ties)
clear
set seed 77
set obs 30
gen double s1 = round(rnormal(2, .5), .25)
* comparator value planted to tie with several reps
count if s1 == 2
local nties = r(N)
di "PLANTED_TIES=" `nties'
rename s1 _bs_1
label var _bs_1 "[a]x"
matrix OBS = (2)
matrix colnames OBS = a:x
qui bstat _bs_1, stat(OBS) n(50) level(90)
scalar z0plain = el(e(z0),1,1)
cap qui bstat _bs_1, stat(OBS) n(50) level(90) ties
di "BSTAT_TIES_RC=" _rc
scalar z0ties = el(e(z0),1,1)
* our formulas
qui sum _bs_1
mata: v = st_data(., "_bs_1"); c = 2
mata: st_numscalar("ourplain", invnormal(sum(v :<= c)/rows(v)))
mata: st_numscalar("ourties", invnormal((sum(v :< c) + 0.5*sum(v :== c))/rows(v)))
di "ET7_TIESPIN plain: bstat=" z0plain " ours=" ourplain " diff=" reldif(z0plain, ourplain)
di "ET7_TIESPIN ties:  bstat=" z0ties " ours=" ourties " diff=" reldif(z0ties, ourties)

if reldif(z0plain, ourplain) > 0 | reldif(z0ties, ourties) > 0 {
  di as error "ET7 FAIL: ties z0 convention diverges from bstat"
  local ++fails
}

* ---- E-T9: NO ANALYTIC across-subgroup equality p-value survives a bootstrap run.
* REGRESSION TEST for a leak found by audit on 2026-09-04 and fixed the same day. The
* blanking loop in wooldid_bootjnt iterated a STATIC six-name item list; the worker also
* posts contreat twins via subinstr(item,"att","ame"), which yields ame_i and ame_itime -
* both missing from that list. So under bootstrap(bayes) + contreat() + subgroups() with
* att_i/att_itime, four ANALYTIC equality p-values printed while the footer asserted they
* were absent and e(boot_sgjt_absent) was set. The disclosure contradicted the output.
*
* This test is written to FAIL LOUDLY RATHER THAN VACUOUSLY: it first checks that the
* ANALYTIC baseline actually produces the ame_* scalars, because a fixture that does not
* would let the bootstrap leg pass while checking nothing. Four earlier attempts at this
* fixture each passed vacuously for a different reason.
*
* NB: `clear', NOT `clear all'. `clear all' drops loaded PROGRAMS, which would discard the
* working-tree wooldid this file `do'-loaded at the top and silently fall through to the
* INSTALLED copy on the adopath - the test would then be exercising a different version of
* the package than the one under test, and its rc would be meaningless.
clear
set more off
set seed 777
set sortseed 777
set obs 12
gen grp = _n
gen ttre = .
replace ttre = 4 if mod(grp,3) == 0
replace ttre = 6 if mod(grp,3) == 1
gen gfe = rnormal()
expand 4
bysort grp: gen person = _n + 100*grp
gen sg = 1 + mod(_n,2)
gen pfe = rnormal()
expand 8
bysort person: gen t = _n
gen w = 0.5 + rchi2(3)
gen c = 0.3 + 1.5*runiform()
gen post = (t >= ttre & ttre != .)
gen eff = cond(sg==1, 1.0, 0.5)
gen y = 2 + gfe + pfe + 0.04*(t-1) + eff*post + 0.6*c*post + rnormal(0, 0.6)
tempfile fxeq
save `fxeq'
local EQB att att_i att_itime cluster(person) subgroups(sg) contreat(c) ///
    contreatcontrols(basic) esprelength(2) espostlength(2) safety(off)
use `fxeq', clear
cap qui wooldid y grp t ttre [aw=w], `EQB'
local eq_rc0 = _rc
local eq_present 0
foreach nm in ame ame_i ame_itime pre_ame pre_ame_i pre_ame_itime att att_i att_itime pre_att pre_att_i pre_att_itime {
  local eqa_`nm' = e(joint_`nm'_allsgssame)
  if e(joint_`nm'_allsgssame) < . & strpos("`nm'", "ame") > 0 local ++eq_present
}
use `fxeq', clear
cap qui wooldid y grp t ttre [aw=w], `EQB' bootstrap(bayes) reps(12) rseed(101) dots(0)
local eq_rc1 = _rc
* RE-PINNED 2026-09-06 (F034): the equality tests are now PRESENT without boottstat as
* fixed-metric rank tests; a "leak" is an equality p that EQUALS its analytic baseline
* value (an add-one share of at most 13 replications cannot coincide with a chi2 tail
* probability except by construction), or an analytic value standing where the bootstrap
* produced none
local eq_leaked 0
local eq_bootpresent 0
foreach nm in ame ame_i ame_itime pre_ame pre_ame_i pre_ame_itime att att_i att_itime pre_att pre_att_i pre_att_itime {
  if e(joint_`nm'_allsgssame) < . {
    local ++eq_bootpresent
    if "`eqa_`nm''" != "" & "`eqa_`nm''" != "." {
      if reldif(e(joint_`nm'_allsgssame), `eqa_`nm'') < 1e-12 {
        local ++eq_leaked
        di as error "ET9 leaked: joint_`nm'_allsgssame = " e(joint_`nm'_allsgssame) " equals the analytic value"
      }
    }
  }
}
di "ET9_rc0=`eq_rc0' rc1=`eq_rc1' analytic_ame_present=`eq_present' boot_present=`eq_bootpresent' leaked=`eq_leaked' family=[`e(boot_jt_family)'] sgjt_absent=[`e(boot_sgjt_absent)']"
if `eq_rc0' != 0 | `eq_rc1' != 0 {
  di as error "ET9 FAIL (VACUOUS): a run errored, so the leak check proved nothing"
  local ++fails
}
else if `eq_present' == 0 {
  di as error "ET9 FAIL (VACUOUS): the analytic baseline posts no ame_* equality scalars, so this fixture cannot exhibit the leak it guards"
  local ++fails
}
else if `eq_leaked' > 0 {
  di as error "ET9 FAIL: `eq_leaked' analytic equality p-value(s) survived a bootstrap run"
  local ++fails
}

* ---------------------------------------------------------------- ET10: the across-subgroup equality test under boottstat is the
* polarised chain Wald (wooldid_bootsgw); with two subgroups it is the squared studentized contrast, so it must be reproducible BY
* NAME from e(bootdist), e(bootdist_pwse), e(pwsgs) (col 2 = se_analytic) and e(b): add-one tail of W_r >= W_obs. Unclustered on
* purpose: under cluster(unit) the analytic pairwise se is degenerate on some designs (pre-existing), which makes this absent.
mkpan
cap qui wooldid y unit t ttre, att esprelength(1) espostlength(1) subgroups(sg) jointtests bootstrap(standard) reps(40) rseed(7) boottstat
local rc10 = _rc
local p10 = e(joint_att_allsgssame)
di "ET10_run rc=`rc10' allsgssame=" `p10'
if `rc10' != 0 | `p10' >= . | `p10' <= 0 | `p10' > 1 {
  di as error "ET10 FAIL: the subgroup-equality test is not posted (in (0,1]) under unclustered boottstat"
  local ++fails
}
else {
  matrix ET10_A = e(bootdist)
  matrix ET10_PS = e(bootdist_pwse)
  matrix ET10_PW = e(pwsgs)
  matrix ET10_OB = e(b)
  local et10cn : colnames ET10_OB
  local et10k0 : list posof "main_SG0_att" in et10cn
  local et10k1 : list posof "main_SG1_att" in et10cn
  local et10pn : colnames ET10_PS
  local et10c1 : word 1 of `et10pn'
  di "ET10_map k0=`et10k0' k1=`et10k1' pwsecol1=[`et10c1'] pwrow1=[`: word 1 of `: rownames ET10_PW'']"
  if `et10k0' == 0 | `et10k1' == 0 | "`et10c1'" != "1vs0.sg" | "`: word 1 of `: rownames ET10_PW''" != "1vs0.sg" {
    di as error "ET10 FAIL (VACUOUS): the subgroup level columns or the contrast column could not be mapped by name"
    local ++fails
  }
  else {
    scalar et10k0 = `et10k0'
    scalar et10k1 = `et10k1'
    mata:
      A=st_matrix("ET10_A"); PS=st_matrix("ET10_PS"); PW=st_matrix("ET10_PW"); OB=st_matrix("ET10_OB")
      k0=st_numscalar("et10k0"); k1=st_numscalar("et10k1")
      dobs=OB[1,k1]-OB[1,k0]; so=PW[1,2]
      d=A[.,k1]-A[.,k0]; s=PS[.,1]
      use=selectindex((d:<.):&(s:<.):&(s:>0))
      W=((d[use]:-dobs):/s[use]):^2
      st_numscalar("et10n", rows(W))
      st_numscalar("et10so", so)
      st_numscalar("et10p", (so < . & so > 0 & rows(W) > 0 ? (1+sum(W:>=(dobs/so)^2))/(rows(W)+1) : .))
    end
    di "ET10_hand usable=" et10n " se_obs=" et10so " p=" %12.8f et10p "  table=" %12.8f `p10' "  reldif=" reldif(et10p, `p10')
    if et10n < 2 | et10so >= . {
      di as error "ET10 FAIL (VACUOUS): fewer than 2 usable replications or no observed pairwise se, nothing was checked"
      local ++fails
    }
    else if reldif(et10p, `p10') > 1e-10 {
      di as error "ET10 FAIL: the subgroup-equality p is not the add-one tail of the polarised replication statistics"
      local ++fails
    }
  }
}

if `fails' == 0 di as result "BOOTJOINT_ALL_PASS"
else di as error "BOOTJOINT_FAILS=`fails'"
