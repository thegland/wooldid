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
mkpan
qui wooldid y unit t ttre, att cluster(unit) esprelength(3) espostlength(1) jointtests bootstrap(bayes) reps(40) rseed(99)
scalar jp = e(joint_horiz)
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
mata: S = (dev' * dev)/(rows(CC)-1)
mata: st_numscalar("jphand", chi2tail(cols(A), b' * invsym(S) * b))
di "ET1_JOINT_RELDIF=" reldif(jp, jphand)
if reldif(jp, jphand) > 1e-12 {
  di as error "ET1 FAIL"
  local ++fails
}

* ---- E-T2: SI-A1 carve-out
mkpan
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(1) jointtests
local t2a = _rc
cap qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(1) jointtests bootstrap(bayes) reps(20) rseed(4)
local t2b = _rc
di "ET2_ANALYTIC_RC=" `t2a' " BOOT_RC=" `t2b' " jp=" e(joint_horiz)
if `t2a' != 198 | `t2b' != 0 | e(joint_horiz) >= . {
  di as error "ET2 FAIL"
  local ++fails
}

* ---- E-T3: pwsgs filled + allsgssame boot Wald + per-sg macro split
mkpan
qui wooldid y unit t ttre, att att_it cluster(unit) esprelength(1) espostlength(1) subgroups(sg) bootstrap(bayes) reps(30) rseed(7) safety(off)
matrix P = e(pwsgs)
local pc : colnames P
di "ET3_PWCOLS=[`pc'] se11=" P[1,2] " allsg=" e(joint_att_allsgssame)
if P[1,2] >= . | P[1,2] <= 0 | e(joint_att_allsgssame) >= . {
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

if `fails' == 0 di as result "BOOTJOINT_ALL_PASS"
else di as error "BOOTJOINT_FAILS=`fails'"
