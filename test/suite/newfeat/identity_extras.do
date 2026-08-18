* identity_extras.do — compact identities for the remaining new features:
*  (a) F4 anticipation(1) == manual onset shift (gen double ttre2 = ttre - 1), identical
*      window options on both sides (the ado applies the shift BEFORE any window or
*      reference-period logic and documents the hand-edit as bit-identical), e(b)/e(V)
*      mreldif < 1e-12; checked with and without an ES window.
*  (b) F8 escomp: recompute nobs / sumw_att (and the count columns, shares, and the
*      ref/trt coverage scalars) directly from the data; rows align 1:1 with ES labels.
*  (c) F7 pretrendpower(80): recompute the GLS pre-trend slope/se/W from the pre-block of
*      e(b_es_att)/e(V_es_att) in mata; verify prep, hbar, mdt_wald via the noncentral-F
*      formula, the NIS size/power via the same mvnormalcvqp quadrature. Internal
*      consistency + formula fidelity, NOT an external package.
*  (d) F9 fastnose parity (one spot config, e(b) at 1e-10) and F3 gteffects aggregation
*      identity (w_att-weighted mean == att at 1e-10).
* Args: 1 = path to wooldid.ado (file or dir); 2 = suite data dir (contains dta_binary.dta).
* Run:  "<stata>" -e do identity_extras.do <adopath> <datadir>
*    -> expect IDENTITY_EXTRAS_ALL_PASS

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 97531
set sortseed 97531

local ADO `1'
if `"`ADO'"' == "" local ADO /tmp/wdsuite_build/newfeat/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'
local DATA `2'
if `"`DATA'"' == "" local DATA /Users/hegland/Documents/Research/wooldid_dev/test/data

global pf 0

**************************** (a) F4 anticipation shift-equivalence ****************************
use `"`DATA'/dta_binary.dta"', clear
qui gen double ttre2 = ttre - 1

* A1: no ES window
cap noi wooldid y i t ttre, att cluster(i) anticipation(1)
if _rc != 0 {
  dis as error "A1 FAIL: anticipation run rc=" _rc
  global pf = $pf + 1
  dis as error "IDENTITY_EXTRAS_FAILS=$pf"
  exit, clear
}
tempname BA VA
matrix `BA' = e(b)
matrix `VA' = e(V)
local na = e(N)
local cna : colnames `BA'
if e(anticipation) != 1 {
  dis as error "A1b FAIL: e(anticipation)=" e(anticipation) " (want 1)"
  global pf = $pf + 1
}
cap noi wooldid y i t ttre2, att cluster(i)
if _rc != 0 {
  dis as error "A1c FAIL: manual-shift run rc=" _rc
  global pf = $pf + 1
}
local cnm : colnames e(b)
if "`cna'" != "`cnm'" {
  dis as error "A1 FAIL: colnames differ: [`cna'] vs [`cnm']"
  global pf = $pf + 1
}
if mreldif(`BA', e(b)) > 1e-12 | mreldif(`VA', e(V)) > 1e-12 | `na' != e(N) {
  dis as error "A1 FAIL: b mreldif=" %12.4e mreldif(`BA', e(b)) " V mreldif=" %12.4e mreldif(`VA', e(V)) " N " `na' " vs " e(N)
  global pf = $pf + 1
}
else dis "A1 info: b mreldif = " %12.4e mreldif(`BA', e(b)) " (0 = bit-identical); V mreldif = " %12.4e mreldif(`VA', e(V))

* A2: ES window (identical window options on both sides)
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod anticipation(1)
if _rc != 0 {
  dis as error "A2 FAIL: anticipation ES run rc=" _rc
  global pf = $pf + 1
  dis as error "IDENTITY_EXTRAS_FAILS=$pf"
  exit, clear
}
tempname BA2 VA2
matrix `BA2' = e(b)
matrix `VA2' = e(V)
local na2 = e(N)
local cna2 : colnames `BA2'
cap noi wooldid y i t ttre2, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod
if _rc != 0 {
  dis as error "A2b FAIL: manual-shift ES run rc=" _rc
  global pf = $pf + 1
}
local cnm2 : colnames e(b)
if "`cna2'" != "`cnm2'" {
  dis as error "A2 FAIL: colnames differ: [`cna2'] vs [`cnm2']"
  global pf = $pf + 1
}
if mreldif(`BA2', e(b)) > 1e-12 | mreldif(`VA2', e(V)) > 1e-12 | `na2' != e(N) {
  dis as error "A2 FAIL: b mreldif=" %12.4e mreldif(`BA2', e(b)) " V mreldif=" %12.4e mreldif(`VA2', e(V)) " N " `na2' " vs " e(N)
  global pf = $pf + 1
}
else dis "A2 info: b mreldif = " %12.4e mreldif(`BA2', e(b)) "; V mreldif = " %12.4e mreldif(`VA2', e(V))

**************************** (b) F8 escomp recomputation ****************************
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod customweightmultiplier(w)
if _rc != 0 {
  dis as error "B0 FAIL: escomp run rc=" _rc
  global pf = $pf + 1
  dis as error "IDENTITY_EXTRAS_FAILS=$pf"
  exit, clear
}
tempname E
matrix `E' = e(escomp)
local ern : rownames `E'
local ebl : colnames e(b_es_att)
local trtN = e(escomp_trtN)
local trtW = e(escomp_trtW)
local refN = e(escomp_refN)
local refW = e(escomp_refW)
local refncoh = e(escomp_refncoh)
local refnunits = e(escomp_refnunits)

* rows align 1:1 with the ES rows (same labels, same order)
if "`ern'" != "`ebl'" {
  dis as error "B1 FAIL: escomp rownames [`ern'] != b_es_att colnames [`ebl']"
  global pf = $pf + 1
}

* estimation-sample tag (fallback to all rows if e(sample) is not posted)
qui gen byte __smp = e(sample)
qui count if __smp
if r(N) == 0 qui replace __smp = 1
qui gen double __rely = t - ttre
qui gen double __wa = w

* per-horizon recomputation from the data: nobs (exact), sumw_att (1e-12), count columns
tempname WSUM
scalar `WSUM' = 0
local nr = rowsof(`E')
forval j = 1/`nr' {
  local r = `E'[`j', 8]
  local lbl : word `j' of `ern'
  * label <-> relyr consistency
  local wantlbl = cond(`r' < 0, "pre" + string(abs(`r')), "post" + string(`r'))
  if "`lbl'" != "`wantlbl'" {
    dis as error "B2 FAIL row `j': label [`lbl'] vs relyr `r' (want `wantlbl')"
    global pf = $pf + 1
  }
  qui count if __smp & ttre < . & __rely == `r'
  if r(N) != `E'[`j', 4] {
    dis as error "B3 FAIL `lbl': nobs recomputed " r(N) " vs escomp " `E'[`j', 4]
    global pf = $pf + 1
  }
  qui sum __wa if __smp & ttre < . & __rely == `r', meanonly
  if reldif(r(sum), `E'[`j', 5]) > 1e-12 {
    dis as error "B4 FAIL `lbl': sumw_att recomputed " %21.0g r(sum) " vs escomp " %21.0g `E'[`j', 5] " (reldif " %12.4e reldif(r(sum), `E'[`j',5]) ")"
    global pf = $pf + 1
  }
  scalar `WSUM' = `WSUM' + `E'[`j', 5]
  qui levelsof ttre if __smp & ttre < . & __rely == `r', local(cohl)
  if `: word count `cohl'' != `E'[`j', 1] {
    dis as error "B5 FAIL `lbl': ncohorts recomputed `: word count `cohl'' vs escomp " `E'[`j', 1]
    global pf = $pf + 1
  }
  qui levelsof i if __smp & ttre < . & __rely == `r', local(unil)
  if `: word count `unil'' != `E'[`j', 2] | `: word count `unil'' != `E'[`j', 3] {
    dis as error "B6 FAIL `lbl': nunits/nclusters recomputed `: word count `unil'' vs escomp " `E'[`j', 2] "/" `E'[`j', 3]
    global pf = $pf + 1
  }
}
* share columns: each row's share of the reported-window weight; the column sums to 1
tempname SSUM
scalar `SSUM' = 0
forval j = 1/`nr' {
  if reldif(`E'[`j', 6], `E'[`j', 5]/`WSUM') > 1e-12 {
    dis as error "B7 FAIL row `j': wshare_att " %21.0g `E'[`j', 6] " vs recomputed " %21.0g `E'[`j', 5]/`WSUM'
    global pf = $pf + 1
  }
  scalar `SSUM' = `SSUM' + `E'[`j', 6]
}
if abs(`SSUM' - 1) > 1e-12 {
  dis as error "B8 FAIL: wshare_att column sums to " %21.0g `SSUM'
  global pf = $pf + 1
}
* coverage scalars: whole treated sample, and the reference pool (relyr == -1, fixed base)
qui count if __smp & ttre < .
if r(N) != `trtN' {
  dis as error "B9 FAIL: trtN recomputed " r(N) " vs " `trtN'
  global pf = $pf + 1
}
qui sum __wa if __smp & ttre < ., meanonly
if reldif(r(sum), `trtW') > 1e-12 {
  dis as error "B10 FAIL: trtW recomputed " %21.0g r(sum) " vs " %21.0g `trtW'
  global pf = $pf + 1
}
qui count if __smp & ttre < . & __rely == -1
if r(N) != `refN' {
  dis as error "B11 FAIL: refN recomputed " r(N) " vs " `refN'
  global pf = $pf + 1
}
qui sum __wa if __smp & ttre < . & __rely == -1, meanonly
if reldif(r(sum), `refW') > 1e-12 {
  dis as error "B12 FAIL: refW recomputed " %21.0g r(sum) " vs " %21.0g `refW'
  global pf = $pf + 1
}
qui levelsof ttre if __smp & ttre < . & __rely == -1, local(rcl)
qui levelsof i if __smp & ttre < . & __rely == -1, local(rul)
if `: word count `rcl'' != `refncoh' | `: word count `rul'' != `refnunits' {
  dis as error "B13 FAIL: refncoh/refnunits recomputed `: word count `rcl''/`: word count `rul'' vs `refncoh'/`refnunits'"
  global pf = $pf + 1
}

**************************** (c) F7 pretrendpower first-principles ****************************
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod pretrendpower(80)
if _rc != 0 {
  dis as error "C0 FAIL: pretrendpower run rc=" _rc
  global pf = $pf + 1
  dis as error "IDENTITY_EXTRAS_FAILS=$pf"
  exit, clear
}
local dfr = e(df_r)
tempname BES VES POW
matrix `BES' = e(b_es_att)
matrix `VES' = e(V_es_att)
matrix `POW' = e(roth_power_att)
local r_slope = e(roth_slope_att)
local r_se    = e(roth_slope_se_att)
local r_ll    = e(roth_slope_ll_att)
local r_ul    = e(roth_slope_ul_att)
local r_prep  = e(roth_prep_att)
local r_hbar  = e(roth_hbar_att)
local r_npre  = e(roth_npre_att)
local r_keff  = e(roth_keff_att)
local r_W     = e(roth_W_att)
local r_nis   = e(roth_nissize_att)
local esrel = -1

* pre-block extraction in b_es column order; m = relyr - reference
local eslbl : colnames `BES'
local k 0
local previdx
local mvals
local j 0
foreach l of local eslbl {
  local ++j
  if substr("`l'", 1, 3) == "pre" {
    local ++k
    local previdx `previdx' `j'
    local mvals `mvals' `=-real(substr("`l'", 4, .)) - (`esrel')'
  }
}
if `k' != `r_npre' {
  dis as error "C1 FAIL: parsed `k' pre labels vs e(roth_npre_att)=" `r_npre'
  global pf = $pf + 1
}
tempname S BP M
matrix `S' = J(`k', `k', 0)
matrix `BP' = J(`k', 1, 0)
matrix `M' = J(`k', 1, 0)
local r 0
foreach jj of local previdx {
  local ++r
  matrix `BP'[`r', 1] = `BES'[1, `jj']
  matrix `M'[`r', 1] = `: word `r' of `mvals''
  local c 0
  foreach jj2 of local previdx {
    local ++c
    matrix `S'[`r', `c'] = `VES'[`jj', `jj2']
  }
}
mata: S = st_matrix(st_local("S")); S = (S + S')/2; Si = invsym(S)
mata: m = st_matrix(st_local("M")); bp = st_matrix(st_local("BP"))
mata: c2 = (m'Si*m)[1,1]
mata: st_numscalar("__c2", c2)
mata: st_numscalar("__slope", (m'Si*bp)[1,1]/c2)
mata: st_numscalar("__se", 1/sqrt(c2))
mata: st_numscalar("__W", (bp'Si*bp)[1,1])
mata: st_numscalar("__keff", rows(S) - sum(diagonal(Si) :== 0))

* GLS slope / se / joint statistic / effective rank
if reldif(__slope, `r_slope') > 1e-8 | reldif(__se, `r_se') > 1e-8 {
  dis as error "C2 FAIL: GLS slope/se recomputed (" %21.0g __slope ", " %21.0g __se ") vs (" %21.0g `r_slope' ", " %21.0g `r_se' ")"
  global pf = $pf + 1
}
else dis "C2 info: slope reldif = " %12.4e reldif(__slope, `r_slope') "; se reldif = " %12.4e reldif(__se, `r_se')
if reldif(__W, `r_W') > 1e-8 {
  dis as error "C3 FAIL: W recomputed " %21.0g __W " vs " %21.0g `r_W'
  global pf = $pf + 1
}
if __keff != `r_keff' {
  dis as error "C3b FAIL: keff recomputed " __keff " vs " `r_keff'
  global pf = $pf + 1
}

* slope CI and joint-test p (this also confirms the residual-df branch: crit/p use e(df_r))
local dfok = (`dfr' > 0 & `dfr' < .)
local crit = cond(`dfok', invttail(`dfr', .025), invnormal(.975))
if reldif(`r_slope' - `crit'*`r_se', `r_ll') > 1e-10 | reldif(`r_slope' + `crit'*`r_se', `r_ul') > 1e-10 {
  dis as error "C4 FAIL: slope CI recomputed (" `r_slope' - `crit'*`r_se' ", " `r_slope' + `crit'*`r_se' ") vs (" `r_ll' ", " `r_ul' ")"
  global pf = $pf + 1
}
local prep_f = cond(`dfok', Ftail(__keff, `dfr', __W/__keff), .)
local prep_c = chi2tail(__keff, __W)
if `dfok' & reldif(`prep_f', `r_prep') <= 1e-10 dis "C5 info: prep matches the F(keff, df_r) branch (df_r = `dfr')"
else if reldif(`prep_c', `r_prep') <= 1e-10 {
  dis as error "C5 FAIL: prep matches the chi2 branch, not the F branch (df plumbing?): F=" `prep_f' " chi2=" `prep_c' " vs " `r_prep'
  global pf = $pf + 1
}
else {
  dis as error "C5 FAIL: prep matches neither branch: F=" `prep_f' " chi2=" `prep_c' " vs " `r_prep'
  global pf = $pf + 1
}

* hbar: att-weighted (here unweighted) mean relative post period over the FULL post support
qui gen byte __smp = e(sample)
qui count if __smp
if r(N) == 0 qui replace __smp = 1
qui gen double __rely = t - ttre
qui sum __rely if __smp & ttre < . & __rely >= 0 & __rely < ., meanonly
if reldif(r(mean), `r_hbar') > 1e-10 {
  dis as error "C6 FAIL: hbar recomputed " %21.0g r(mean) " vs " %21.0g `r_hbar'
  global pf = $pf + 1
}

* mdt_wald via the noncentral-F inverse at the requested power; bias = mdt*(hbar - ref)
if `POW'[1, 1] != 80 {
  dis as error "C7 FAIL: power level " `POW'[1,1] " (want 80)"
  global pf = $pf + 1
}
local lam = cond(`dfok', npnF(__keff, `dfr', invF(__keff, `dfr', .95), 1 - 80/100), npnchi2(__keff, invchi2(__keff, .95), 1 - 80/100))
if reldif(sqrt(`lam'/__c2), `POW'[1, 2]) > 1e-8 {
  dis as error "C7 FAIL: mdt_wald recomputed " %21.0g sqrt(`lam'/__c2) " vs " %21.0g `POW'[1, 2]
  global pf = $pf + 1
}
else dis "C7 info: mdt_wald reldif = " %12.4e reldif(sqrt(`lam'/__c2), `POW'[1,2]) " (lambda* = " %21.0g `lam' ")"
if reldif(`POW'[1, 2]*(`r_hbar' - `esrel'), `POW'[1, 4]) > 1e-10 {
  dis as error "C8 FAIL: bias_wald recomputed " %21.0g `POW'[1,2]*(`r_hbar' - `esrel') " vs " %21.0g `POW'[1, 4]
  global pf = $pf + 1
}

* NIS size and power at the reported mdt_nis, via the same deterministic quadrature
scalar __crit = `crit'
scalar __gnis = `POW'[1, 3]
mata: sd = sqrt(diagonal(S))'
mata: cr = st_numscalar("__crit")
mata: st_numscalar("__nis0", 1 - mvnormalcvqp(-cr*sd, cr*sd, J(1, cols(sd), 0), vech(S)', 1000))
mata: st_numscalar("__pwnis", 1 - mvnormalcvqp(-cr*sd, cr*sd, st_numscalar("__gnis")*(m'), vech(S)', 1000))
if reldif(__nis0, `r_nis') > 1e-10 {
  dis as error "C9 FAIL: NIS size recomputed " %21.0g __nis0 " vs " %21.0g `r_nis'
  global pf = $pf + 1
}
if abs(__pwnis - 0.80) > 2e-4 {
  dis as error "C10 FAIL: NIS rejection probability at reported mdt_nis = " %21.0g __pwnis " (want 0.80 within 2e-4)"
  global pf = $pf + 1
}
else dis "C10 info: NIS power at reported mdt_nis = " %21.0g __pwnis
if reldif(`POW'[1, 3]*(`r_hbar' - `esrel'), `POW'[1, 5]) > 1e-10 {
  dis as error "C11 FAIL: bias_nis recomputed " %21.0g `POW'[1,3]*(`r_hbar' - `esrel') " vs " %21.0g `POW'[1, 5]
  global pf = $pf + 1
}

**************************** (d) F9 fastnose parity + F3 gteffects identity ****************************
use `"`DATA'/dta_binary.dta"', clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
if _rc != 0 {
  dis as error "D0 FAIL: full run rc=" _rc
  global pf = $pf + 1
  dis as error "IDENTITY_EXTRAS_FAILS=$pf"
  exit, clear
}
tempname BF
matrix `BF' = e(b)
local cnf : colnames `BF'
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) fastnose gteffects
if _rc != 0 {
  dis as error "D1 FAIL: fastnose run rc=" _rc
  global pf = $pf + 1
  dis as error "IDENTITY_EXTRAS_FAILS=$pf"
  exit, clear
}
local cnn : colnames e(b)
if "`cnf'" != "`cnn'" {
  dis as error "D2 FAIL: colnames differ: [`cnf'] vs [`cnn']"
  global pf = $pf + 1
}
if mreldif(`BF', e(b)) > 1e-10 {
  dis as error "D3 FAIL: full vs fastnose e(b) mreldif = " %12.4e mreldif(`BF', e(b))
  global pf = $pf + 1
}
else dis "D3 info: full vs fastnose e(b) mreldif = " %12.4e mreldif(`BF', e(b))

* F3: w_att-weighted mean of the gteffects estimates over the post cells == att
tempname ATTFN G
scalar `ATTFN' = _b[main_att]
matrix `G' = e(gteffects)
local gcn : colnames `G'
if "`: word 4 of `gcn''" != "relyr" | "`: word 5 of `gcn''" != "estimate" | "`: word 12 of `gcn''" != "w_att" {
  dis as error "D4 FAIL: unexpected gteffects columns: [`gcn']"
  global pf = $pf + 1
}
mata: G = st_matrix(st_local("G"))
mata: idx = selectindex((G[., 4] :>= 0) :& (G[., 4] :< .) :& (G[., 12] :< .))
mata: st_numscalar("__gtatt", sum(G[idx, 5] :* G[idx, 12]) / sum(G[idx, 12]))
mata: idp = selectindex((G[., 4] :< 0) :& (G[., 12] :< .))
mata: st_numscalar("__gtpre", sum(G[idp, 5] :* G[idp, 12]) / sum(G[idp, 12]))
if reldif(__gtatt, `ATTFN') > 1e-10 {
  dis as error "D5 FAIL: gteffects w_att aggregate " %21.0g __gtatt " vs att " %21.0g `ATTFN' " (reldif " %12.4e reldif(__gtatt, `ATTFN') ")"
  global pf = $pf + 1
}
else dis "D5 info: gteffects att aggregation reldif = " %12.4e reldif(__gtatt, `ATTFN')
* informational only (pre aggregate; uncounted)
dis "D5 info (uncounted): pre-cell w_att aggregate = " %21.0g __gtpre " vs pretreatment_att " %21.0g _b[pretreatment_att] " (reldif " %12.4e reldif(__gtpre, _b[pretreatment_att]) ")"

if $pf == 0 dis "IDENTITY_EXTRAS_ALL_PASS"
else dis as error "IDENTITY_EXTRAS_FAILS=$pf"
exit, clear
