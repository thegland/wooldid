* Targeted suite for the NFADV-01/02/05 adversarial recalibration (wooldid_advstat):
* hand-reconstruction of both p-values from e(b_es_att)/e(V_es_att) + the known
* geometry, bitwise slope invariance vs the parent build, the intercept_rt0
* re-expression, esf scoping (pooled Caution, no p, line retained), the degenerate
* covariance Caution under timetrends, per-subgroup reconstruction, and a null Monte
* Carlo documenting the primary statistic's calibration (and the old frozen-line
* statistic's miscalibration) on a heteroskedastic covariance.
* Run from test/:  "<stata>" -e do test_advstat.do   -> expect ADVSTAT_ALL_PASS
*
* Parent pins: produced from the pre-fix build (git show cb1996c:wooldidinstall/
* wooldid.ado > scratch copy, run separately) on the same two configs. They also
* match the pre-fix battery baselines (test/baseline/esf_adv.csv, adv_altb.csv).
*   esf_adv  (primary): slope +1.5808d87098235X-003  _cons -1.06362b4a54d8fX+000
*   adv_altb (altmeth): slope -1.76046eab851a1X-004  _cons +1.1daaa2f50c7fdX-001
*   min_t = 8 on dta_binary for both configs
clear all
set more off
set maxvar 20000
set graphics off
set seed 4242
set sortseed 4242
local DATA "data"
do "../wooldidinstall/wooldid.ado"

local fails 0

* parent pins as exact decimal literals (17 significant digits round-trip to the
* doubles bit-for-bit; the %21x transcription is re-verified below before use)
scalar advt_pslope1 = 0.16798562139139653
scalar advt_pcons1 = -1.0242640549157362
scalar advt_pslope2 = -0.09131282072514192
scalar advt_pcons2 = 0.5579424785739772
foreach pin in "pslope1 +1.5808d87098235X-003" "pcons1 -1.06362b4a54d8fX+000" "pslope2 -1.76046eab851a1X-004" "pcons2 +1.1daaa2f50c7fdX-001" {
  local pn : word 1 of `pin'
  local px : word 2 of `pin'
  local chk : di %21x scalar(advt_`pn')
  if "`chk'" != "`px'" {
    dis as error "PIN FAIL: advt_`pn' transcribes to `chk', expected `px'"
    local ++fails
  }
}

* independent reconstructions (deliberately a different route than wooldid_advx:
* pinv + GLS normal equations here vs whitened eigenspace least squares there)
mata:
void advt_primary(string scalar bn, string scalar vn, string scalar rn, real scalar dfr)
{
    real colvector b, rel, g, d
    real matrix V, X, Vp
    real rowvector sv
    real scalar rV, q, W, p
    b = st_matrix(bn)'
    V = st_matrix(vn)
    rel = st_matrix(rn)
    X = J(rows(rel), 1, 1), rel
    V = (V + V')/2
    sv = svdsv(V)
    rV = sum(sv :> max(sv)*1e-12)
    st_local("advt_rankV", strofreal(rV))
    Vp = pinv(V)
    g = invsym(X'Vp*X)*(X'Vp*b)
    d = b - X*g
    W = d'Vp*d
    q = rV - 2
    st_local("advt_q", strofreal(q))
    if (q <= 0) {
        st_local("advt_p", "")
        return
    }
    if (dfr < .) p = Ftail(q, dfr, W/q)
    else p = chi2tail(q, W)
    st_numscalar("advt_ps", p)
    st_local("advt_p", strofreal(p, "%21x"))
}

void advt_alt(string scalar bn, string scalar vn, string scalar rn,
              real scalar esrel, real scalar dfr)
{
    real colvector b, rel, idxP, idxQ, bP, relP, d
    real matrix V, ZP, ZQ, Aop, VPP, VPQ, VQQ, S, Sp
    real rowvector sv
    real scalar npr, nq, np, q, W, p
    b = st_matrix(bn)'
    V = st_matrix(vn)
    rel = st_matrix(rn)
    idxP = selectindex(rel :<= -1)
    idxQ = selectindex(rel :>= 0)
    npr = rows(idxP)
    nq = rows(idxQ)
    bP = b[idxP] \ 0
    relP = rel[idxP] \ esrel
    np = rows(relP)
    ZP = J(np, 1, 1), relP
    ZQ = J(nq, 1, 1), rel[idxQ]
    Aop = ZQ*invsym(ZP'ZP)*ZP'
    VPP = J(np, np, 0)
    VPP[|1,1 \ npr,npr|] = V[idxP, idxP]
    VPQ = J(np, nq, 0)
    VPQ[|1,1 \ npr,nq|] = V[idxP, idxQ]
    VQQ = V[idxQ, idxQ]
    d = b[idxQ] - Aop*bP
    S = VQQ + Aop*VPP*Aop' - Aop*VPQ - VPQ'*Aop'
    S = (S + S')/2
    sv = svdsv(S)
    q = sum(sv :> max(sv)*1e-12)
    st_local("advt_q", strofreal(q))
    Sp = pinv(S)
    W = d'Sp*d
    if (dfr < .) p = Ftail(q, dfr, W/q)
    else p = chi2tail(q, W)
    st_numscalar("advt_ps", p)
    st_local("advt_p", strofreal(p, "%21x"))
}
end

* pulls e(b_es_att[_sSG])/e(V_es_att[_sSG]) + relyrs parsed from the pre#/post#
* labels into the (persistent) matrices advt_bb/advt_vv/advt_rr, and df_r into
* the caller local advt_dfr
cap program drop advt_pull
program define advt_pull
  args suffix
  matrix advt_bb = e(b_es_att`suffix')
  matrix advt_vv = e(V_es_att`suffix')
  local lbl : colnames advt_bb
  local k = colsof(advt_bb)
  matrix advt_rr = J(`k', 1, 0)
  local i 0
  foreach l of local lbl {
    local ++i
    if substr("`l'", 1, 3) == "pre" matrix advt_rr[`i', 1] = -real(substr("`l'", 4, .))
    else matrix advt_rr[`i', 1] = real(substr("`l'", 5, .))
  }
  local d .
  if "`e(df_r)'" != "" & "`e(df_r)'" != "." local d `e(df_r)'
  c_local advt_dfr `d'
end

* ---- T1: primary hand-reconstruction on the esf_adv battery config ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod esrelativeto(-2) adversarial
if _rc != 0 | "`e(adversarial_pval_att)'" == "" {
  dis as error "T1 FAIL: rc=" _rc " pval=[`e(adversarial_pval_att)'] (esf primary must produce a p)"
  local ++fails
}
else {
  advt_pull ""
  mata: advt_primary("advt_bb", "advt_vv", "advt_rr", `advt_dfr')
  if `advt_rankV' != colsof(advt_bb) {
    dis as error "T1 FAIL: event-study covariance not full rank (rank `advt_rankV' of `=colsof(advt_bb)') - reconstruction geometry ambiguous"
    local ++fails
  }
  if `advt_q' != `=colsof(advt_bb)' - 2 {
    dis as error "T1 FAIL: residual rank `advt_q' != k_real - 2"
    local ++fails
  }
  if "`advt_p'" == "" {
    dis as error "T1 FAIL: reconstruction refused (q<=0) but wooldid posted a p"
    local ++fails
    exit 9
  }
  if reldif(scalar(advt_ps), e(adversarial_pval_att)) > 1e-10 {
    dis as error "T1 FAIL: reconstructed p `advt_p' vs returned " %21x e(adversarial_pval_att) " (reldif " %8.2e reldif(scalar(advt_ps), e(adversarial_pval_att)) ")"
    local ++fails
  }
  else dis as result "T1: primary p reconstructed at reldif " %8.2e reldif(scalar(advt_ps), e(adversarial_pval_att)) " (q=`advt_q', p=" %12.6e e(adversarial_pval_att) ")"
}
* the parent-build bitwise pins (see header)
tempname L1
matrix `L1' = e(adversarial_line_att)
local news : di %21x `L1'[1,1]
if "`news'" != "+1.5808d87098235X-003" {
  dis as error "T3 FAIL: esf_adv slope `news' != parent +1.5808d87098235X-003 (bitwise)"
  local ++fails
}
tempvar rel1 smp1
qui gen double `smp1' = e(sample)
qui gen double `rel1' = t - ttre if `smp1' & ttre != .
qui sum `rel1', meanonly
local mint = abs(r(min))
if `mint' != 8 {
  dis as error "T3 FAIL: min_t `mint' != 8 on dta_binary"
  local ++fails
}
local wanti : di %21x scalar(advt_pcons1) + (`mint'+1)*scalar(advt_pslope1)
local newi : di %21x `L1'[1,2]
if "`newi'" != "`wanti'" {
  dis as error "T3 FAIL: esf_adv intercept_rt0 `newi' != parent _cons + (min_t+1)*slope = `wanti'"
  local ++fails
}
local cn : colfullnames `L1'
local rn : rowfullnames `L1'
if "`cn'" != "slope intercept_rt0" | "`rn'" != "y1" {
  dis as error "T3 FAIL: line stripe cols=[`cn'] rows=[`rn'] (want slope intercept_rt0 / y1)"
  local ++fails
}

* ---- T2: altmethod hand-reconstruction on the adv_altb battery config ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-3) esprelength(2) espostlength(2) adversarial altmethod_adversarial
if _rc != 0 | "`e(adversarial_pval_att)'" == "" {
  dis as error "T2 FAIL: rc=" _rc " pval=[`e(adversarial_pval_att)'] (esf altmethod must produce a p)"
  local ++fails
}
else {
  advt_pull ""
  mata: advt_alt("advt_bb", "advt_vv", "advt_rr", -3, `advt_dfr')
  if `advt_q' != 2 {
    dis as error "T2 FAIL: forecast-covariance rank `advt_q' != 2"
    local ++fails
  }
  if reldif(scalar(advt_ps), e(adversarial_pval_att)) > 1e-10 {
    dis as error "T2 FAIL: reconstructed p `advt_p' vs returned " %21x e(adversarial_pval_att)
    local ++fails
  }
  else dis as result "T2: altmethod p reconstructed at reldif " %8.2e reldif(scalar(advt_ps), e(adversarial_pval_att)) " (q=`advt_q', p=" %12.6e e(adversarial_pval_att) ")"
  * the reviewer-measured corrected value for this committed cell (F(2,79))
  if reldif(e(adversarial_pval_att), 3.6886e-06) > 1e-3 {
    dis as error "T2 FAIL: p " %12.6e e(adversarial_pval_att) " far from the reviewer-measured 3.6886e-06"
    local ++fails
  }
  * parent-build bitwise pins (see header; min_t + 1 = 9 on dta_binary)
  tempname L2
  matrix `L2' = e(adversarial_line_att)
  local news : di %21x `L2'[1,1]
  local wants : di %21x scalar(advt_pslope2)
  if "`news'" != "`wants'" {
    dis as error "T3 FAIL: adv_altb slope `news' != parent `wants' (bitwise)"
    local ++fails
  }
  local wanti : di %21x scalar(advt_pcons2) + 9*scalar(advt_pslope2)
  local newi : di %21x `L2'[1,2]
  if "`newi'" != "`wanti'" {
    dis as error "T3 FAIL: adv_altb intercept_rt0 `newi' != parent _cons + (min_t+1)*slope = `wanti'"
    local ++fails
  }
}

* ---- T4: pooled reference regime -> Caution, no p, line retained ----
use "`DATA'/dta_binary", clear
qui log using advstat_t4.txt, text replace name(advst)
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) adversarial
local rc4 = _rc
qui log close advst
if `rc4' != 0 {
  dis as error "T4 FAIL: rc=" `rc4'
  local ++fails
}
else {
  if "`e(adversarial_pval_att)'" != "" {
    dis as error "T4 FAIL: adversarial_pval_att posted under the pooled reference regime"
    local ++fails
  }
  cap confirm matrix e(adversarial_line_att)
  if _rc != 0 {
    dis as error "T4 FAIL: the display line must survive the pooled-regime scoping"
    local ++fails
  }
  if strpos(fileread("advstat_t4.txt"), "Caution: the adversarial pre-trend p-value requires esfixedbaseperiod") == 0 {
    dis as error "T4 FAIL: no pooled-regime Caution in the console output"
    local ++fails
  }
}

* ---- T5: timetrends + esf -> degenerate-covariance Caution, no p, line retained ----
* (previously the frozen-line `test' silently dropped all constraints here and the
* posted p was missing-valued)
use "`DATA'/dta_binary", clear
qui log using advstat_t5.txt, text replace name(advst)
cap noi wooldid y i t ttre, att cluster(i) timetrends esfixedbaseperiod esrelativeto(-1) esprelength(3) espostlength(3) adversarial
local rc5 = _rc
qui log close advst
if `rc5' != 0 {
  dis as error "T5 FAIL: rc=" `rc5'
  local ++fails
}
else {
  if "`e(adversarial_pval_att)'" != "" {
    dis as error "T5 FAIL: adversarial_pval_att posted on a degenerate (timetrends) covariance"
    local ++fails
  }
  cap confirm matrix e(adversarial_line_att)
  if _rc != 0 {
    dis as error "T5 FAIL: the display line must survive the degenerate-covariance refusal"
    local ++fails
  }
  if strpos(fileread("advstat_t5.txt"), "the composite-null covariance is degenerate") == 0 {
    dis as error "T5 FAIL: no degenerate-covariance Caution in the console output"
    local ++fails
  }
}

* ---- T6: per-subgroup esf p-values finite and reconstructable the same way ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) esfixedbaseperiod esrelativeto(-2) esprelength(2) espostlength(2) adversarial
if _rc != 0 {
  dis as error "T6 FAIL: rc=" _rc
  local ++fails
}
else {
  local sgn 0
  local sc : e(scalars)
  foreach s of local sc {
    if strpos("`s'", "adversarial_pval_att_s") == 1 {
      local ++sgn
      * "adversarial_pval_att" is 20 characters; the "_s<sg>" suffix starts at 21
      local sgsfx = substr("`s'", 21, .)
      if missing(e(`s')) {
        dis as error "T6 FAIL: e(`s') is missing"
        local ++fails
        continue
      }
      advt_pull "`sgsfx'"
      mata: advt_primary("advt_bb", "advt_vv", "advt_rr", `advt_dfr')
      if "`advt_p'" == "" {
        dis as error "T6 FAIL: reconstruction refused for `s' but wooldid posted a p"
        local ++fails
      }
      else if reldif(scalar(advt_ps), e(`s')) > 1e-10 {
        dis as error "T6 FAIL: `s' reconstructed `advt_p' vs returned " %21x e(`s')
        local ++fails
      }
      else dis as result "T6: `s' reconstructed at reldif " %8.2e reldif(scalar(advt_ps), e(`s'))
    }
  }
  if `sgn' == 0 {
    dis as error "T6 FAIL: no per-subgroup adversarial p-values posted"
    local ++fails
  }
}

* ---- T7: null Monte Carlo (outside wooldid, pure mata, known covariance) ----
* geometry: real horizons (-3,-2,0) with planted reference (-1,0); heteroskedastic
* V = diag(25, .01, .01); true effects all 0 (an exact line through the reference).
* 200 draws at nominal 5%: the recalibrated primary statistic (q = 1) must reject at
* ~5% ([0.02, 0.09] Monte-Carlo band); the OLD frozen-line statistic (fitted values
* frozen as constants, full V^-1 weighting, q = 3) must reject far outside the band
* (reviewer-measured ~0.81 on this covariance) - documenting the cure.
mata:
rseed(20260819)
rel = (-3 \ -2 \ 0)
relf = rel \ -1                       // + planted reference row
sd = (5 \ .1 \ .1)
V = diag(sd:^2)
X = J(3,1,1), rel
Xf = J(4,1,1), relf
Vp = invsym(V)                        // full rank here; GLS weight for the new stat
rejnew = 0
rejold = 0
reps = 200
for (s = 1; s <= reps; s++) {
    b = sd :* rnormal(3, 1, 0, 1)
    // new: minimum-distance under V^{-1}, q = rank(V) - 2 = 1
    g = invsym(X'Vp*X)*(X'Vp*b)
    d = b - X*g
    Wn = d'Vp*d
    if (chi2tail(1, Wn) < .05) rejnew++
    // old: unweighted OLS on real + planted rows, fitted values frozen, full V^{-1}
    bf = b \ 0
    go = invsym(Xf'Xf)*(Xf'bf)
    do_ = b - X*go
    Wo = do_'Vp*do_
    if (chi2tail(3, Wo) < .05) rejold++
}
st_local("mcnew", strofreal(rejnew/reps))
st_local("mcold", strofreal(rejold/reps))
end
dis as result "T7: null MC rejection at 5% - recalibrated primary = `mcnew', old frozen-line = `mcold' (200 sims)"
if `mcnew' < 0.02 | `mcnew' > 0.09 {
  dis as error "T7 FAIL: recalibrated primary rejection `mcnew' outside [0.02, 0.09]"
  local ++fails
}
if `mcold' >= 0.02 & `mcold' <= 0.09 {
  dis as error "T7 FAIL: old frozen-line rejection `mcold' inside the nominal band - no cure documented"
  local ++fails
}

if `fails' == 0 dis "ADVSTAT_ALL_PASS"
else dis as error "ADVSTAT_FAILS=`fails'"
