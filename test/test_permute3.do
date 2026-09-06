* Targeted suite for the 2026-09-05 permute() completeness build (chunk 1): every
* reported p-value and interval under permute() is a permutation quantity, coverage by
* subtraction (blank, never analytic), the shift intervals, the fast-mode default and its
* opt-out, the fixed-metric joint/equality/adversarial tests, contrasts, cells, plots,
* permsave(), level(), refusals, replay. Every kernel number is reconstructed in Mata
* from the posted draws (never through a local: a decimal round-trip loses bits).
* Run from test/:  "<stata>" -e do test_permute3.do   -> expect PERMUTE3_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 24680
set sortseed 24680
local DATA "data"
do "../wooldidinstall/wooldid.ado"
global pf 0

cap program drop xrc
program define xrc
  args label want got
  if `got' != `want' {
    dis as error "`label' FAIL: rc=`got' (expected `want')"
    global pf = $pf + 1
  }
end

* fixture: nclus clusters (units), the first ntreat treated at t = 4, 8 periods, a
* cluster-level subgroup and a unit-level dose; effect eff, larger in subgroup 1
cap program drop mkperm3
program define mkperm3
  args nclus ntreat eff
  clear
  set seed 4242
  set sortseed 4242
  qui set obs `nclus'
  qui gen i = _n
  qui gen ttre = cond(_n <= `ntreat', 4, .)
  qui gen sg = mod(i, 2)
  qui gen dose = 0.5 + runiform()
  qui expand 8
  qui bysort i: gen t = _n
  qui gen y = `eff'*(t >= ttre & ttre < .)*(1 + 0.5*sg) + 0.1*i + rnormal()
end

* Mata: recompute every element's p and interval from the posted draws and compare
* BITWISE to e(perm_p)/e(perm_ci); returns the number of mismatches
mata:
real scalar wdp3_q(real colvector v, real scalar pr)
{
  real scalar n, x
  n = rows(v)
  x = n*pr
  if (x == floor(x) & x >= 1) return((v[x] + v[min((x+1, n))])/2)
  return(v[max((min((ceil(x), n)), 1))])
}
real scalar wdp3_recon(real scalar exact, real scalar lvl) return(wdp3_reconx(exact, lvl, 0))
// absp = 1: the permabsp convention - p from |s|, the symmetric dual interval; the null
// band e(perm_nullband) must hold the SIGNED level quantiles in both conventions
real scalar wdp3_reconx(real scalar exact, real scalar lvl, real scalar absp)
{
  real matrix A, S, CI, NB
  real rowvector cm, cs, P, bb, N
  real colvector v, s
  real scalar k, r, n, plo, phi, p, qlo, qhi, qa, bad, so, sobs, ts, M, al
  A = st_matrix("e(permdist_b)")
  cm = st_matrix("e(perm_comparator_b)")
  P = st_matrix("e(perm_p)")
  CI = st_matrix("e(perm_ci)")
  NB = st_matrix("e(perm_nullband)")
  bb = st_matrix("e(b)")
  N = st_matrix("e(perm_n)")
  M = st_numscalar("e(perm_reps)")
  al = (100 - lvl)/200
  ts = (st_global("e(perm_stat)") == "t")
  if (ts) {
    S = st_matrix("e(permdist_se)")
    cs = st_matrix("e(perm_comparator_se)")
  }
  bad = 0
  for (k = 1; k <= cols(A); k++) {
    if (!ts) {
      s = A[., k]
      sobs = cm[k]
      so = 1
    }
    else {
      s = J(rows(A), 1, .)
      for (r = 1; r <= rows(A); r++) {
        if (A[r, k] < . & S[r, k] < . & S[r, k] > 0) s[r] = A[r, k]/S[r, k]
      }
      so = cs[k]
      sobs = (so < . & so > 0 ? cm[k]/so : .)
    }
    v = select(s, s :< .)
    n = rows(v)
    if (n != N[k]) {
      bad++
      continue
    }
    if (n < 2 | sobs >= .) {
      if (P[k] < .) bad++
      continue
    }
    v = sort(v, 1)
    if (absp) {
      if (exact & n == M) p = sum(abs(v) :>= abs(sobs))/n
      else p = (1 + sum(abs(v) :>= abs(sobs)))/(n + 1)
    }
    else {
      if (exact & n == M) {
        plo = sum(v :<= sobs)/n
        phi = sum(v :>= sobs)/n
      }
      else {
        plo = (1 + sum(v :<= sobs))/(n + 1)
        phi = (1 + sum(v :>= sobs))/(n + 1)
      }
      p = min((1, 2*min((plo, phi))))
    }
    if (p != P[k]) bad++
    qlo = wdp3_q(v, al)
    qhi = wdp3_q(v, 1 - al)
    if (NB[1, k] != qlo | NB[2, k] != qhi) bad++
    if (qlo > 0 | qhi < 0) {
      if (CI[1, k] < . | CI[2, k] < .) bad++
      continue
    }
    if (absp) {
      qa = wdp3_q(sort(abs(v), 1), 1 - 2*al)
      if (bb[k] - qa*so != CI[1, k] | bb[k] + qa*so != CI[2, k]) bad++
    }
    else if (bb[k] - qhi*so != CI[1, k] | bb[k] - qlo*so != CI[2, k]) bad++
  }
  return(bad)
}
// upper-tail p of an observed statistic among replication statistics (add-one)
real scalar wdp3_ut(real colvector W, real scalar Wobs, real scalar exact, real scalar M)
{
  real colvector v
  real scalar n
  v = select(W, W :< .)
  n = rows(v)
  if (exact & n == M) return(sum(v :>= Wobs)/n)
  return((1 + sum(v :>= Wobs))/(n + 1))
}
end

* the table rows <-> e(perm_p) columns are positional (main rows first, then the
* subgroup tables in e(wdid_sglist) order): check that every displayed p/interval
* cell equals the posted vectors bitwise and that the se/t columns are as expected
cap program drop chk3tab
program define chk3tab, rclass
  syntax , label(string) [ase ts]
  local bad 0
  tempname W
  matrix `W' = e(wdidmainresults)
  local cn : colnames `W'
  if "`: word 2 of `cn''" != "se_analytic" | "`: word 3 of `cn''" != "t_perm" | "`: word 4 of `cn''" != "p_perm" local ++bad
  local kk 0
  forval r = 1/`=rowsof(`W')' {
    local ++kk
    if "`ase'`ts'" == "" & `W'[`r', 2] < . local ++bad
    if "`ts'" == "" & `W'[`r', 3] < . local ++bad
    if "`ts'" != "" & `W'[`r', 4] < . & `W'[`r', 3] >= . local ++bad
    if `W'[`r', 4] != e(perm_p)[1, `kk'] local ++bad
    if `W'[`r', 5] != e(perm_ci)[1, `kk'] | `W'[`r', 6] != e(perm_ci)[2, `kk'] local ++bad
    if `W'[`r', 1] < . & `W'[`r', 4] >= . local ++bad
    if `W'[`r', 4] < . & e(perm_kernel)[10, `kk'] != 1 & (`W'[`r', 5] >= . | `W'[`r', 6] >= .) local ++bad
  }
  foreach sg in `e(wdid_sglist)' {
    cap confirm matrix e(wdidsgresults_`sg')
    if _rc != 0 continue
    tempname Ws
    matrix `Ws' = e(wdidsgresults_`sg')
    forval r = 1/`=rowsof(`Ws')' {
      local ++kk
      if `Ws'[`r', 4] != e(perm_p)[1, `kk'] local ++bad
      if `Ws'[`r', 5] != e(perm_ci)[1, `kk'] | `Ws'[`r', 6] != e(perm_ci)[2, `kk'] local ++bad
    }
  }
  if `kk' != colsof(e(perm_p)) local ++bad
  if `bad' != 0 {
    dis as error "`label' FAIL: `bad' table-vs-kernel check(s)"
    global pf = $pf + 1
  }
  return scalar bad = `bad'
end

* ================= Q1: sampled shuffle, fast-mode default: provenance, blanks, header,
* footer wording, plots, permsave ====================================================
mkperm3 8 4 3
set seed 101
qui log using perm3_q1.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) makeplots saveplots(perm3_q1) permsave(perm3_q1) replace
local rc = _rc
qui log close p3
xrc Q1a 0 `rc'
if `rc' == 0 {
  mata: st_numscalar("__bad", wdp3_recon(0, 95))
  if __bad != 0 {
    dis as error "Q1b FAIL: " __bad " draw-reconstruction mismatch(es)"
    global pf = $pf + 1
  }
  chk3tab , label(Q1c)
  local bad 0
  cap confirm matrix e(V)
  if _rc == 0 local ++bad
  if "`e(wdid_fastnose)'" != "" local ++bad
  if strpos("`e(wdid_vcespec)'", "randomization inference") == 0 local ++bad
  if e(perm_completed) != 30 | e(perm_reps) != 30 | e(perm_level) != 95 | e(perm_exact) != 0 local ++bad
  if "`e(perm_stat)'" != "b" local ++bad
  * the legacy scalars are the SAME numbers the tables carry
  if e(perm_p_att) != e(wdidmainresults)[1, 4] local ++bad
  if e(perm_lb_att) != e(perm_nullband)[1, 1] | e(perm_ub_att) != e(perm_nullband)[2, 1] local ++bad
  if e(perm_sd_att) != e(perm_sd)[1, 1] local ++bad
  * the comparator equals the fast-mode observed estimate bitwise
  if mreldif(e(perm_comparator_b), e(b)) != 0 local ++bad
  if `bad' {
    dis as error "Q1d FAIL: `bad' e()/header check(s)"
    global pf = $pf + 1
  }
  * log scan: no analytic wording, the footer present
  local lg = fileread("perm3_q1.txt")
  if strpos(`"`lg'"', "fast mode") | strpos(`"`lg'"', "on a single line") | strpos(`"`lg'"', "Standard Errors: cluster(") {
    dis as error "Q1e FAIL: analytic wording in the permute log"
    global pf = $pf + 1
  }
  if strpos(`"`lg'"', "Randomization Inference") == 0 | strpos(`"`lg'"', "shift interval") == 0 | strpos(`"`lg'"', "Permutation replications (30)") == 0 {
    dis as error "Q1f FAIL: footer/progress wording missing"
    global pf = $pf + 1
  }
  * the ES plot data carries the table's shift intervals bitwise
  tempname W
  matrix `W' = e(wdidmainresults)
  preserve
  use perm3_q1_ES.dta, clear
  qui drop if rownamesvar == "reference"
  local nb 0
  local nl 0
  forval r = 1/`=_N' {
    local rl = relyr[`r']
    forval j = 1/`=rowsof(`W')' {
      if `W'[`j', 7] == `rl' {
        local ++nl
        if lb_95ci[`r'] != `W'[`j', 5] | ub_95ci[`r'] != `W'[`j', 6] | estimate[`r'] != `W'[`j', 1] local ++nb
      }
    }
  }
  restore
  if `nl' != 4 | `nb' != 0 {
    dis as error "Q1g FAIL: ES plot data vs table: `nl' rows matched, `nb' bad"
    global pf = $pf + 1
  }
  cap confirm file perm3_q1_ES_att.gph
  if _rc != 0 {
    dis as error "Q1h FAIL: ES plot file not written"
    global pf = $pf + 1
  }
  * permsave(): rows, the observed row, the onset block, and a p recomputed from the file
  tempname PB CB PV
  matrix `PB' = e(b)
  matrix `CB' = e(perm_comparator_b)
  matrix `PV' = e(perm_p)
  local nrep = e(perm_reps)
  preserve
  use perm3_q1.dta, clear
  local bad 0
  if _N != `nrep' + 1 local ++bad
  if _rep[1] != 0 local ++bad
  if main_att[1] != `PB'[1, 1] local ++bad
  cap confirm variable onset_1 onset_8
  if _rc != 0 local ++bad
  qui count if _rep == 0 & onset_1 == 4 & onset_8 >= .
  if r(N) != 1 local ++bad
  mata: v = st_data(., "main_att"); v = select(v[2::rows(v)], v[2::rows(v)] :< .); n = rows(v); c0 = st_matrix(st_local("CB"))[1, 1]; p = min((1, 2*min(((1 + sum(v :<= c0))/(n + 1), (1 + sum(v :>= c0))/(n + 1))))); st_numscalar("__pf", p)
  if __pf != `PV'[1, 1] local ++bad
  restore
  if `bad' {
    dis as error "Q1i FAIL: `bad' permsave check(s)"
    global pf = $pf + 1
  }
  * replay prints the footer and the rewired table
  qui log using perm3_q1r.txt, text replace name(p3)
  cap noi wooldid
  local rcr = _rc
  qui log close p3
  if `rcr' != 0 | strpos(fileread("perm3_q1r.txt"), "Randomization Inference") == 0 | strpos(fileread("perm3_q1r.txt"), "p_perm") == 0 {
    dis as error "Q1j FAIL: replay rc=`rcr' or footer/table missing"
    global pf = $pf + 1
  }
}

* ================= Q2: leak check against the analytic twin; permanalyticse ==========
mkperm3 8 4 3
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
tempname AN
matrix `AN' = e(wdidmainresults)
mkperm3 8 4 3
set seed 102
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30)
xrc Q2a 0 _rc
if _rc == 0 {
  local bad 0
  forval r = 1/`=rowsof(`AN')' {
    * a permutation p equal to the analytic p is a leak; a blank se is required
    if e(wdidmainresults)[`r', 4] == `AN'[`r', 4] local ++bad
    if e(wdidmainresults)[`r', 2] < . local ++bad
    if e(wdidmainresults)[`r', 5] == `AN'[`r', 5] | e(wdidmainresults)[`r', 6] == `AN'[`r', 6] local ++bad
  }
  if `bad' {
    dis as error "Q2b FAIL: `bad' analytic cell(s) leaked into the permute table"
    global pf = $pf + 1
  }
}
mkperm3 8 4 3
set seed 102
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) permanalyticse
xrc Q2c 0 _rc
if _rc == 0 {
  local bad 0
  forval r = 1/`=rowsof(`AN')' {
    if e(wdidmainresults)[`r', 2] != `AN'[`r', 2] local ++bad
    if e(wdidmainresults)[`r', 3] < . local ++bad
    if e(wdidmainresults)[`r', 4] == `AN'[`r', 4] local ++bad
  }
  cap confirm matrix e(V)
  if _rc != 0 local ++bad
  if "`e(perm_analyticse)'" != "1" local ++bad
  mata: st_numscalar("__bad", wdp3_recon(0, 95))
  if __bad != 0 local ++bad
  chk3tab , label(Q2d) ase
  if `bad' {
    dis as error "Q2e FAIL: `bad' permanalyticse check(s)"
    global pf = $pf + 1
  }
}

* ================= Q3: two seeds move every non-saturated p and bound ================
mkperm3 8 4 3
set seed 103
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30)
tempname S1
matrix `S1' = e(wdidmainresults)
mkperm3 8 4 3
set seed 104
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30)
local moved 0
local same 0
forval r = 1/`=rowsof(`S1')' {
  if e(wdidmainresults)[`r', 5] != `S1'[`r', 5] | e(wdidmainresults)[`r', 6] != `S1'[`r', 6] local ++moved
  else local ++same
}
* every interval bound must move with the seed (p can saturate at the add-one floor)
if `same' != 0 {
  dis as error "Q3 FAIL: `same' row(s) had identical intervals under two seeds"
  global pf = $pf + 1
}

* ================= Q4: exact enumeration, mirror symmetry: p = 2/20 and the shift lower
* bound is exactly 0 (the comparator IS the maximum draw) ==============================
mkperm3 6 3 3
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
xrc Q4a 0 _rc
if _rc == 0 {
  local bad 0
  if abs(e(perm_p_att) - 0.10) > 1e-12 local ++bad
  if e(wdidmainresults)[1, 4] != e(perm_p_att) local ++bad
  if e(perm_exact) != 1 | e(perm_completed) != 20 local ++bad
  if e(wdidmainresults)[1, 5] != 0 local ++bad
  if e(wdidmainresults)[1, 6] <= _b[main_att] local ++bad
  mata: st_numscalar("__bad", wdp3_recon(1, 95))
  if __bad != 0 local ++bad
  if `bad' {
    dis as error "Q4b FAIL: `bad' exact-mode check(s): p=" e(perm_p_att) " lb=" e(wdidmainresults)[1, 5]
    global pf = $pf + 1
  }
}

* ================= Q5: subgroups + jointtests + adversarial: the fixed-metric tests are
* reconstructable from e(perm_*W)/e(perm_*Wobs); contrasts from the draws =============
mkperm3 8 4 3
set seed 105
qui log using perm3_q5.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) subgroups(sg) jointtests adversarial esfixedbaseperiod esrelativeto(-1)
local rc = _rc
qui log close p3
xrc Q5a 0 `rc'
if `rc' == 0 {
  local bad 0
  chk3tab , label(Q5b)
  mata: st_numscalar("__bad", wdp3_recon(0, 95))
  if __bad != 0 local ++bad
  * joint pre-trend test: overall and per subgroup
  foreach en in joint_horiz joint_horiz_s0 joint_horiz_s1 {
    cap local jc = colnumb(e(perm_jointW), "`en'")
    if _rc != 0 | `jc' >= . {
      local ++bad
      continue
    }
    mata: st_numscalar("__pj", wdp3_ut(st_matrix("e(perm_jointW)")[., `jc'], st_matrix("e(perm_jointWobs)")[1, `jc'], 0, .))
    if "`en'" == "joint_horiz" {
      if e(joint_horiz) != __pj local ++bad
    }
    else if reldif(real("`e(`en')'"), __pj) > 1e-13 local ++bad
    * (the per-subgroup joint p's are ereturn LOCALS by the analytic convention, so they
    * carry a 16-digit decimal round-trip; the scalar ones are compared bitwise)
  }
  * across-subgroup equality
  cap local sc = colnumb(e(perm_sgW), "joint_att_allsgssame")
  if _rc != 0 | `sc' >= . local ++bad
  else {
    mata: st_numscalar("__ps", wdp3_ut(st_matrix("e(perm_sgW)")[., `sc'], st_matrix("e(perm_sgWobs)")[1, `sc'], 0, .))
    if e(joint_att_allsgssame) != __ps local ++bad
  }
  * adversarial (esf): overall and per subgroup, sharp-null wording in the log
  foreach en in adversarial_pval_att adversarial_pval_att_s0 adversarial_pval_att_s1 {
    cap local ac = colnumb(e(perm_advW), "`en'")
    if _rc != 0 | `ac' >= . {
      local ++bad
      continue
    }
    mata: st_numscalar("__pa", wdp3_ut(st_matrix("e(perm_advW)")[., `ac'], st_matrix("e(perm_advWobs)")[1, `ac'], 0, .))
    if e(`en') != __pa local ++bad
  }
  cap confirm matrix e(adversarial_line_att)
  if _rc != 0 local ++bad
  local lg = fileread("perm3_q5.txt")
  if strpos(`"`lg'"', "on a single line") | strpos(`"`lg'"', "sharp null of no effect at any horizon") == 0 local ++bad
  if strpos(`"`lg'"', "joint test that all pre-treatment att horizon estimates are 0") == 0 local ++bad
  * pairwise contrast row att:1vs0.sg from the draws
  cap confirm matrix e(pwsgs)
  if _rc != 0 local ++bad
  else {
    local prn : rowfullnames e(pwsgs)
    local rr : list posof "att:1vs0.sg" in prn
    if `rr' == 0 local ++bad
    else {
      local ka = colnumb(e(permdist_b), "main_SG1_att")
      local kb = colnumb(e(permdist_b), "main_SG0_att")
      mata: A = st_matrix("e(permdist_b)")[., (`ka', `kb')]; ok = rowmissing(A) :== 0; d = select(A[., 1] :- A[., 2], ok); n = rows(d); cd = st_matrix("e(perm_comparator_b)"); c0 = cd[`ka'] - cd[`kb']; p = min((1, 2*min(((1 + sum(d :<= c0))/(n + 1), (1 + sum(d :>= c0))/(n + 1))))); st_numscalar("__pc", p); v = sort(d, 1); st_numscalar("__lo", st_matrix("e(pwsgs)")[`rr', 1] - wdp3_q(v, 0.975)); st_numscalar("__hi", st_matrix("e(pwsgs)")[`rr', 1] - wdp3_q(v, 0.025))
      mata: st_numscalar("__off", (wdp3_q(v, 0.025) > 0 | wdp3_q(v, 0.975) < 0))
      if e(pwsgs)[`rr', 4] != __pc local ++bad
      if __off == 0 & (e(pwsgs)[`rr', 5] != __lo | e(pwsgs)[`rr', 6] != __hi) local ++bad
      if __off == 1 & (e(pwsgs)[`rr', 5] < . | e(pwsgs)[`rr', 6] < .) local ++bad
      if e(pwsgs)[`rr', 2] < . | e(pwsgs)[`rr', 3] < . local ++bad
      local pwc : colnames e(pwsgs)
      if "`: word 4 of `pwc''" != "p_perm" local ++bad
    }
  }
  if "`e(perm_precaveat)'" != "1" local ++bad
  if `bad' {
    dis as error "Q5c FAIL: `bad' joint/equality/adversarial/contrast check(s)"
    global pf = $pf + 1
  }
}

* ================= Q6: permtstat: se_analytic is the analytic se, t_perm = b/se, the
* t-kernel reconstructs from e(permdist_b)/e(permdist_se) =============================
mkperm3 8 4 3
set seed 106
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) permtstat
xrc Q6a 0 _rc
if _rc == 0 {
  local bad 0
  if "`e(perm_stat)'" != "t" local ++bad
  forval r = 1/`=rowsof(`AN')' {
    if e(wdidmainresults)[`r', 2] != `AN'[`r', 2] local ++bad
    if e(wdidmainresults)[`r', 4] == `AN'[`r', 4] local ++bad
    if e(wdidmainresults)[`r', 3] < . & e(wdidmainresults)[`r', 3] != e(perm_comparator_b)[1, `r']/e(perm_comparator_se)[1, `r'] local ++bad
  }
  mata: st_numscalar("__bad", wdp3_recon(0, 95))
  if __bad != 0 local ++bad
  chk3tab , label(Q6b) ts
  cap confirm matrix e(permdist_se)
  if _rc != 0 local ++bad
  if `bad' {
    dis as error "Q6c FAIL: `bad' permtstat check(s)"
    global pf = $pf + 1
  }
}

* ================= Q7: refusals and screens ==========================================
mkperm3 8 4 3
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) esfixedbaseperiod pretrendpower(80)
xrc Q7a 198 _rc
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) subgroups(sg) suppresssgprimaryeffects
xrc Q7b 198 _rc
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) level(5)
xrc Q7c 198 _rc
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permsave(x)
xrc Q7d 198 _rc
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permanalyticse
xrc Q7e 198 _rc
* a negative subgroup code on a NEVER-TREATED unit passes the observed run's screen and
* would kill the replications that treat it: refused before the loop
qui replace sg = -1 if i == 8
qui log using perm3_q7.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) subgroups(sg)
local rc = _rc
qui log close p3
xrc Q7f 411 `rc'
if strpos(fileread("perm3_q7.txt"), "Permutation replications") != 0 {
  dis as error "Q7g FAIL: the subgroup screen fired after the replication loop started"
  global pf = $pf + 1
}
mkperm3 8 4 3
* permsave() without replace onto an existing file: refused at parse (before the loop)
qui log using perm3_q7b.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) permsave(perm3_q1)
local rc = _rc
qui log close p3
xrc Q7h 602 `rc'
if strpos(fileread("perm3_q7b.txt"), "Permutation replications") != 0 {
  dis as error "Q7i FAIL: permsave existence check fired after the replication loop started"
  global pf = $pf + 1
}
* a user-typed fastnose with permtstat stays refused; jointtests under permute is exempt
* from the analytic clustered-SE screen (cluster at i)
cap wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) permtstat fastnose
xrc Q7j 198 _rc

* ================= Q8: gteffects blanked; level(90) column names, plots at level 90 ====
mkperm3 8 4 3
set seed 108
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) gteffects level(90) makeplots saveplots(perm3_q8) replace
xrc Q8a 0 _rc
if _rc == 0 {
  local bad 0
  if "`e(perm_gt_blank)'" != "1" local ++bad
  tempname G
  matrix `G' = e(gteffects)
  forval r = 1/`=rowsof(`G')' {
    forval cc = 6/9 {
      if `G'[`r', `cc'] < . local ++bad
    }
  }
  local cn : colnames e(wdidmainresults)
  if "`: word 5 of `cn''" != "lb_90ci" | "`: word 6 of `cn''" != "ub_90ci" local ++bad
  if e(perm_level) != 90 local ++bad
  mata: st_numscalar("__bad", wdp3_recon(0, 90))
  if __bad != 0 local ++bad
  cap confirm file perm3_q8_ES_att.gph
  if _rc != 0 local ++bad
  tempname W
  matrix `W' = e(wdidmainresults)
  preserve
  use perm3_q8_ES.dta, clear
  qui drop if rownamesvar == "reference"
  local nb 0
  forval r = 1/`=_N' {
    local rl = relyr[`r']
    forval j = 1/`=rowsof(`W')' {
      if `W'[`j', 7] == `rl' {
        if lb_95ci[`r'] != `W'[`j', 5] | ub_95ci[`r'] != `W'[`j', 6] local ++nb
      }
    }
  }
  restore
  if `nb' != 0 local ++bad
  if `bad' {
    dis as error "Q8b FAIL: `bad' gteffects/level(90) check(s)"
    global pf = $pf + 1
  }
}

* ================= Q9: continuous treatment: cfx cells permuted for latticetype all,
* blanked for baseTreated; plots redrawn from the rewired matrices =====================
mkperm3 8 4 3
set seed 109
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) contreat(dose) contreatcontrols(basic) lattice(all 2) makecfxplots permute(shuffle) permreps(30) saveplots(perm3_q9) replace
xrc Q9a 0 _rc
if _rc == 0 {
  local bad 0
  local nm : word 1 of `e(wdid_cfxmats)'
  cap confirm matrix e(`nm')
  if _rc != 0 local ++bad
  else {
    mata: st_numscalar("__bad", wdp3_recon(0, 95))
    if __bad != 0 local ++bad
    chk3tab , label(Q9b)
    tempname CF
    matrix `CF' = e(`nm')
    cap confirm matrix e(permdist_`nm')
    if _rc != 0 local ++bad
    else {
      forval r = 1/`=rowsof(`CF')' {
        mata: v = st_matrix("e(permdist_`nm')")[., `r']; v = select(v, v :< .); n = rows(v); c0 = st_matrix("e(permcmp_`nm')")[1, `r']; p = (n >= 2 & c0 < . ? min((1, 2*min(((1 + sum(v :<= c0))/(n + 1), (1 + sum(v :>= c0))/(n + 1))))) : .); st_numscalar("__pc", p); v = sort(v, 1); st_numscalar("__lo", (n >= 2 ? st_matrix("`CF'")[`r', 1] - wdp3_q(v, 0.975) : .)); st_numscalar("__hi", (n >= 2 ? st_matrix("`CF'")[`r', 1] - wdp3_q(v, 0.025) : .))
        mata: st_numscalar("__off", (n >= 2 ? (wdp3_q(v, 0.025) > 0 | wdp3_q(v, 0.975) < 0) : 0))
        if `CF'[`r', 4] != __pc local ++bad
        if __pc < . & __off == 0 & (`CF'[`r', 5] != __lo | `CF'[`r', 6] != __hi) local ++bad
        if __off == 1 & (`CF'[`r', 5] < . | `CF'[`r', 6] < .) local ++bad
        if `CF'[`r', 2] < . | `CF'[`r', 3] < . local ++bad
      }
    }
    if "`e(perm_cfx_blank)'" != "" local ++bad
    * the redrawn plot data equals the rewired matrix
    preserve
    cap use perm3_q9_cfx.dta, clear
    if _rc != 0 local ++bad
    else {
      qui keep if matname == "`nm'"
      local nb 0
      forval r = 1/`=_N' {
        if b[`r'] != `CF'[`r', 1] | p[`r'] != `CF'[`r', 4] local ++nb
        if ll[`r'] != `CF'[`r', 5] | ul[`r'] != `CF'[`r', 6] local ++nb
      }
      if `nb' != 0 | _N != rowsof(`CF') local ++bad
    }
    restore
  }
  if `bad' {
    dis as error "Q9c FAIL: `bad' cfx cell check(s)"
    global pf = $pf + 1
  }
}
mkperm3 8 4 3
set seed 110
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) contreat(dose) contreatcontrols(basic) lattice(baseTreated 2) makecfxplots permute(shuffle) permreps(30)
xrc Q9d 0 _rc
if _rc == 0 {
  local bad 0
  if "`e(perm_cfx_blank)'" != "1" local ++bad
  local nm : word 1 of `e(wdid_cfxmats)'
  tempname CF
  matrix `CF' = e(`nm')
  forval r = 1/`=rowsof(`CF')' {
    forval cc = 4/6 {
      if `CF'[`r', `cc'] < . local ++bad
    }
  }
  if `bad' {
    dis as error "Q9e FAIL: `bad' drifting-lattice blank check(s)"
    global pf = $pf + 1
  }
}

* ================= Q10: permute(probabilities) and permstrata run through the new engine
mkperm3 8 4 3
set seed 111
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(30) permstrata(sg) subgroups(sg)
xrc Q10a 0 _rc
if _rc == 0 {
  mata: st_numscalar("__bad", wdp3_recon(0, 95))
  if __bad != 0 {
    dis as error "Q10b FAIL: " __bad " reconstruction mismatch(es) under permstrata"
    global pf = $pf + 1
  }
}
mkperm3 8 4 3
set seed 112
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(probabilities) permreps(30) permprobs((ttre 4 .) (1 .5 .5) (2 .5 .5) (3 .5 .5) (4 .5 .5) (5 .5 .5) (6 .5 .5) (7 .5 .5) (8 .5 .5)) dots(0)
xrc Q10c 0 _rc
if _rc == 0 {
  mata: st_numscalar("__bad", wdp3_recon(0, 95))
  if __bad != 0 {
    dis as error "Q10d FAIL: " __bad " reconstruction mismatch(es) under permute(probabilities)"
    global pf = $pf + 1
  }
}

* ================= Q12: permabsp - the absolute-value p convention with its symmetric dual
* interval, reconstructed bitwise; it must DIFFER from the equal-tailed p on the same draws
* (non-vacuity), coincide on the mirror-symmetric exact fixture (2/20), coincide with the
* permtci interval under permtstat, and be refused without permute() ==================
mkperm3 8 4 3
set seed 1212
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) permute(shuffle) permreps(30)
xrc Q12a 0 _rc
tempname ETP ETD ETB
if _rc == 0 {
  matrix `ETP' = e(perm_p)
  matrix `ETD' = e(permdist_b)
  matrix `ETB' = e(perm_nullband)
  if "`e(perm_pconv)'" != "equaltailed" | e(perm_absp) != 0 {
    dis as error "Q12b FAIL: e(perm_pconv)/e(perm_absp) wrong on the equal-tailed run"
    global pf = $pf + 1
  }
}
mkperm3 8 4 3
set seed 1212
qui log using perm3_q12.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) permute(shuffle) permreps(30) permabsp
local rc12 = _rc
qui log close p3
xrc Q12c 0 `rc12'
if `rc12' == 0 {
  local bad 0
  if "`e(perm_pconv)'" != "absolute" | e(perm_absp) != 1 local ++bad
  * same seed, same draws; the band (signed level quantiles) is convention-free
  if mreldif(e(permdist_b), `ETD') != 0 local ++bad
  if mreldif(e(perm_nullband), `ETB') != 0 local ++bad
  * the kernel's own quantile rows are the symmetric +/- q under permabsp
  tempname K12
  matrix `K12' = e(perm_kernel)
  forval k = 1/`=colsof(`K12')' {
    if `K12'[13, `k'] == 1 & (`K12'[4, `k'] != -`K12'[5, `k'] | `K12'[5, `k'] < 0) local ++bad
  }
  mata: st_numscalar("__bad", wdp3_reconx(0, 95, 1))
  if __bad != 0 local ++bad
  chk3tab , label(Q12d)
  * non-vacuity: at least one live p differs between the two conventions on these draws
  local ndiff 0
  local nlive 0
  forval k = 1/`=colsof(`ETP')' {
    if `ETP'[1, `k'] < . & e(perm_p)[1, `k'] < . {
      local ++nlive
      if `ETP'[1, `k'] != e(perm_p)[1, `k'] local ++ndiff
    }
  }
  if `nlive' < 3 | `ndiff' == 0 local ++bad
  * the legacy scalar carries the absolute-value p too, and the footer names the convention
  if e(perm_p_att) != e(wdidmainresults)[1, 4] local ++bad
  local lg = fileread("perm3_q12.txt")
  if strpos("`lg'", "absolute-value permutation p-values (permabsp)") == 0 local ++bad
  if strpos("`lg'", "equal-tailed permutation p-values") != 0 local ++bad
  if strpos("`lg'", "95% null band") == 0 local ++bad
  if `bad' {
    dis as error "Q12e FAIL: `bad' permabsp check(s) (live=`nlive' differing=`ndiff')"
    global pf = $pf + 1
  }
}
* exact enumeration on the Q4 fixture: the observed assignment is the unique maximum in
* |s| (the fixture is not an exact mirror - the 0.1*i unit trend breaks it), so the
* absolute-value p is the 1/M floor where the equal-tailed p was 2/20
mkperm3 6 3 3
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact permabsp
xrc Q12f 0 _rc
if _rc == 0 {
  local bad 0
  if abs(e(perm_p_att) - 0.05) > 1e-12 local ++bad
  if e(perm_exact) != 1 | e(perm_completed) != 20 local ++bad
  mata: st_numscalar("__bad", wdp3_reconx(1, 95, 1))
  if __bad != 0 local ++bad
  * symmetric about the estimate
  if reldif(_b[main_att] - e(wdidmainresults)[1, 5], e(wdidmainresults)[1, 6] - _b[main_att]) > 1e-12 local ++bad
  if `bad' {
    dis as error "Q12g FAIL: `bad' exact permabsp check(s): p=" e(perm_p_att)
    global pf = $pf + 1
  }
}
* permtstat + permtci + permabsp: the permtci interval IS the table's interval
mkperm3 8 4 3
set seed 1213
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(30) permtstat permtci permabsp
xrc Q12h 0 _rc
if _rc == 0 {
  local bad 0
  if "`e(perm_stat)'" != "t" | "`e(perm_pconv)'" != "absolute" local ++bad
  mata: st_numscalar("__bad", wdp3_reconx(0, 95, 1))
  if __bad != 0 local ++bad
  if e(perm_tlb_att) >= . local ++bad
  else if reldif(e(perm_tlb_att), e(wdidmainresults)[1, 5]) > 1e-12 | reldif(e(perm_tub_att), e(wdidmainresults)[1, 6]) > 1e-12 local ++bad
  if `bad' {
    dis as error "Q12i FAIL: `bad' permtstat+permtci+permabsp check(s)"
    global pf = $pf + 1
  }
}
* level(90) under permtci: the inverted-|t| interval follows level()
mkperm3 8 4 3
set seed 1213
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(30) permtstat permtci permabsp level(90)
xrc Q12j 0 _rc
if _rc == 0 {
  local bad 0
  mata: st_numscalar("__bad", wdp3_reconx(0, 90, 1))
  if __bad != 0 local ++bad
  if e(perm_tlb_att) >= . local ++bad
  else if reldif(e(perm_tlb_att), e(wdidmainresults)[1, 5]) > 1e-12 | reldif(e(perm_tub_att), e(wdidmainresults)[1, 6]) > 1e-12 local ++bad
  if `bad' {
    dis as error "Q12k FAIL: `bad' level(90) permtci check(s)"
    global pf = $pf + 1
  }
}
* stray
mkperm3 8 4 3
cap noi wooldid y i t ttre, att cluster(i) permabsp
xrc Q12l 198 _rc

* ================= Q13: sharp-null caveat direction (review 2026-09-05). In a staggered
* design the reassigned onsets move actual post-treatment outcomes into the draws'
* pre-treatment horizons, so the pre-period null distribution WIDENS with the true
* effect while the observed pre-period statistic (no pre-trend in the DGP) does not:
* the pre-period p rises and the null SD grows. Same seed, same draws, effect 0 vs 30.
* The footer and help used to say the p-values FALL; the footer text is scanned =========
cap program drop mkperm3s
program define mkperm3s
  args eff
  clear
  set seed 4242
  set sortseed 4242
  qui set obs 12
  qui gen i = _n
  qui gen ttre = cond(_n <= 4, 4, cond(_n <= 8, 6, .))
  qui expand 9
  qui bysort i: gen t = _n
  qui gen y = `eff'*(t >= ttre & ttre < .) + 0.1*i + rnormal()
end
foreach eff in 0 30 {
  mkperm3s `eff'
  set seed 113
  qui log using perm3_q13.txt, text replace name(p3)
  cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(80) jointtests
  local rc = _rc
  qui log close p3
  xrc Q13a_`eff' 0 `rc'
  if `rc' == 0 {
    local rn : roweq e(wdidmainresults)
    local k 0
    local q13k 0
    foreach r of local rn {
      local ++k
      if "`r'" == "ES_att" & e(wdidmainresults)[`k', 7] == -2 local q13k `k'
    }
    if `q13k' == 0 {
      dis as error "Q13b_`eff' FAIL: no ES_att row at relyr -2"
      global pf = $pf + 1
    }
    else {
      local q13p`eff' = e(perm_p)[1, `q13k']
      local q13sd`eff' = e(perm_sd)[1, `q13k']
      local q13j`eff' = e(joint_horiz)
    }
  }
}
if "`q13sd0'" != "" & "`q13sd30'" != "" {
  if !(`q13sd30' > 3*`q13sd0' & `q13p30' >= `q13p0' & `q13j30' >= `q13j0') {
    dis as error "Q13c FAIL: pre-period null SD `q13sd0' -> `q13sd30', p `q13p0' -> `q13p30', joint p `q13j0' -> `q13j30' (expected the null band to widen and the p-values to rise with the true effect)"
    global pf = $pf + 1
  }
  local lg = fileread("perm3_q13.txt")
  if strpos(`"`lg'"', "rise toward 1 as the true effect grows") == 0 | strpos(`"`lg'"', "p-values fall as the true post-treatment effect grows") > 0 {
    dis as error "Q13d FAIL: the footer caveat does not state the corrected (rising p-value) direction"
    global pf = $pf + 1
  }
}

* ================= Q14: rseed() under permute() (accepted since 2026-09-05): the same
* seed reproduces the draws bitwise whatever the session's RNG state, e(perm_seed) records
* it, and the integer validation now runs under permute() too =========================
mkperm3 8 4 3
set seed 1
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(25) rseed(77)
xrc Q14a 0 _rc
if _rc == 0 {
  tempname q14a
  matrix `q14a' = e(permdist_b)
  if "`e(perm_seed)'" != "77" {
    dis as error "Q14b FAIL: e(perm_seed) = `e(perm_seed)'"
    global pf = $pf + 1
  }
  set seed 999
  qui gen double q14junk = runiform()
  drop q14junk
  cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(25) rseed(77)
  xrc Q14c 0 _rc
  if _rc == 0 {
    if mreldif(e(permdist_b), `q14a') != 0 {
      dis as error "Q14d FAIL: rseed(77) did not reproduce the permutation draws"
      global pf = $pf + 1
    }
  }
}
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(25) rseed(1.5)
xrc Q14e 198 _rc
cap noi wooldid y i t ttre, att cluster(i) rseed(77)
xrc Q14f 198 _rc

* ================= Q15: the cluster() variable must nest i under permute() (gate added
* 2026-09-05): a cluster finer than the unit would give a unit two onsets in a
* replication - refused up front (459) with the reason; cluster == i and a coarser
* cluster (onset constant within it) still run =========================================
mkperm3 8 4 3
qui gen c2 = 2*i - (t <= 4)
qui log using perm3_q15.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(c2) permute(shuffle) permreps(10)
local rc = _rc
qui log close p3
xrc Q15a 459 `rc'
local lg = fileread("perm3_q15.txt")
if strpos(`"`lg'"', "requires the cluster() variable to nest i") == 0 {
  dis as error "Q15b FAIL: the finer-than-i refusal does not name the cause"
  global pf = $pf + 1
}
qui gen s = ceil(i/2)
cap noi wooldid y i t ttre, att cluster(s) permute(shuffle) permreps(10)
xrc Q15c 0 _rc
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(10)
xrc Q15d 0 _rc

* ================= Q16: exact enumeration and the fixed metric (fix 2026-09-05): when
* every assignment is enumerated the observed assignment is one of the draws, so the
* permutation covariance is formed from the draws alone (appending the comparator counted
* it twice). Reconstruct W_obs and the plain-share p of joint_horiz from the draws ======
mkperm3 6 3 3
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permexact jointtests
xrc Q16a 0 _rc
if _rc == 0 {
  local bad 0
  if e(perm_exact) != 1 | e(perm_completed) != e(perm_npossible) local ++bad
  local rn : roweq e(wdidmainresults)
  local k 0
  local q16cols
  foreach r of local rn {
    local ++k
    if "`r'" == "ES_att" & e(wdidmainresults)[`k', 7] < 0 local q16cols `q16cols' `k'
  }
  if wordcount("`q16cols'") < 2 local ++bad
  cap local jc = colnumb(e(perm_jointW), "joint_horiz")
  if _rc != 0 | `jc' >= . local ++bad
  if `bad' == 0 {
    mata: idx = strtoreal(tokens(st_local("q16cols"))); A = st_matrix("e(permdist_b)")[., idx']; ok = rowmissing(A) :== 0; A = select(A, ok); n = rows(A); d = st_matrix("e(perm_comparator_b)")[1, idx']'; dev = A :- mean(A); S = (dev' * dev) / (n - 1); W = d' * invsym(S) * d; Wr = st_matrix("e(perm_jointW)")[., `jc']; Wo = st_matrix("e(perm_jointWobs)")[1, `jc']; st_numscalar("__q16n", n); st_numscalar("__q16w", W); st_numscalar("__q16wo", Wo); st_numscalar("__q16p", sum(select(Wr, Wr :< .) :>= W)/n)
    if __q16n != e(perm_npossible) local ++bad
    if reldif(__q16w, __q16wo) > 1e-8 local ++bad
    if reldif(__q16p, e(joint_horiz)) > 1e-12 local ++bad
  }
  if `bad' {
    dis as error "Q16b FAIL: exact-mode joint test not reconstructable from the draws-only metric (n=`=__q16n', W=`=__q16w' vs posted `=__q16wo', p=`=__q16p' vs `=e(joint_horiz)')"
    global pf = $pf + 1
  }
}

* ================= Q11: a plain run after all of this is untouched by the engine ======
mkperm3 8 4 3
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
xrc Q11a 0 _rc
if _rc == 0 {
  if mreldif(e(wdidmainresults), `AN') != 0 {
    dis as error "Q11b FAIL: the analytic table moved after the permute runs"
    global pf = $pf + 1
  }
}

cap erase perm3_q1.dta
foreach f in perm3_q1.txt perm3_q1r.txt perm3_q5.txt perm3_q7.txt perm3_q7b.txt perm3_q12.txt perm3_q13.txt perm3_q14.txt perm3_q15.txt perm3_q16.txt perm3_q17.txt {
  cap erase `f'
}
cap erase perm3_q1_ES.dta
cap erase perm3_q8_ES.dta
cap erase perm3_q9_cfx.dta
local gl : dir . files "perm3_q*.gph"
foreach f of local gl {
  cap erase `f'
}

if $pf == 0 dis as result "PERMUTE3_ALL_PASS"
else dis as error "PERMUTE3_FAILS=$pf"
