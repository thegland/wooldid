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
cap program drop usegate3
program define usegate3
  * a permsave() write that fails is a caution since 2026-09-21, so rc 0 no longer proves the
  * file exists: name the miss as a FAIL and load a stand-in with the columns the checks
  * touch, so they fail on content rather than abort the suite r(601)
  args file label
  cap confirm file `file'.dta
  if _rc {
    dis as error "`label' FAIL: permsave() wrote no `file'.dta (e(perm_save_rc)=`e(perm_save_rc)')"
    global pf = $pf + 1
    clear
    qui set obs 1
    gen _rep = 0
    gen main_att = .
    gen onset_1 = .
    gen onset_8 = .
    exit
  }
  use `file', clear
end

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
  usegate3 perm3_q1 Q1s
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
* t-kernel reconstructs from e(permdist_b)/e(permdist_se). NOTE (audit 2026-09-13): on
* this fixture every column's signed tail counts are identical under the t and the b
* ranking, so the reconstruction here has NO power against a b-ranked engine; Q12i
* (permtstat + permabsp, seed 1213) is the sampled-mode leg that discriminates, and it
* carries the guard. ================================================================
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
* a user-typed fastnose with permanalyticse is refused as well: permanalyticse exists to
* compute the analytic standard errors that fastnose drops
mkperm3 8 4 3
qui log using perm3_q7k.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) permanalyticse fastnose
local rc = _rc
qui log close p3
xrc Q7k 198 `rc'
if strpos(fileread("perm3_q7k.txt"), "Error: fastnose is not compatible with permtstat, permtci or permanalyticse.") == 0 {
  dis as error "Q7k2 FAIL: the permanalyticse fastnose refusal did not name itself"
  global pf = $pf + 1
}
cap erase perm3_q7k.txt
* permsave() takes a file name only: permsave(name, replace) is refused before the loop,
* where it used to write a file literally called "name, replace.dta"
cap erase "perm3_q7l, replace.dta"
cap erase perm3_q7l.dta
qui log using perm3_q7l.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) permsave(perm3_q7l, replace)
local rc = _rc
qui log close p3
xrc Q7l 198 `rc'
if strpos(fileread("perm3_q7l.txt"), "Error: permsave() should only contain a file name.") == 0 | strpos(fileread("perm3_q7l.txt"), "Permutation replications") != 0 {
  dis as error "Q7l2 FAIL: the permsave() comma refusal did not name itself before the replication loop"
  global pf = $pf + 1
}
cap confirm file "perm3_q7l, replace.dta"
if _rc == 0 {
  dis as error "Q7l3 FAIL: permsave(perm3_q7l, replace) wrote a file named with the comma"
  global pf = $pf + 1
}
cap erase "perm3_q7l, replace.dta"
cap erase perm3_q7l.dta
cap erase perm3_q7l.txt

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
        * the saved data carry the table's own names (lb_90ci; saved as lb_95ci until
        * the 2026-09-05 fix)
        if lb_90ci[`r'] != `W'[`j', 5] | ub_90ci[`r'] != `W'[`j', 6] local ++nb
      }
    }
  }
  cap confirm variable lb_95ci
  if _rc == 0 local ++nb
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
  * non-vacuity (audit 2026-09-13): this is the sampled-mode leg with power against a b-ranked
  * engine - on at least one column the |t| tail count must differ from the |b| tail count
  * (2026-09-13: pretreatment_att and ES_att_pre1 give 8 vs 10); a re-seed that loses it fails here
  mata:
    A_ = st_matrix("e(permdist_b)"); S_ = st_matrix("e(permdist_se)")
    cb_ = st_matrix("e(perm_comparator_b)"); cs_ = st_matrix("e(perm_comparator_se)")
    nd_ = 0
    for (k_ = 1; k_ <= cols(A_); k_++) {
      if (cs_[k_] >= . | cs_[k_] <= 0 | cb_[k_] >= .) continue
      ok_ = (A_[., k_] :< .) :& (S_[., k_] :< .) :& (S_[., k_] :> 0)
      tt_ = select(A_[., k_] :/ S_[., k_], ok_)
      bb_ = select(A_[., k_], ok_)
      if (sum(abs(tt_) :>= abs(cb_[k_]/cs_[k_])) != sum(abs(bb_) :>= abs(cb_[k_]))) nd_++
    }
    st_numscalar("__q12idisc", nd_)
  end
  if __q12idisc == 0 {
    dis as error "Q12i FAIL (VACUOUS): the |t| and |b| tail counts agree on every column, so this leg could not tell a b-ranked engine from the t-ranked one; re-seed"
    local ++bad
  }
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
  * the footer line wraps at linesize, so the scan runs on the newline-collapsed log
  local lg = fileread("perm3_q13.txt")
  local lg = subinstr(`"`lg'"', char(13), "", .)
  local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
  local lg = subinstr(`"`lg'"', char(10), " ", .)
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
cap noi wooldid y i t ttre, att cluster(c2) permute(shuffle) permreps(20)
local rc = _rc
qui log close p3
xrc Q15a 459 `rc'
local lg = fileread("perm3_q15.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local lg = subinstr(`"`lg'"', char(10), " ", .)
if strpos(`"`lg'"', "strict subsets of your cluster() variable") == 0 {
  dis as error "Q15b FAIL: the finer-than-i refusal does not name the cause"
  global pf = $pf + 1
}
qui gen s = ceil(i/2)
cap noi wooldid y i t ttre, att cluster(s) permute(shuffle) permreps(20)
xrc Q15c 0 _rc
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20)
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

* ================= Q17: a subgroup value of 9+ digits would push the hidden analytic
* twin's name past 32 characters; the twin is skipped and the run completes (guard added
* 2026-09-05; e(wdidsgresults_<sg>) itself is 23 characters and unaffected) ============
mkperm3 8 4 3
* (double: a 9-digit value is not representable as the fixture's float sg)
qui gen double sg2 = 100000000 + sg
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) subgroups(sg2)
xrc Q17a 0 _rc
if _rc == 0 {
  cap confirm matrix e(wdidsgresults_100000001)
  local bad = _rc != 0
  cap confirm matrix e(perm_analytic_sgresults_100000001)
  if _rc == 0 local ++bad
  if `bad' {
    dis as error "Q17b FAIL: long subgroup value - visible table absent or the >32-char twin posted"
    global pf = $pf + 1
  }
}

* ================= Q18: dots(#) under permute(): one mark every # replications ========
mkperm3 8 4 3
qui log using perm3_q18.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) dots(10)
local rc = _rc
qui log close p3
xrc Q18a 0 `rc'
local lg = fileread("perm3_q18.txt")
local p0 = strpos(`"`lg'"', "Permutation replications (20): ")
if `p0' == 0 {
  dis as error "Q18b FAIL: no progress line in the log"
  global pf = $pf + 1
}
else {
  local seg = substr(`"`lg'"', `p0' + 31, 300)
  local seg = substr(`"`seg'"', 1, strpos(`"`seg'"', "done") - 1)
  local nd = length(`"`seg'"') - length(subinstr(`"`seg'"', ".", "", .))
  if `nd' != 2 {
    dis as error "Q18c FAIL: dots(10) with permreps(20) printed `nd' marks (expected 2)"
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
foreach f in perm3_q1.txt perm3_q1r.txt perm3_q5.txt perm3_q7.txt perm3_q7b.txt perm3_q12.txt perm3_q13.txt perm3_q14.txt perm3_q15.txt perm3_q16.txt perm3_q17.txt perm3_q18.txt {
  cap erase `f'
}
cap erase perm3_q1_ES.dta
cap erase perm3_q8_ES.dta
cap erase perm3_q9_cfx.dta
local gl : dir . files "perm3_q*.gph"
foreach f of local gl {
  cap erase `f'
}

* ================= Q19: the permreps resolution caution (F036, 2026-09-06): below 39
* replications the default equal-tailed p cannot reach .05, so a caution names the floor;
* it is silent at 40 and silent under permabsp at 20 (1/21 < .05) ======================
mkperm3 8 4 3
qui log using perm3_q19.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) rseed(19)
local rc = _rc
qui log close p3
xrc Q19a 0 `rc'
local lg = fileread("perm3_q19.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
local lg = subinstr(`"`lg'"', char(10), " ", .)
if strpos(`"`lg'"', "smallest attainable per-coefficient permutation p-value is 0.095 (2/(permreps + 1))") == 0 {
  dis as error "Q19b FAIL: no resolution caution at permreps(20)"
  global pf = $pf + 1
}
mkperm3 8 4 3
qui log using perm3_q19.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(40) rseed(19)
local rc = _rc
qui log close p3
xrc Q19c 0 `rc'
local lg = fileread("perm3_q19.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "smallest attainable per-coefficient permutation p-value") > 0 {
  dis as error "Q19d FAIL: the resolution caution fired at permreps(40)"
  global pf = $pf + 1
}
mkperm3 8 4 3
qui log using perm3_q19.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permabsp rseed(19)
local rc = _rc
qui log close p3
xrc Q19e 0 `rc'
local lg = fileread("perm3_q19.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "smallest attainable per-coefficient permutation p-value") > 0 {
  dis as error "Q19f FAIL: the resolution caution fired under permabsp at permreps(20)"
  global pf = $pf + 1
}
cap erase perm3_q19.txt

* ================= Q20: the permtstat-at-cluster(i) caution (F049, 2026-09-06): the
* studentized statistic is not protected when the cluster-robust se cannot see the treated
* clusters; the caution fires at cluster() = i with permtstat, and is silent with a coarser
* cluster or without permtstat ==========================================================
mkperm3 8 4 3
qui log using perm3_q20.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permtstat rseed(20)
local rc = _rc
qui log close p3
xrc Q20a 0 `rc'
local lg = fileread("perm3_q20.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "permtstat with cluster() at the i level") == 0 {
  dis as error "Q20b FAIL: no caution for permtstat at cluster() = i"
  global pf = $pf + 1
}
mkperm3 8 4 3
qui gen s = ceil(i/2)
qui log using perm3_q20.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(s) permute(shuffle) permreps(20) permtstat rseed(20)
local rc = _rc
qui log close p3
xrc Q20c 0 `rc'
local lg = fileread("perm3_q20.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "permtstat with cluster() at the i level") > 0 {
  dis as error "Q20d FAIL: the caution fired with a cluster coarser than i"
  global pf = $pf + 1
}
mkperm3 8 4 3
qui log using perm3_q20.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) rseed(20)
local rc = _rc
qui log close p3
xrc Q20e 0 `rc'
local lg = fileread("perm3_q20.txt")
local lg = subinstr(`"`lg'"', char(13), "", .)
local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
if strpos(`"`lg'"', "permtstat with cluster() at the i level") > 0 {
  dis as error "Q20f FAIL: the caution fired without permtstat"
  global pf = $pf + 1
}
cap erase perm3_q20.txt

* ================= Q21: the permutation reproducibility record (F059, 2026-09-06):
* e(perm_rngstate) is posted, replays the draws exactly, and permsave() carries the
* engine/seed/rngstate characteristics; the projection line prints ===================
mkperm3 8 4 3
set seed 2121
qui log using perm3_q21.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permsave(perm3_q21) replace dots(1)
local rc = _rc
qui log close p3
xrc Q21a 0 `rc'
local q21st `e(perm_rngstate)'
matrix Q21D1 = e(permdist_b)
local lg = fileread("perm3_q21.txt")
if strlen("`q21st'") < 100 | strpos(`"`lg'"', "projected for 20 replications)") == 0 {
  dis as error "Q21b FAIL: e(perm_rngstate) not posted or no runtime projection"
  global pf = $pf + 1
}
preserve
usegate3 perm3_q21 Q21s
local q21eng `_dta[wooldid_engine]'
local q21rl = strlen("`_dta[wooldid_rngstate]'")
restore
if "`q21eng'" != "permute" | `q21rl' < 100 {
  dis as error "Q21c FAIL: permsave() characteristics missing (engine=`q21eng' rnglen=`q21rl')"
  global pf = $pf + 1
}
mkperm3 8 4 3
set rngstate `q21st'
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20)
xrc Q21d 0 _rc
if mreldif(Q21D1, e(permdist_b)) != 0 {
  dis as error "Q21e FAIL: set rngstate e(perm_rngstate) does not replay the permutation draws"
  global pf = $pf + 1
}
cap erase perm3_q21.txt
cap erase perm3_q21.dta

* ================= Q22 (F052/F064, 2026-09-06): the p-value conventions and the
* Monte-Carlo precision record. From the same draws: permrawp = c/n |s| share with no
* add-one; permside(left|right) = the one-sided add-one shares, which e(perm_p_lo) /
* e(perm_p_hi) carry under every convention and whose 2*min is the default p; e(perm_p_se)
* and e(perm_p_ci) follow the documented formulas; the conventions are mutually exclusive
mkperm3 8 4 3
qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(40) rseed(22)
local q22n = e(perm_completed)
matrix Q22D = e(permdist_b)
matrix Q22C = e(perm_comparator_b)
local q22p0 = e(perm_p)[1, 1]
local q22lo = e(perm_p_lo)[1, 1]
local q22hi = e(perm_p_hi)[1, 1]
local q22se = e(perm_p_se)[1, 1]
local q22cl = e(perm_p_ci)[1, 1]
local q22ch = e(perm_p_ci)[2, 1]
mata: D = st_matrix("Q22D")[., 1]; c = st_matrix("Q22C")[1, 1]; D = select(D, D :< .); n = rows(D); m0 = sum(D :<= c); m1 = sum(D :>= c); ma = sum(abs(D) :>= abs(c)); st_numscalar("q22_n", n); st_numscalar("q22_plo", (1 + m0)/(n + 1)); st_numscalar("q22_phi", (1 + m1)/(n + 1)); st_numscalar("q22_raw", ma/n); cc = min((m0, m1)); q = cc/n; st_numscalar("q22_se", 2*sqrt(q*(1 - q)/n)); st_numscalar("q22_cl", min((1, 2*(cc <= 0 ? 0 : invibeta(cc, n - cc + 1, .025))))); st_numscalar("q22_ch", min((1, 2*(cc >= n ? 1 : invibeta(cc + 1, n - cc, .975)))))
di "Q22 n=`q22n'/" q22_n " p=" %10.8f `q22p0' " 2min=" %10.8f min(1, 2*min(q22_plo, q22_phi)) " lo=" %10.8f `q22lo' "/" %10.8f q22_plo " hi=" %10.8f `q22hi' "/" %10.8f q22_phi " se=" %10.8f `q22se' "/" %10.8f q22_se " ci=[" %8.6f `q22cl' "," %8.6f `q22ch' "] hand=[" %8.6f q22_cl "," %8.6f q22_ch "] pconv=`e(perm_pconv)'"
if reldif(`q22p0', min(1, 2*min(q22_plo, q22_phi))) > 1e-12 | reldif(`q22lo', q22_plo) > 1e-12 | reldif(`q22hi', q22_phi) > 1e-12 ///
   | reldif(`q22se', q22_se) > 1e-10 | reldif(`q22cl', q22_cl) > 1e-8 | reldif(`q22ch', q22_ch) > 1e-8 | "`e(perm_pconv)'" != "equaltailed" {
  dis as error "Q22a FAIL: the tail-share / precision record does not reconstruct from the draws"
  global pf = $pf + 1
}
mkperm3 8 4 3
qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(40) rseed(22) permrawp
local q22r = e(perm_p)[1, 1]
di "Q22 permrawp p=" %10.8f `q22r' " hand c/n=" %10.8f q22_raw " pconv=`e(perm_pconv)' absp=" e(perm_absp)
if reldif(`q22r', q22_raw) > 1e-12 | "`e(perm_pconv)'" != "raw" | e(perm_absp) != 0 {
  dis as error "Q22b FAIL: permrawp is not the plain c/n absolute share"
  global pf = $pf + 1
}
mkperm3 8 4 3
qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(40) rseed(22) permside(left)
local q22l = e(perm_p)[1, 1]
local q22lconv `e(perm_pconv)'
mkperm3 8 4 3
qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(40) rseed(22) permside(right)
local q22rt = e(perm_p)[1, 1]
di "Q22 left p=" %10.8f `q22l' " (" %10.8f q22_plo ") right p=" %10.8f `q22rt' " (" %10.8f q22_phi ") convs=`q22lconv'/`e(perm_pconv)'"
if reldif(`q22l', q22_plo) > 1e-12 | reldif(`q22rt', q22_phi) > 1e-12 | "`q22lconv'" != "left" | "`e(perm_pconv)'" != "right" {
  dis as error "Q22c FAIL: permside() does not report the one-sided add-one share"
  global pf = $pf + 1
}
mkperm3 8 4 3
cap qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permside(up)
xrc Q22d 198 _rc
mkperm3 8 4 3
cap qui wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permrawp permabsp
xrc Q22e 198 _rc
mkperm3 8 4 3
cap qui wooldid y i t ttre, att cluster(i) permside(left)
xrc Q22f 198 _rc
* the bootstrap analogue: hidden e(boot_ptails) reconstructs from e(bootdist)
mkperm3 8 4 3
qui wooldid y i t ttre, att cluster(i) bootstrap(standard) reps(40) rseed(22)
matrix Q22B = e(bootdist)
matrix Q22T = e(boot_ptails)
mata: D = st_matrix("Q22B")[., 1]; D = select(D, D :< .); n = rows(D); m0 = sum(D :<= 0); m1 = sum(D :>= 0); cc = min((m0, m1)); q = cc/n; st_numscalar("q22b_lo", m0/n); st_numscalar("q22b_hi", m1/n); st_numscalar("q22b_se", 2*sqrt(q*(1 - q)/n)); st_numscalar("q22b_cl", min((1, 2*(cc <= 0 ? 0 : invibeta(cc, n - cc + 1, .025))))); st_numscalar("q22b_ch", min((1, 2*(cc >= n ? 1 : invibeta(cc + 1, n - cc, .975)))))
di "Q22 boot_ptails: lo=" Q22T[1,1] "/" q22b_lo " hi=" Q22T[2,1] "/" q22b_hi " se=" Q22T[3,1] "/" q22b_se " ci=[" Q22T[4,1] "," Q22T[5,1] "] hand=[" q22b_cl "," q22b_ch "]"
if reldif(Q22T[1,1], q22b_lo) > 1e-12 | reldif(Q22T[2,1], q22b_hi) > 1e-12 | reldif(Q22T[3,1], q22b_se) > 1e-10 | reldif(Q22T[4,1], q22b_cl) > 1e-8 | reldif(Q22T[5,1], q22b_ch) > 1e-8 {
  dis as error "Q22g FAIL: e(boot_ptails) does not reconstruct from e(bootdist)"
  global pf = $pf + 1
}

* ================= Q23 (F065, 2026-09-06): the max-|studentized| step-down over the post
* horizons reconstructs by name from e(permdist_b) and e(perm_comparator_b): fixed-metric
* studentization (SD of draws + comparator), raw per-horizon |t| shares, single-step max
* and the monotone step-down in decreasing |t_obs| order, all add-one; the footer names it
* permstepdown added 2026-09-10 (PS-45, author ruling 2026-09-07): the procedure used to
* run on every permute() run with two or more post horizons and is now opt-in. The numbers
* it produces did not change, so this leg still pins them
mkperm3 8 4 3
qui log using perm3_q23.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(2) permute(shuffle) permreps(40) rseed(23) permstepdown
local rc = _rc
qui log close p3
xrc Q23a 0 `rc'
cap confirm matrix e(perm_stepdown)
xrc Q23b 0 _rc
if _rc == 0 {
  matrix Q23S = e(perm_stepdown)
  matrix Q23D = e(permdist_b)
  matrix Q23C = e(perm_comparator_b)
  local q23sc : colnames Q23S
  local q23dc : colnames Q23D
  local q23idx
  foreach c of local q23sc {
    local q23idx `q23idx' `: list posof "`c'" in q23dc'
  }
  mata:
  idx = strtoreal(tokens(st_local("q23idx")))
  A = st_matrix("Q23D")[., idx']
  cm = st_matrix("Q23C")[1, idx']
  ok = rowmissing(A) :== 0
  D = select(A, ok)
  n = rows(D)
  CC = D \ cm
  sdv = sqrt(diagonal(variance(CC)))'
  TT = abs(D) :/ sdv
  tob = abs(cm) :/ sdv
  H = cols(TT)
  praw = J(1, H, .)
  pss = J(1, H, .)
  psd = J(1, H, .)
  for (h = 1; h <= H; h++) praw[h] = (1 + sum(TT[., h] :>= tob[h]))/(n + 1)
  mx = rowmax(TT)
  for (h = 1; h <= H; h++) pss[h] = (1 + sum(mx :>= tob[h]))/(n + 1)
  ord = order(-tob', 1)'
  pp = 0
  for (j = 1; j <= H; j++) {
    mx = rowmax(TT[., ord[j..H]])
    pj = (1 + sum(mx :>= tob[ord[j]]))/(n + 1)
    if (pj < pp) pj = pp
    psd[ord[j]] = pj
    pp = pj
  }
  st_matrix("Q23H", tob \ praw \ pss \ psd)
  st_numscalar("q23_n", n)
  end
  di "Q23 H=" colsof(Q23S) " n=" q23_n " maxdiff=" mreldif(Q23S, Q23H) " stepdown=[" Q23S[4,1] " " Q23S[4,2] "] raw=[" Q23S[2,1] " " Q23S[2,2] "]"
  local lg = fileread("perm3_q23.txt")
  local lg = subinstr(`"`lg'"', char(13), "", .)
  local lg = subinstr(`"`lg'"', char(10) + "> ", "", .)
  if colsof(Q23S) != 2 | mreldif(Q23S, Q23H) > 1e-12 | Q23S[4,1] < Q23S[2,1] | Q23S[4,2] < Q23S[2,2] | strpos(`"`lg'"', "max-|studentized| step-down") == 0 {
    dis as error "Q23c FAIL: the step-down record does not reconstruct by name, is not >= the raw p, or is not printed"
    global pf = $pf + 1
  }
}
cap erase perm3_q23.txt

* ================= Q24 (CS-87, 2026-09-06): makeplots with no event-study window under
* permute() no longer dies r(2000) in the deferred plot call ===========================
mkperm3 8 4 3
set graphics off
cap qui wooldid y i t ttre, att cluster(i) makeplots permute(shuffle) permreps(20) rseed(24)
xrc Q24a 0 _rc
cap confirm matrix e(wdidmainresults)
xrc Q24b 0 _rc


* ================= Q25-Q28 (bug sweep BS-20/22/24/25, 2026-09-06): the follow-ups
cap program drop scan3
program define scan3, rclass
  args lf needle
  tempname fh
  local found 0
  file open `fh' using "`lf'", read text
  file read `fh' line
  while r(eof) == 0 {
    if strpos(`"`macval(line)'"', `"`needle'"') != 0 local ++found
    file read `fh' line
  }
  file close `fh'
  return scalar found = `found'
end
* Q25 (BS-20): the resolution caution follows the convention in force - at permreps(20)
* the equal-tailed floor 2/21 = .095 draws it, the one-sided/absolute floor 1/21 does
* not, and the plain share has no floor; the "cannot reject" claim and the printed p
* must never disagree in one output
mkperm3 8 4 3
tempfile q25
save `q25'
local q25n
foreach conv in "" "permabsp" "permside(right)" "permrawp" {
  use `q25', clear
  qui log using perm3_q25.txt, text replace name(p3)
  cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) rseed(25) `conv'
  local rc = _rc
  qui log close p3
  xrc "Q25 `conv'" 0 `rc'
  scan3 perm3_q25.txt "smallest attainable per-coefficient permutation p-value"
  local q25n `q25n' `r(found)'
  local q25p_`=cond("`conv'" == "", "eq", cond("`conv'" == "permabsp", "abs", cond("`conv'" == "permrawp", "raw", "right")))' = e(perm_p)[1,1]
}
cap erase perm3_q25.txt
dis "Q25 caution counts (eq abs right raw) = `q25n'; p: eq=`q25p_eq' abs=`q25p_abs' right=`q25p_right' raw=`q25p_raw'"
if "`q25n'" != "1 0 0 0" | `q25p_right' >= .05 | `q25p_raw' >= .05 {
  dis as error "Q25 FAIL: the floor caution does not follow the convention (counts `q25n'; right p `q25p_right', raw p `q25p_raw')"
  global pf = $pf + 1
}
* Q26 (BS-22): under permexact the shares are exact - e(perm_p_se) and e(perm_p_ci) are
* missing and the footer says so; a sampled companion reports both (the SE is 0 when no
* draw reaches the observed statistic, as Stata's permute reports), and the interval
* brackets the p
use `q25', clear
qui log using perm3_q26.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact
local rc = _rc
qui log close p3
xrc Q26a 0 `rc'
matrix Q26S = e(perm_p_se)
matrix Q26C = e(perm_p_ci)
matrix Q26P = e(perm_p)
mata: st_local("q26miss", strofreal(allof(st_matrix("Q26S"), .) & allof(st_matrix("Q26C"), .)))
* BS-23 (author ruling 2026-09-07) re-pinned the wording: the footer used to say the
* enumeration was "exhaustive", which is false on an incomplete one, and now says the
* assignments are enumerated rather than sampled - a claim that holds either way
scan3 perm3_q26.txt "enumerated rather than sampled, so the p-values are exact shares"
local q26f = r(found)
scan3 perm3_q26.txt "Monte-Carlo precision of the p-values"
local q26g = r(found)
cap erase perm3_q26.txt
use `q25', clear
qui log using perm3_q26.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(40) rseed(26)
local rc = _rc
qui log close p3
xrc Q26b 0 `rc'
matrix Q26S2 = e(perm_p_se)
matrix Q26C2 = e(perm_p_ci)
matrix Q26P2 = e(perm_p)
scan3 perm3_q26.txt "Monte-Carlo precision of the p-values"
local q26h = r(found)
cap erase perm3_q26.txt
dis "Q26 exact: p=" Q26P[1,1] " se/ci all missing=`q26miss' footer exact=`q26f' mc=`q26g' | sampled: p=" Q26P2[1,1] " se=" Q26S2[1,1] " ci=[" Q26C2[1,1] "," Q26C2[2,1] "] footer mc=`q26h'"
if `q26miss' != 1 | `q26f' != 1 | `q26g' != 0 | Q26P[1,1] >= . | Q26S2[1,1] >= . | Q26S2[1,1] < 0 | Q26C2[1,1] > Q26P2[1,1] | Q26C2[2,1] < Q26P2[1,1] | `q26h' != 1 {
  dis as error "Q26 FAIL: Monte-Carlo precision under permexact / sampled draws"
  global pf = $pf + 1
}
* Q27 (BS-24): the step-down footer names its own metric (p_raw), and the kernel row is
* the documented p_raw
use `q25', clear
qui log using perm3_q27.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(2) permute(shuffle) permreps(40) rseed(23) permstepdown
local rc = _rc
qui log close p3
xrc Q27a 0 `rc'
scan3 perm3_q27.txt "horizon: p_raw -> step-down p, where p_raw is the two-sided |t| share"
local q27a = r(found)
scan3 perm3_q27.txt "raw p -> step-down"
local q27b = r(found)
cap erase perm3_q27.txt
local q27rn : rownames e(perm_stepdown)
dis "Q27 footer new=`q27a' old=`q27b' rows=[`q27rn']"
if `q27a' != 1 | `q27b' != 0 | !`: list posof "p_raw" in q27rn' {
  dis as error "Q27 FAIL: step-down footer / kernel row"
  global pf = $pf + 1
}
* Q28 (BS-25): a plain-share p of exactly 0 (permrawp, sampled draws, a strong effect on
* a large assignment space) draws the resolution note; the add-one conventions and an
* exact enumeration do not (an exact share is never below 1/M since the observed
* assignment is a member)
mkperm3 20 8 3
tempfile q28
save `q28'
qui log using perm3_q28.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(99) rseed(28) permrawp
local rc = _rc
qui log close p3
xrc Q28a 0 `rc'
local q28p = e(perm_p)[1,1]
scan3 perm3_q28.txt "p = 0 above means p <"
local q28n = r(found)
local q28x
foreach conv in "" "permabsp" {
  use `q28', clear
  qui log using perm3_q28.txt, text replace name(p3)
  cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(99) rseed(28) `conv'
  local rc = _rc
  qui log close p3
  xrc "Q28 `conv'" 0 `rc'
  scan3 perm3_q28.txt "p = 0 above means p <"
  local q28x `q28x' `r(found)'
}
use `q25', clear
qui log using perm3_q28.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permexact permrawp
local rc = _rc
qui log close p3
xrc Q28d 0 `rc'
mata: st_local("q28min", strofreal(min(st_matrix("e(perm_p)"))))
local q28M = e(perm_completed)
scan3 perm3_q28.txt "p = 0 above means p <"
local q28e = r(found)
cap erase perm3_q28.txt
dis "Q28 permrawp sampled: p=`q28p' note=`q28n' | default/permabsp notes=`q28x' | exact permrawp: min p=`q28min' completed=`q28M' note=`q28e'"
if `q28p' != 0 | `q28n' != 1 | "`q28x'" != "0 0" | `q28min' < 1/`q28M' - 1e-6 | `q28e' != 0 {
  dis as error "Q28 FAIL: the p = 0 resolution note under permrawp"
  global pf = $pf + 1
}

* ---- Q29 (BS-23): an INCOMPLETE permexact enumeration keeps the exact share ----
* Units 7-8 are observed only from t = 4, so any assignment that treats one of them loses
* pre-horizons and that element's reference set shrinks. Before BS-23 those elements fell
* back to the add-one form (1 + m)/(n + 1); the observed assignment is one of the
* enumerated columns whether or not others failed, so the share stays exact over the
* assignments that completed. What incompleteness costs is resolution, not the convention.
mkperm3 8 2 3
qui drop if inlist(i, 7, 8) & t < 4
tempfile q29d
qui save `q29d'
use `q29d', clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(2) permute(shuffle) permexact jointtests
local q29rc = _rc
matrix Q29D = e(permdist_b)
matrix Q29C = e(perm_comparator_b)
matrix Q29P = e(perm_p)
matrix Q29LO = e(perm_p_lo)
matrix Q29HI = e(perm_p_hi)
matrix Q29N = e(perm_n)
matrix Q29SE = e(perm_p_se)
matrix Q29CI = e(perm_p_ci)
local q29inc = "`e(perm_incomplete)'"
local q29add = "`e(perm_addone)'"
local q29np = e(perm_npossible)
local q29jh = e(joint_horiz)
matrix Q29W = e(perm_jointW)
matrix Q29WOBS = e(perm_jointWobs)
* T1 membership: some replication row reproduces the comparator bitwise
mata: st_local("q29t1", strofreal(sum(rowsum(st_matrix("Q29D") :== st_matrix("Q29C")) :== cols(st_matrix("Q29D"))) >= 1))
* T2 the convention, recomputed from the posted draws
mata {
  D = st_matrix("Q29D"); C = st_matrix("Q29C")
  bad = 0; minn = .; anyshort = 0
  for (k = 1; k <= cols(D); k++) {
    v = select(D[., k], D[., k] :< .)
    nk = rows(v)
    if (nk < 2 | C[1, k] >= .) continue
    m1 = sum(v :>= C[1, k]); m0 = sum(v :<= C[1, k])
    if (abs(st_matrix("Q29HI")[1, k] - m1/nk) > 1e-12) bad = bad + 1
    if (abs(st_matrix("Q29LO")[1, k] - m0/nk) > 1e-12) bad = bad + 1
    if (abs(st_matrix("Q29P")[1, k] - min((1, 2*min((m0, m1))/nk))) > 1e-12) bad = bad + 1
    if (st_matrix("Q29N")[1, k] != nk) bad = bad + 1
    if (nk < strtoreal(st_local("q29np"))) anyshort = anyshort + 1
    if (m1 == 1 & abs(st_matrix("Q29HI")[1, k] - 1/nk) > 1e-12) bad = bad + 1
    if (nk < minn) minn = nk
  }
  st_local("q29bad", strofreal(bad))
  st_local("q29short", strofreal(anyshort))
  st_local("q29minn", strofreal(minn))
}
* T5 the enumeration carries no Monte-Carlo error, incomplete columns included
mata: st_local("q29t5", strofreal(allof(st_matrix("Q29SE"), .) & allof(st_matrix("Q29CI"), .)))
* T7 the joint test is the plain share over the reduced set
mata {
  w = select(st_matrix("Q29W")[., 1], st_matrix("Q29W")[., 1] :< .)
  nj = rows(w)
  // read Wobs from the matrix, never through a Stata local: a local stringifies at the
  // default precision, and the observed assignment's own row carries W == Wobs exactly,
  // so a rounded comparator moves it across the >= boundary and costs the count one row
  nex = sum(w :>= st_matrix("Q29WOBS")[1, 1])
  st_local("q29nj", strofreal(nj))
  st_local("q29jp", strofreal(nex/nj, "%21x"))
}
* T8 no leakage into sampled mode: M is missing there, so the add-one form must stand
use `q29d', clear
cap qui wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(2) permute(shuffle) permreps(200) rseed(29)
local q29rc8 = _rc
matrix Q29D8 = e(permdist_b)
matrix Q29C8 = e(perm_comparator_b)
matrix Q29P8 = e(perm_p)
matrix Q29SE8 = e(perm_p_se)
mata {
  D = st_matrix("Q29D8"); C = st_matrix("Q29C8")
  bad8 = 0
  for (k = 1; k <= cols(D); k++) {
    v = select(D[., k], D[., k] :< .)
    nk = rows(v)
    if (nk < 2 | C[1, k] >= .) continue
    m1 = sum(v :>= C[1, k]); m0 = sum(v :<= C[1, k])
    if (abs(st_matrix("Q29P8")[1, k] - min((1, 2*(min((m0, m1)) + 1)/(nk + 1)))) > 1e-12) bad8 = bad8 + 1
  }
  st_local("q29bad8", strofreal(bad8))
  st_local("q29se8", strofreal(allof(st_matrix("Q29SE8"), .)))
}
dis "Q29 BS-23: rc=`q29rc' npossible=`q29np' shortest n=`q29minn' short cols=`q29short' recompute mismatches=`q29bad' | incomplete=[`q29inc'] addone=[`q29add'] member=`q29t1' no-MC-error=`q29t5' | joint n=`q29nj' p=" %10.8f `q29jh' " vs share " %10.8f `=`q29jp'' " | sampled rc=`q29rc8' add-one mismatches=`q29bad8' se all missing=`q29se8'"
if `q29rc' != 0 | `q29t1' != 1 | `q29bad' != 0 | `q29short' == 0 | "`q29inc'" == "" | "`q29add'" != "" | `q29t5' != 1 {
  dis as error "Q29 FAIL: an incomplete permexact enumeration did not keep the exact share"
  global pf = $pf + 1
}
if `q29nj' >= `q29np' | abs(`q29jh' - `=`q29jp'') > 1e-12 {
  dis as error "Q29b FAIL: the joint test is not the plain share over the reduced reference set (n=`q29nj', p=`q29jh')"
  global pf = $pf + 1
}
if `q29rc8' != 0 | `q29bad8' != 0 | `q29se8' != 0 {
  dis as error "Q29c FAIL: membership leaked into sampled mode (mismatches=`q29bad8', se all missing=`q29se8')"
  global pf = $pf + 1
}

* ================= Q30 (PS-45, author ruling 2026-09-07): the step-down is OPT-IN, and the
* pre-treatment horizons are their own family - adjusted among themselves, reported and
* stored apart, never pooled with the post-treatment ones (pooling would change every
* post-treatment adjusted p by enlarging its family)
mkperm3 8 4 3
tempfile q30
qui save `q30'
qui log using perm3_q30a.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(40) rseed(30)
local q30rc = _rc
qui log close p3
xrc Q30a 0 `q30rc'
local q30off 0
foreach m in perm_stepdown perm_stepdown_pre perm_stepdown_i perm_Cstepdown {
  cap confirm matrix e(`m')
  if _rc == 0 local ++q30off
}
scan3 perm3_q30a.txt "Multiplicity across the"
local q30foot = r(found)
cap erase perm3_q30a.txt
if `q30off' != 0 | `q30foot' != 0 {
  dis as error "Q30a FAIL: the step-down ran without permstepdown (matrices=`q30off', footer=`q30foot')"
  global pf = $pf + 1
}
use `q30', clear
qui log using perm3_q30b.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(40) rseed(30) permstepdown
local q30rcb = _rc
qui log close p3
xrc Q30b 0 `q30rcb'
cap confirm matrix e(perm_stepdown)
local q30hp = _rc == 0
cap confirm matrix e(perm_stepdown_pre)
local q30hr = _rc == 0
local q30pc
local q30rcn
if `q30hp' local q30pc : colnames e(perm_stepdown)
if `q30hr' local q30rcn : colnames e(perm_stepdown_pre)
* every post family column must be a post horizon and every pre family column a pre one,
* and the two lists must not overlap - that is what "separate families" means here
local q30bad 0
foreach c of local q30pc {
  if strpos("`c'", "post") == 0 local ++q30bad
  if `: list posof "`c'" in q30rcn' local ++q30bad
}
foreach c of local q30rcn {
  if strpos("`c'", "pre") == 0 local ++q30bad
}
scan3 perm3_q30b.txt "across the post-treatment att horizons"
local q30fp = r(found)
scan3 perm3_q30b.txt "across the pre-treatment att horizons"
local q30fr = r(found)
cap erase perm3_q30b.txt
dis "Q30 post=[`q30pc'] pre=[`q30rcn'] footers post=`q30fp' pre=`q30fr' bad=`q30bad'"
if `q30hp' != 1 | `q30hr' != 1 | `q30bad' != 0 | `: word count `q30pc'' < 2 | `: word count `q30rcn'' < 2 | `q30fp' != 1 | `q30fr' != 1 {
  dis as error "Q30b FAIL: post matrix=`q30hp' pre matrix=`q30hr' mixed columns=`q30bad' footers=`q30fp'/`q30fr'"
  global pf = $pf + 1
}
* the pre family is adjusted among ITS OWN members: with two pre horizons the largest
* adjusted p equals the single-step p of that family, which cannot be read off the post one
if `q30hr' {
  matrix Q30R = e(perm_stepdown_pre)
  local q30mono 1
  forval c = 1/`=colsof(Q30R)' {
    if Q30R[4, `c'] < Q30R[2, `c'] - 1e-12 local q30mono 0
  }
  if `q30mono' != 1 {
    dis as error "Q30c FAIL: a pre-treatment adjusted p is below its own raw p"
    global pf = $pf + 1
  }
}

* ---- Q31 P7 accident 1 (author ruling 2026-09-08): permute() has its OWN exchangeability
* caution, not the bootstrap's few-cluster ones ported across. Two contracts. (a) With
* cluster() at i and the DEFAULT statistic the run says the coefficient carries no protection
* against a treated cluster of unusual variance and names permtstat; under permtstat that line
* is replaced by the permtstat-specific one, never both. (b) The few-ever-treated-cluster line
* is about how many DISTINCT assignments exist, so it must NOT fire on a fixture with many
* treated clusters - a caution that always prints teaches users to ignore it.
use "`DATA'/dta_binary", clear
qui log using perm3_q31a.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20) rseed(31) dots(0)
local q31rca = _rc
qui log close p3
qui log using perm3_q31b.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) esprelength(1) espostlength(1) permute(shuffle) permreps(20) rseed(31) permtstat dots(0)
local q31rcb = _rc
qui log close p3
tempname fh
local q31a_def 0
local q31a_ts 0
local q31a_few 0
file open `fh' using perm3_q31a.txt, read text
file read `fh' line
while r(eof) == 0 {
  if strpos(`"`line'"', "the permuted statistic is the coefficient itself") > 0 local q31a_def 1
  if strpos(`"`line'"', "Caution: permtstat with cluster() at the i level") > 0 local q31a_ts 1
  if strpos(`"`line'"', "of the resampling clusters are ever treated") > 0 local q31a_few 1
  file read `fh' line
}
file close `fh'
local q31b_def 0
local q31b_ts 0
file open `fh' using perm3_q31b.txt, read text
file read `fh' line
while r(eof) == 0 {
  if strpos(`"`line'"', "the permuted statistic is the coefficient itself") > 0 local q31b_def 1
  if strpos(`"`line'"', "Caution: permtstat with cluster() at the i level") > 0 local q31b_ts 1
  file read `fh' line
}
file close `fh'
if `q31rca' != 0 | `q31rcb' != 0 | `q31a_def' != 1 | `q31a_ts' != 0 | `q31a_few' != 0 | `q31b_def' != 0 | `q31b_ts' != 1 {
  dis as error "Q31 FAIL: default run (coef caution=`q31a_def' ts caution=`q31a_ts' few=`q31a_few'), permtstat run (coef=`q31b_def' ts=`q31b_ts')"
  global pf = $pf + 1
}
else dis as result "Q31 PASS: the coefficient-statistic caution and the permtstat one are alternatives, and the few-treated-cluster line stays silent on an 80-cluster fixture"

* ---- Q32 PS-46 (author ruling 2026-09-10): under the fast-mode default the run posts its
* permutation inference under the derivable names estout reads AHEAD of any covariance, so a
* coefficient table built with the standard tools carries this run's real numbers. The
* contract, both halves: e(p) and e(ci_l)/e(ci_u) are posted and EQUAL the perm_ matrices;
* and nothing is invented - no e(V), no e(se), so test/lincom/_se[] still refuse. Under
* permanalyticse an e(V) exists and the tools work from it, so the derivable names are not
* posted there and the two routes cannot disagree. ----
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(20) rseed(32) dots(0)
local q32bad 0
foreach nm in p ci_l ci_u {
  cap confirm matrix e(`nm')
  if _rc != 0 local ++q32bad
}
foreach nm in V se {
  cap confirm matrix e(`nm')
  if _rc == 0 local ++q32bad
}
if `q32bad' == 0 {
  tempname Q32P Q32L Q32U Q32PP Q32CI
  matrix `Q32P' = e(p)
  matrix `Q32L' = e(ci_l)
  matrix `Q32U' = e(ci_u)
  matrix `Q32PP' = e(perm_p)
  matrix `Q32CI' = e(perm_ci)
  if colsof(`Q32P') != colsof(`Q32PP') local ++q32bad
  else {
    forval c = 1/`=colsof(`Q32P')' {
      if `Q32P'[1,`c'] != `Q32PP'[1,`c'] local ++q32bad
      if `Q32L'[1,`c'] != `Q32CI'[1,`c'] local ++q32bad
      if `Q32U'[1,`c'] != `Q32CI'[2,`c'] local ++q32bad
    }
  }
}
cap test main_att = 0
local q32trc = _rc
if `q32trc' == 0 local ++q32bad
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(20) rseed(32) permanalyticse dots(0)
cap confirm matrix e(p)
local q32ase = (_rc == 0)
cap confirm matrix e(V)
local q32asev = (_rc == 0)
if `q32ase' != 0 | `q32asev' != 1 local ++q32bad
if `q32bad' != 0 {
  dis as error "Q32 FAIL: fast-mode posting bad=`q32bad' (test rc=`q32trc'; permanalyticse e(p)=`q32ase' e(V)=`q32asev')"
  global pf = $pf + 1
}
else dis as result "Q32 PASS: fast mode posts e(p)/e(ci_l)/e(ci_u) equal to the perm_ matrices, posts no e(V) or e(se), and test still refuses (rc=`q32trc')"

* ------------------------------------------- Q33: D-1, cfxplotnose under permute() (2026-09-11)
* cfxplotnose used to truncate the cfx matrices to one column unconditionally, while
* fastnose carved permute() out. So under permute() every cfx matrix came back one column
* wide, wooldid_permcellsx skipped all of them on its six-column test, ndone stayed 0 - and
* therefore the drift caution never printed and e(perm_cfx_blank) was never posted: analytic
* plot cells under a permutation header with nothing said. cfxplotnose now takes the same
* carve-out, so the cells carry permutation inference.
clear
set obs 20
gen i = _n
gen ttre = cond(i <= 8, 5, cond(i <= 15, 7, .))
expand 10
bysort i: gen t = _n
set seed 4411
gen dose = runiform() + 0.5
gen y = 1.1*(t >= ttre & ttre < .)*dose + 0.08*i + 0.04*t + rnormal()
qui wooldid y i t ttre, att cluster(i) contreat(dose) contreatcontrols(basic) ///
  makecfxplots lattice(all 3) cfxplotnose permute(shuffle) permreps(20) rseed(77) ///
  esprelength(1) espostlength(1) dots(0)
local q33nm : word 1 of `e(wdid_cfxmats)'
local q33c = .
cap local q33c = colsof(e(`q33nm'))
cap confirm matrix e(permdist_`q33nm')
local q33pd = (_rc == 0)
* the rewiring must have filled the inference columns from the replications
local q33p = .
cap local q33p = e(`q33nm')[1, 4]
if "`q33nm'" == "" | `q33c' != 6 | `q33pd' != 1 | `q33p' >= . {
  dis as error "Q33 FAIL: permute() + cfxplotnose must keep the 6-column cfx shape and permute it - matrix=[`q33nm'], colsof=`q33c', permdist posted=`q33pd', first cell p=`q33p'"
  global pf = $pf + 1
}
else dis as result "Q33 PASS: cfxplotnose no longer truncates under permute() - `q33nm' has `q33c' columns, a permutation p (`q33p') and its own permdist_"

* ================= Q34: aweights under permute() - the recon identity on a weighted
* replication fit (nothing in the tree runs this pair; mc_ri R7 last did on 2026-08-24).
* Placebo fixture + rseed(1342) so the observed statistic sits INSIDE the sampled
* distribution (2026-09-13: tails 26/5, p .387); with a planted effect it was the extreme
* draw and the p check only saw the floor ====
mkperm3 8 4 0
qui gen double w = 0.5 + runiform()
set seed 134
cap noi wooldid y i t ttre [aw=w], att cluster(i) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(1342) dots(0)
local rc = _rc
xrc Q34a 0 `rc'
if `rc' == 0 {
  mata: st_numscalar("__bad", wdp3_recon(0, 95))
  if __bad != 0 {
    dis as error "Q34b FAIL: " __bad " element(s) whose p/interval is not the convention applied to the posted draws under [aw=w]"
    global pf = $pf + 1
  }
  else dis as result "Q34b PASS: weighted permutation p/intervals reconstruct bitwise"
  local q34m = e(perm_completed)
  if `q34m' != 30 {
    dis as error "Q34c FAIL: perm_completed=`q34m' (weighted replications failed)"
    global pf = $pf + 1
  }
  * the observed att is interior: its smaller tail among the 30 draws holds at least two
  mata: A_ = st_matrix("e(permdist_b)"); k_ = selectindex(st_matrixcolstripe("e(permdist_b)")[., 2] :== "main_att"); v_ = select(A_[., k_], A_[., k_] :< .); c_ = st_matrix("e(perm_comparator_b)")[1, k_]; st_numscalar("__q34tail", min((sum(v_ :<= c_), sum(v_ :>= c_))))
  if __q34tail < 2 {
    dis as error "Q34d FAIL: smaller tail = " __q34tail " of 30 - the observed statistic is the extreme draw, so the p identity only saw the floor; change the fixture or rseed"
    global pf = $pf + 1
  }
  else dis "Q34d info: smaller tail " __q34tail " of 30 (p=" e(perm_p_att) ")"
}

* ================= Q35: poisson under permute() - the same identity on ppmlhdfe replications
* (att only; the fixture is interior as it stands: 2026-09-13 tails 24/8, p .581) ====
mkperm3 8 4 3
qui gen ycount = round(exp(0.3*(t >= ttre & ttre < .) + 0.05*i))
set seed 135
cap noi wooldid ycount i t ttre, att cluster(i) poisson permute(shuffle) permreps(30) rseed(1350) dots(0)
local rc = _rc
xrc Q35a 0 `rc'
if `rc' == 0 {
  mata: st_numscalar("__bad", wdp3_recon(0, 95))
  if __bad != 0 {
    dis as error "Q35b FAIL: " __bad " element(s) whose p/interval is not the convention applied to the posted draws under poisson"
    global pf = $pf + 1
  }
  else dis as result "Q35b PASS: poisson permutation p/intervals reconstruct bitwise"
  if e(perm_completed) != 30 {
    dis as error "Q35c FAIL: perm_completed=" e(perm_completed) " (ppmlhdfe replications failed)"
    global pf = $pf + 1
  }
  mata: A_ = st_matrix("e(permdist_b)"); k_ = selectindex(st_matrixcolstripe("e(permdist_b)")[., 2] :== "main_att"); v_ = select(A_[., k_], A_[., k_] :< .); c_ = st_matrix("e(perm_comparator_b)")[1, k_]; st_numscalar("__q35tail", min((sum(v_ :<= c_), sum(v_ :>= c_))))
  if __q35tail < 2 {
    dis as error "Q35d FAIL: smaller tail = " __q35tail " of 30 - the observed statistic is the extreme draw, so the p identity only saw the floor; change the fixture or rseed"
    global pf = $pf + 1
  }
  else dis "Q35d info: smaller tail " __q35tail " of 30 (p=" e(perm_p_att) ")"
}

* ================= Q36: the hidden analytic twin e(perm_analytic_pwsgs) ===============
* The pairwise block snapshots e(pwsgs) BEFORE the permutation rewrites its inference
* columns and posts the snapshot as the hidden matrix e(perm_analytic_pwsgs) - the only
* place a user can still read the analytic contrast inference the displayed table replaced.
* Until 2026-09-13 nothing under test/ read it except the battery's own %21x dump of
* perm_sgjt_adv, and on that fixture every inference cell of it is missing, so a regression
* that stopped filling the live columns would have diffed nothing at all.
* This leg runs it where the columns ARE live. permanalyticse is required: the permute()
* default is fast mode (the shell sets fastnose itself unless permtstat or permanalyticse),
* and fast mode computes no analytic standard error, so the default run's twin is all
* missing by construction. What is asserted: the twin carries the ANALYTIC column names
* while the displayed table carries the permutation ones; the estimate and the analytic SE
* are shared bit for bit; column 3 is live in the twin and BLANK in the displayed table
* (contrasts are permuted as coefficient differences in both modes, so there is no permuted
* statistic to report); the twin's own columns are internally consistent as an analytic
* table (t = estimate/se, the interval symmetric about the estimate, its half width the
* t critical value at e(df_r) and e(level) times the SE); and the analytic p does NOT
* appear in the displayed p column, which is the leak the hidden twin exists to make
* checkable. The fixture is dta_binary, where both contrast rows carry a live analytic SE.
* Note what is NOT claimed: the subgroup is nested in cluster(i), and the BS-15 screen that
* would blank these rows on the analytic path does not run under permute(), so the raw
* analytic cells are what the twin carries - the leg asserts the snapshot identity and the
* analytic algebra, not the statistical quality of that standard error.
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) permute(shuffle) permreps(20) permanalyticse rseed(4242) dots(0)
local rc36 = _rc
xrc Q36a 0 `rc36'
if `rc36' == 0 {
  local bad 0
  cap confirm matrix e(perm_analytic_pwsgs)
  if _rc != 0 {
    dis as error "Q36b FAIL: e(perm_analytic_pwsgs) is not posted under permute + subgroups + permanalyticse"
    global pf = $pf + 1
  }
  else {
    tempname PA36 PW36
    matrix `PA36' = e(perm_analytic_pwsgs)
    matrix `PW36' = e(pwsgs)
    local ca36 : colnames `PA36'
    local cw36 : colnames `PW36'
    local ra36 : rowfullnames `PA36'
    local rw36 : rowfullnames `PW36'
    local df36 = e(df_r)
    local lv36 = e(level)
    local tc36 = invttail(`df36', (100 - `lv36')/200)
    dis "Q36 twin cols=[`ca36'] display cols=[`cw36'] rows=" rowsof(`PA36') " df_r=`df36' tcrit=" %12.9f `tc36' " analyticse=[`e(perm_analyticse)']"
    if rowsof(`PA36') != rowsof(`PW36') | colsof(`PA36') != 6 | "`ra36'" != "`rw36'" local ++bad
    if "`e(perm_analyticse)'" != "1" local ++bad
    * the twin keeps the ANALYTIC header, the display carries the permutation header
    if "`ca36'" != "estimate se t p lb_95ci ub_95ci" local ++bad
    if "`cw36'" != "estimate se_analytic t_perm p_perm lb_95ci ub_95ci" local ++bad
    local live36 0
    forval r = 1/`=rowsof(`PA36')' {
      * estimate and analytic SE are the same object in both tables
      forval c = 1/2 {
        if (`PA36'[`r',`c'] < . | `PW36'[`r',`c'] < .) & reldif(`PA36'[`r',`c'], `PW36'[`r',`c']) > 0 local ++bad
      }
      if `PA36'[`r',2] < . & `PA36'[`r',2] > 0 {
        local ++live36
        * the analytic statistic survives only in the twin; the display blanks it
        if `PA36'[`r',3] >= . local ++bad
        if `PW36'[`r',3] < . local ++bad
        * and the twin's columns are an analytic table, not leftovers
        if reldif(`PA36'[`r',3], `PA36'[`r',1]/`PA36'[`r',2]) > 1e-12 local ++bad
        if `PA36'[`r',5] >= . | `PA36'[`r',6] >= . local ++bad
        else {
          if reldif((`PA36'[`r',5] + `PA36'[`r',6])/2, `PA36'[`r',1]) > 1e-12 local ++bad
          if reldif((`PA36'[`r',6] - `PA36'[`r',5])/(2*`PA36'[`r',2]), `tc36') > 1e-6 local ++bad
        }
        * the leak the twin exists to make checkable: the analytic p is not the reported p
        if `PA36'[`r',4] < . & `PW36'[`r',4] < . & `PA36'[`r',4] == `PW36'[`r',4] local ++bad
      }
    }
    dis "Q36 live analytic rows=`live36' of " rowsof(`PA36') " bad=`bad'"
    if `live36' == 0 {
      dis as error "Q36c FAIL (VACUOUS): every analytic SE in the twin is missing on this fixture, so nothing was asserted - change the fixture"
      global pf = $pf + 1
    }
    else if `bad' {
      dis as error "Q36d FAIL: `bad' check(s) on the analytic pairwise twin"
      global pf = $pf + 1
    }
    else dis as result "Q36 PASS: e(perm_analytic_pwsgs) is the pre-permutation analytic snapshot of e(pwsgs)"
  }
}

* ================= Q37 (2026-09-21 ruling, built ed6a0da; gate added 2026-09-22): a permsave()
* write that fails after the loop is a caution, not a refusal. Forced with a read-only target:
* rc 0, results intact, the caution names the file, e(perm_save_rc) carries the code, the
* previous file is untouched, and a bare wooldid replays ===============================
mkperm3 8 4 3
cap erase perm3_q37.dta
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permsave(perm3_q37) rseed(37)
xrc Q37a 0 _rc
cap confirm file perm3_q37.dta
xrc Q37b 0 _rc
! chmod 444 perm3_q37.dta
mkperm3 8 4 3
qui log using perm3_q37.txt, text replace name(p3)
cap noi wooldid y i t ttre, att cluster(i) permute(shuffle) permreps(20) permsave(perm3_q37) replace rseed(38)
local rc = _rc
qui log close p3
xrc Q37c 0 `rc'
local q37src `e(perm_save_rc)'
local lg = fileread("perm3_q37.txt")
if strpos(`"`lg'"', "Caution: permsave(perm3_q37) could not produce the replication dataset") == 0 | "`q37src'" == "" | "`q37src'" == "0" {
  dis as error "Q37d FAIL: the failed write did not caution, or did not post e(perm_save_rc) (got `q37src')"
  global pf = $pf + 1
}
* Q37h the caution prints below the result tables, after the part(perm) footer, with the
* run's other end-of-run cautions rather than above the deferred display
local q37pf = strpos(`"`lg'"', "Randomization Inference (permute(")
local q37pc = strpos(`"`lg'"', "Caution: permsave(perm3_q37) could not produce the replication dataset")
if `q37pf' == 0 | `q37pc' < `q37pf' {
  dis as error "Q37h FAIL: the failed-write caution (at `q37pc') does not print below the permutation footer (at `q37pf')"
  global pf = $pf + 1
}
if e(perm_completed) != 20 {
  dis as error "Q37e FAIL: the results are not intact after the failed write"
  global pf = $pf + 1
}
cap noi wooldid
xrc Q37f 0 _rc
preserve
usegate3 perm3_q37 Q37s
local q37seed `_dta[wooldid_seed]'
qui count if _rep > 0
local q37n = r(N)
restore
if `q37n' != 20 | "`q37seed'" != "37" {
  dis as error "Q37g FAIL: the read-only file was altered (rows=`q37n' seed=`q37seed'; want 20 and 37)"
  global pf = $pf + 1
}
! chmod 644 perm3_q37.dta
cap erase perm3_q37.dta
cap erase perm3_q37.txt

if $pf == 0 dis as result "PERMUTE3_ALL_PASS"
else dis as error "PERMUTE3_FAILS=$pf"
