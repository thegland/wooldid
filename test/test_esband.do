* test_esband.do - PS-45 (author ruling 2026-09-07): simultaneous (sup-t / uniform)
* confidence bands over an event study's horizons, behind simultaneousbands.
* Run from test/: stata-mp -e do test_esband.do  -> expect ESBAND_ALL_PASS
* B1 nothing is posted without the option (it is defaulted off, and no shipped figure or
*    table changes); B2 the analytic band is b +/- crit*se on the family's own covariance,
*    with crit the sup-t quantile of the multivariate t with df_r; B3 it is reproducible
*    and leaves the caller's
*    random-number stream where it found it; B4 the bootstrap engine builds the band from
*    the replications; B5 permute() refuses the combination; B6 a run that cannot build one
*    says so instead of staying silent; B7 subgroups get their own bands; B8 the plot layer
*    draws without error; B11 with few clusters (df_r = 6) the critical value follows the
*    multivariate t with df_r and not the normal law; B12 the critical value is the
*    simulated quantile itself, with no floor at the pointwise t quantile; B13 the dashed
*    band on the figure pinches to 0 at the reference period, as the shading does.
clear all
set more off
set graphics off
local REPO ..
do "`REPO'/wooldidinstall/wooldid.ado"
local fails 0

cap prog drop mkb
prog def mkb
  clear
  set obs 24
  gen unit = _n
  gen sg = mod(unit, 2)
  gen ttre = cond(unit <= 9, 5, cond(unit <= 17, 7, .))
  expand 10
  bysort unit: gen t = _n
  set seed 71
  gen y = 1.4*(t >= ttre & ttre < .) + 0.4*(t >= ttre & ttre < .)*sg + 0.1*unit + 0.05*t + rnormal()
end

local SPEC att cluster(unit) esprelength(2) espostlength(2)

* ------------------------------------------------------------------ B1: off by default
mkb
qui wooldid y unit t ttre, `SPEC'
local b1n 0
foreach m in esband_att esband_att_it {
  cap confirm matrix e(`m')
  if _rc == 0 local ++b1n
}
local b1e = "`e(esband_engine)'"
di "B1 matrices=`b1n' engine=[`b1e']"
if `b1n' != 0 | "`b1e'" != "" {
  di as error "B1 FAIL: a band was posted without simultaneousbands (matrices=`b1n', engine=[`b1e'])"
  local ++fails
}

* ------------------------------------------------- B2: the analytic band and its geometry
mkb
qui wooldid y unit t ttre, `SPEC' simultaneousbands
local b2e = "`e(esband_engine)'"
local b2f = "`e(esband_families)'"
local b2l = e(esband_level)
local b2d = e(esband_ndraw)
cap confirm matrix e(esband_att)
local b2has = _rc == 0
if !`b2has' {
  di as error "B2 FAIL: no e(esband_att) on an analytic run that asked for one"
  local ++fails
}
else {
  matrix B2 = e(esband_att)
  local b2rn : rownames B2
  local b2cn : colnames B2
  local b2c = e(esband_crit_att)
  mata:
    B = st_matrix("B2")
    V = st_matrix("e(V_es_att)")
    b = st_matrix("e(b_es_att)")
    vs = st_matrixcolstripe("e(V_es_att)")
    bs = st_matrixcolstripe("B2")
    cc = st_numscalar("e(esband_crit_att)")
    bad = 0
    asym = 0
    for (j = 1; j <= cols(B); j++) {
      a = 0
      for (r = 1; r <= rows(vs); r++) if (vs[r, 2] == bs[j, 2]) a = r
      if (a == 0) {
        bad = bad + 1
        continue
      }
      se = sqrt(V[a, a])
      if (reldif(B[1, j], b[a] - cc*se) > 1e-12) bad = bad + 1
      if (reldif(B[2, j], b[a] + cc*se) > 1e-12) bad = bad + 1
      if (reldif((B[1, j] + B[2, j])/2, b[a]) > 1e-10) asym = asym + 1
    }
    st_numscalar("b2bad", bad)
    st_numscalar("b2asym", asym)
  end
  * independent sup-t: 200k draws from the family correlation, each draw's max|z| divided by
  * sqrt(chi2(df_r)/df_r) - the multivariate t with the df of the pointwise intervals - and its
  * level() quantile taken as it is, as the engine does. The engine's 10k-draw quantile has MC
  * sd 0.0215 on this fixture (df_r = 23; 200 seeds, 2026-09-23), so 0.065 is a ~3-sd band. On
  * that measurement the kernel at 200k draws gives 2.6541, the hand 2.6513, and the fixed-seed
  * engine 2.6099 (2 sd low); the normal-law value is 2.4694, outside the band, so this check
  * also fails an engine that draws from the normal law
  local b2df = e(df_r)
  mata: V = st_matrix("e(V_es_att)"); sd = sqrt(diagonal(V)); R = (V :/ (sd*sd')); R = (R + R')/2; rseed(20260912); df = st_numscalar("e(df_r)"); Z = rnormal(200000, rows(R), 0, 1) * cholesky(R)'; ms = sort(rowmax(abs(Z)) :/ sqrt(rchi2(200000, 1, df) :/ df), 1); st_numscalar("b2hand", ms[ceil((1 - (100 - st_numscalar("e(esband_level)"))/100)*rows(ms))])
  di "B2 hand sup-t crit=" %8.5f b2hand " Bonferroni bound=" %8.5f invttail(`b2df', (100 - `b2l')/(200*colsof(B2))) " df_r=`b2df'"
  di "B2 engine=`b2e' families=[`b2f'] level=`b2l' draws=`b2d' crit=" %8.5f `b2c' " t=" %8.5f invttail(`b2df', (100 - `b2l')/200) " rows=[`b2rn'] cols=[`b2cn']"
  if "`b2e'" != "analytic-simulated" | "`b2f'" == "" | `b2l' != 95 | `b2d' != 10000 | !(`b2df' > 0 & `b2df' < .) ///
     | "`b2rn'" != "lb ub relyr" | b2bad != 0 | b2asym != 0 | `b2c' > invttail(`b2df', (100 - `b2l')/(200*colsof(B2))) | abs(`b2c' - b2hand) > 0.065 {
    di as error "B2 FAIL: engine=`b2e' rows=[`b2rn'] mismatched=" b2bad " asymmetric=" b2asym " crit=`b2c' hand=" b2hand
    local ++fails
  }
}

* -------------------------------------- B3: reproducible, and the RNG stream is restored
mkb
set seed 999
qui wooldid y unit t ttre, `SPEC'
local b3ctl = runiform()
mkb
set seed 999
qui wooldid y unit t ttre, `SPEC' simultaneousbands
local b3after = runiform()
matrix B3A = e(esband_att)
mkb
qui wooldid y unit t ttre, `SPEC' simultaneousbands
matrix B3B = e(esband_att)
di "B3 stream ctl=" %12.10f `b3ctl' " after band=" %12.10f `b3after' "  band reldif=" mreldif(B3A, B3B)
if `b3ctl' != `b3after' {
  di as error "B3 FAIL: the band's simulation moved the caller's random-number stream"
  local ++fails
}
if mreldif(B3A, B3B) != 0 {
  di as error "B3 FAIL: two identical runs produced different bands (fixed internal seed expected)"
  local ++fails
}

* ------------------------------------------------------------ B4: the bootstrap engine
mkb
qui wooldid y unit t ttre, `SPEC' bootstrap(standard) reps(60) rseed(4) simultaneousbands
local b4e = "`e(esband_engine)'"
local b4n = e(esband_n_att)
local b4c = e(esband_crit_att)
local b4done = e(boot_completed)
cap confirm matrix e(esband_att)
local b4has = _rc == 0
local b4bad 0
if `b4has' {
  mata:
    B = st_matrix("e(esband_att)")
    bo = st_matrix("e(b_es_att)")
    vs = st_matrixcolstripe("e(V_es_att)")
    bs = st_matrixcolstripe("e(esband_att)")
    bb = 0
    for (j = 1; j <= cols(B); j++) {
      a = 0
      for (r = 1; r <= rows(vs); r++) if (vs[r, 2] == bs[j, 2]) a = r
      if (a == 0) {
        bb = bb + 1
        continue
      }
      if (reldif((B[1, j] + B[2, j])/2, bo[a]) > 1e-10) bb = bb + 1
      if (B[2, j] <= B[1, j]) bb = bb + 1
    }
    st_numscalar("b4m", bb)
  end
  local b4bad = b4m
  * the critical value itself: max_h |b_rh - b_h| / sd_h over the listwise-complete draws,
  * at the band level on the averaging order-statistic rule (2026-09-12 smoke: exact on HEAD)
  mata: BD = st_matrix("e(bootdist)"); nm = st_matrixcolstripe("e(bootdist)")[., 2]; idx = selectindex(strpos(nm, "ES_att_") :== 1); A = select(BD[., idx], rowmissing(BD[., idx]) :== 0); bn = st_matrixcolstripe("e(b)")[., 2]; b0 = J(1, rows(idx), .); for (j = 1; j <= rows(idx); j++) b0[j] = st_matrix("e(b)")[1, selectindex(bn :== nm[idx[j]])]; sd = sqrt(diagonal(variance(A)))'; v = sort(rowmax(abs(A :- b0) :/ sd), 1); n = rows(v); x = n*(st_numscalar("e(esband_level)")/100); st_numscalar("b4hand", (x == floor(x) & x >= 1 ? (v[x] + v[min((x+1, n))])/2 : v[max((min((ceil(x), n)), 1))])); st_numscalar("b4n_hand", n)
  di "B4 hand crit=" %21x b4hand " engine=" %21x `b4c' " n_hand=" b4n_hand
  if reldif(b4hand, `b4c') > 1e-12 | b4n_hand != `b4n' local ++b4bad
}
di "B4 engine=`b4e' n=`b4n'/`b4done' crit=" %8.5f `b4c' " problems=`b4bad'"
if !`b4has' | "`b4e'" != "bootstrap" | `b4n' > `b4done' | `b4n' < 2 | `b4c' <= 0 | `b4bad' != 0 {
  di as error "B4 FAIL: has=`b4has' engine=`b4e' n=`b4n' completed=`b4done' crit=`b4c' problems=`b4bad'"
  local ++fails
}

* ------------------------------------------------------- B5: permute() refuses the option
mkb
cap noi wooldid y unit t ttre, `SPEC' permute(shuffle) permreps(20) rseed(5) simultaneousbands
local b5rc = _rc
di "B5 rc=`b5rc'"
if `b5rc' != 198 {
  di as error "B5 FAIL: simultaneousbands under permute() must be refused with r(198), got rc=`b5rc'"
  local ++fails
}

* --------------------------------------------- B6: a request that cannot be met says so
mkb
qui log using esband_b6.txt, text replace name(eb)
cap noi wooldid y unit t ttre, `SPEC' fastnose simultaneousbands
local b6rc = _rc
qui log close eb
local b6lg = fileread("esband_b6.txt")
local b6lg = subinstr(`"`b6lg'"', char(13), "", .)
local b6lg = subinstr(`"`b6lg'"', char(10) + "> ", "", .)
cap confirm matrix e(esband_att)
local b6has = _rc == 0
local b6c = strpos(`"`b6lg'"', "simultaneousbands was requested but no simultaneous band could be built") > 0
cap erase esband_b6.txt
di "B6 rc=`b6rc' matrix=`b6has' caution=`b6c'"
if `b6rc' != 0 | `b6has' != 0 | `b6c' != 1 {
  di as error "B6 FAIL: fastnose rc=`b6rc' matrix posted=`b6has' caution printed=`b6c'"
  local ++fails
}

* ----------------------------------------------------------- B7: per-subgroup bands
mkb
qui wooldid y unit t ttre, `SPEC' subgroups(sg) simultaneousbands
local b7f = "`e(esband_families)'"
local b7n 0
foreach s in 0 1 {
  cap confirm matrix e(esband_att_s`s')
  if _rc == 0 local ++b7n
}
di "B7 families=[`b7f'] subgroup matrices=`b7n'"
if `b7n' != 2 | strpos("`b7f'", "att_s0") == 0 | strpos("`b7f'", "att_s1") == 0 {
  di as error "B7 FAIL: subgroup bands=`b7n' families=[`b7f']"
  local ++fails
}

* ------------------------------------------------------ B8: the plot layer runs cleanly
mkb
cap noi wooldid y unit t ttre, `SPEC' subgroups(sg) simultaneousbands makeplots
local b8rc = _rc
di "B8 makeplots rc=`b8rc'"
if `b8rc' != 0 {
  di as error "B8 FAIL: the simultaneous-band plot layer errored (rc=`b8rc')"
  local ++fails
}

* ---- B9 PS-45 (author ruling 2026-09-10): under boottstat the band is a genuine
* percentile-t - each replication's own standard error in the denominator, the observed
* analytic ones in the endpoints - so it sits on the SAME scale as the pointwise intervals
* printed beside it. The reason the ruling went this way: on the old coefficient-draw-SD
* construction the two scales differ, so a SIMULTANEOUS band could come out NARROWER than
* the pointwise interval it accompanies, which is visibly wrong. Two contracts: the engine
* marker says which construction ran, and the band is at least as wide as the pointwise
* interval at every horizon. ----
mkb
qui wooldid y unit t ttre, `SPEC' bootstrap(standard) reps(60) boottstat rseed(91) dots(0) simultaneousbands
local b9eng "`e(esband_engine)'"
tempname B9B B9M
matrix `B9B' = e(esband_att)
matrix `B9M' = e(wdidmainresults)
* the studentized critical value itself (audit 2026-09-13; kernel = wooldid_bootband): over the
* replications complete in BOTH e(bootdist) and e(bootdist_se) at every horizon, max_h |b_rh - b_h| / se_rh,
* at the band level on the averaging order-statistic rule; the endpoints are b_h -/+ crit * the
* OBSERVED analytic se (hidden e(boot_ase_obs)) - probe 2026-09-13: reldif 0 on HEAD
mata:
  BD = st_matrix("e(bootdist)"); SD = st_matrix("e(bootdist_se)")
  nm = st_matrixcolstripe("e(bootdist)")[., 2]
  idx = selectindex(strpos(nm, "ES_att_") :== 1)
  bn = st_matrixcolstripe("e(b)")[., 2]
  b0 = J(1, rows(idx), .)
  for (j = 1; j <= rows(idx); j++) b0[j] = st_matrix("e(b)")[1, selectindex(bn :== nm[idx[j]])]
  A = BD[., idx]; S = SD[., idx]
  ok = (rowmissing(A) :== 0) :& (rowmissing(S) :== 0) :& (rowmin(S) :> 0)
  A = select(A, ok); S = select(S, ok)
  v = sort(rowmax(abs(A :- b0) :/ S), 1)
  n = rows(v)
  x = n*(st_numscalar("e(esband_level)")/100)
  c = (x == floor(x) & x >= 1 ? (v[x] + v[min((x+1, n))])/2 : v[max((min((ceil(x), n)), 1))])
  // endpoints against the observed analytic se, matched by name
  AO = st_matrix("e(boot_ase_obs)"); an = st_matrixcolstripe("e(boot_ase_obs)")[., 2]
  EB_ = st_matrix("e(esband_att)"); es = st_matrixcolstripe("e(esband_att)")[., 2]
  ebad = 0
  for (j = 1; j <= cols(EB_); j++) {
    ja = selectindex(an :== ("ES_att_" + es[j]))
    jb = selectindex(nm[idx] :== ("ES_att_" + es[j]))
    // Mata rejects `if (c) { x++; continue }` on ONE line ("'}' found where nothing
    // expected"); the block has to span several lines.
    if (rows(ja) != 1 | rows(jb) != 1) {
      ebad++
      continue
    }
    if (reldif(EB_[1, j], b0[jb] - c*AO[1, ja]) > 1e-12 | reldif(EB_[2, j], b0[jb] + c*AO[1, ja]) > 1e-12) ebad++
  }
  st_numscalar("b9hand", c); st_numscalar("b9n_hand", n); st_numscalar("b9ebad", ebad)
end
local b9n = e(esband_n_att)
local b9c = e(esband_crit_att)
di "B9 hand crit=" %21x b9hand " engine=" %21x `b9c' " n=`b9n'/" b9n_hand " endpoint mismatches=" b9ebad
local b9bad 0
if reldif(b9hand, `b9c') > 1e-12 | b9n_hand != `b9n' | b9ebad != 0 local b9bad 1
local b9rn : rowfullnames `B9M'
local b9narrow 0
local b9seen 0
forval c = 1/`=colsof(`B9B')' {
  local b9ry = `B9B'[3, `c']
  forval r = 1/`=rowsof(`B9M')' {
    local b9nm : word `r' of `b9rn'
    if strpos("`b9nm'", "ES_att") == 1 & `B9M'[`r', 7] == `b9ry' {
      local b9pw = (`B9M'[`r', 6] - `B9M'[`r', 5])/2
      local b9bw = (`B9B'[2, `c'] - `B9B'[1, `c'])/2
      local ++b9seen
      if `b9bw' < `b9pw' - 1e-12 local ++b9narrow
    }
  }
}
mkb
qui wooldid y unit t ttre, `SPEC' bootstrap(standard) reps(60) rseed(91) dots(0) simultaneousbands
local b9eng2 "`e(esband_engine)'"
if "`b9eng'" != "bootstrap-studentized" | "`b9eng2'" != "bootstrap" | `b9seen' < 2 | `b9narrow' != 0 | `b9bad' != 0 {
  di as error "B9 FAIL: boottstat engine=[`b9eng'] (want bootstrap-studentized), plain engine=[`b9eng2'] (want bootstrap), horizons checked=`b9seen', narrower-than-pointwise=`b9narrow', crit/endpoint reconstruction failed=`b9bad'"
  local ++fails
}
else di as result "B9 PASS: boottstat gives a percentile-t band (`b9seen' horizons, none narrower than its pointwise interval); without boottstat the draw-SD construction stands"

* ------------------------------------------------- B10: D-9, per-family seed (2026-09-11)
* the analytic band's critical value must not depend on which OTHER families this run
* happens to band. The seed used to be set once before the family loop, so each family
* consumed successive draws and adding an estimand or subgroups() moved an unrelated
* family's band in its low-order digits. Measured before the fix on this fixture:
* crit_att_it was 2.474072413106145 with att banded first and 2.492870068598342 when
* att_it drew first.
mkb
qui wooldid y unit t ttre, att att_it cluster(unit) esprelength(2) espostlength(2) simultaneousbands
local b10x = e(esband_crit_att_it)
mkb
qui wooldid y unit t ttre, att_it cluster(unit) esprelength(2) espostlength(2) simultaneousbands
local b10y = e(esband_crit_att_it)
* and still reproducible run to run
mkb
qui wooldid y unit t ttre, att att_it cluster(unit) esprelength(2) espostlength(2) simultaneousbands
local b10x2 = e(esband_crit_att_it)
if `b10x' >= . | `b10y' >= . | `b10x' != `b10y' | `b10x' != `b10x2' {
  di as error "B10 FAIL: crit_att_it must not depend on whether att drew before it - with att=`b10x', first=`b10y', repeat=`b10x2'"
  local ++fails
}
else di as result "B10 PASS: the band's critical value is keyed to the family, not to the draw order (crit_att_it=`b10x' either way, and reproducible)"

* ---- B11 (author ruling 2026-09-23): few clusters. The pointwise intervals are t(df_r)
* intervals, and the analytic band is drawn from the multivariate t with the same df_r. 7
* clusters, df_r = 6, horizon correlation .237. The engine's 10k-draw critical value has MC sd
* 0.036 on this fixture (200 seeds, 2026-09-23), so 0.11 is a ~3-sd band around an independent
* 200k-draw multivariate-t value (2.8986; the fixed-seed engine gives 2.8812). The normal law
* gives 2.2320 here (a 10k-draw normal-law engine gave 2.2146), far outside that band.
cap prog drop mkb7
prog def mkb7
  clear
  set obs 7
  gen unit = _n
  gen ttre = cond(unit <= 3, 5, cond(unit <= 5, 7, .))
  expand 10
  bysort unit: gen t = _n
  set seed 71
  gen y = 1.4*(t >= ttre & ttre < .) + 0.1*unit + 0.05*t + rnormal()
end
mkb7
qui wooldid y unit t ttre, att cluster(unit) esprelength(1) espostlength(1) simultaneousbands
local b11df = e(df_r)
local b11c = e(esband_crit_att)
local b11l = e(esband_level)
scalar b11hand = .
scalar b11norm = .
cap confirm matrix e(V_es_att)
if _rc == 0 & `b11df' > 0 & `b11df' < . {
  mata: V = st_matrix("e(V_es_att)"); sd = sqrt(diagonal(V)); R = (V :/ (sd*sd')); R = (R + R')/2; rseed(20260912); df = st_numscalar("e(df_r)"); Z = rnormal(200000, rows(R), 0, 1) * cholesky(R)'; mz = rowmax(abs(Z)); q = ceil((1 - (100 - st_numscalar("e(esband_level)"))/100)*200000); mt = sort(mz :/ sqrt(rchi2(200000, 1, df) :/ df), 1); mn = sort(mz, 1); st_numscalar("b11hand", mt[q]); st_numscalar("b11norm", mn[q])
}
di "B11 df_r=`b11df' crit=" %8.5f `b11c' " hand multivariate t=" %8.5f b11hand " normal law=" %8.5f b11norm " pointwise t quantile=" %8.5f invttail(`b11df', (100 - `b11l')/200)
if !(`b11df' > 0 & `b11df' < 15) | `b11c' >= . | b11hand >= . | abs(`b11c' - b11hand) > 0.11 | abs(`b11c' - b11norm) <= 0.11 {
  di as error "B11 FAIL: df_r=`b11df' crit=`b11c' hand multivariate t=" b11hand " normal law=" b11norm
  local ++fails
}
else di as result "B11 PASS: with df_r=`b11df' the analytic band's critical value follows the multivariate t with df_r, not the normal law"

* ---- B12 (author ruling 2026-09-23): the analytic critical value is the simulated level()
* quantile itself. It is not floored at the pointwise t quantile, so where simulation noise
* leaves it below that quantile the band is narrower than the pointwise interval, and the run
* shows it. The kernel is called directly on a rank-one family - two horizons with
* correlation 1 - where the sup over the horizons is a single |t_6|, so the 10k-draw quantile
* sits below invttail(6, .025) = 2.4469 for about half of all seeds; seed 7 gives 2.4131. The
* hand value replays the kernel's own draws on the same seed, so the two must agree, and a
* floor would have returned 2.4469 instead.
cap mata: mata drop b12hq()
mata:
real scalar b12hq(real matrix V0, real scalar lev, real scalar nd, real scalar df)
{
  real matrix V, U, L, Z
  real rowvector lam
  real colvector sd, M
  real scalar m, j, r
  m = rows(V0)
  V = V0
  sd = sqrt(diagonal(V))
  for (j = 1; j <= m; j++) for (r = 1; r <= m; r++) V[j, r] = V[j, r]/(sd[j]*sd[r])
  V = (V + V')/2
  U = .
  lam = .
  symeigensystem(V, U, lam)
  for (j = 1; j <= m; j++) if (lam[j] < 0 | lam[j] >= .) lam[j] = 0
  L = U * diag(sqrt(lam))
  Z = L * rnormal(m, nd, 0, 1)
  M = J(nd, 1, .)
  for (r = 1; r <= nd; r++) M[r] = max(abs(Z[., r]))
  if (df > 0 & df < .) M = M :/ sqrt(rchi2(nd, 1, df) :/ df)
  M = sort(M, 1)
  return(wooldid_bootqs(M, lev/100))
}
end
matrix B12V = (1, 1 \ 1, 1)
matrix colnames B12V = pre1 post0
matrix rownames B12V = pre1 post0
matrix B12B = (0.5, 1)
matrix colnames B12B = pre1 post0
local b12tq = invttail(6, .025)
local b12c
local b12m
set seed 7
cap mata: st_local("b12m", strofreal(wooldid_esbanda("B12V", "B12B", 95, 10000, "B12O", "b12c", 6)))
local b12rc = _rc
set seed 7
mata: st_numscalar("b12hand", b12hq(st_matrix("B12V"), 95, 10000, 6))
scalar b12crit = .
if "`b12c'" != "" scalar b12crit = `b12c'
di "B12 rc=`b12rc' horizons=`b12m' engine crit=" %21x b12crit " hand=" %21x b12hand " (" %8.5f b12hand ") pointwise t quantile=" %8.5f `b12tq'
if `b12rc' != 0 | "`b12m'" != "2" | b12hand >= `b12tq' | b12crit >= . | reldif(b12crit, b12hand) > 1e-12 {
  di as error "B12 FAIL: rc=`b12rc' horizons=`b12m' engine crit=" b12crit " simulated quantile=" b12hand " pointwise t quantile=`b12tq' - the critical value must be the simulated quantile, with no floor"
  local ++fails
}
else di as result "B12 PASS: the critical value is the simulated quantile (" %8.5f b12crit "), below the pointwise t quantile and not floored at it"

* ---- B13 esband-rline-across-reference (change audit 2026-09-23): the pointwise shading draws
* the reference period as a zero-width point, and the dashed band must pinch to it there too.
* The stored band holds the estimated horizons only, so the plot layer left the reference row
* missing and twoway's line joined the two neighbouring horizons straight across it: on an
* esfixedbaseperiod figure the dashed pair crossed the reference at full width while the
* shading closed to 0. The figure is read back from its saveplots() .gph, whose serset holds
* every plotted point (the band pair comes back as its last two variables, renamed). A row
* whose relyr is a band column must carry that column's endpoints; any other row is a
* reference row and must be 0/0. B13a the fixed base period (reference -1, between pre2 and
* post0); B13b the pooled reference (the * point at the left end)
cap prog drop b13chk
prog def b13chk, rclass
  args stem
  local bad 0
  local nref 0
  local nhit 0
  local ncol 0
  cap confirm matrix e(esband_att)
  if _rc == 0 {
    tempname B
    matrix `B' = e(esband_att)
    local ncol = colsof(`B')
    cap graph drop _all
    serset clear
    cap graph use "`stem'_ES_att.gph", name(b13g, replace)
    if _rc == 0 {
      serset set 0
      serset use, clear
      qui ds
      local sv `r(varlist)'
      if wordcount("`sv'") == 6 {
        local lo : word 5 of `sv'
        local hi : word 6 of `sv'
        forval r = 1/`=_N' {
          local m 0
          forval c = 1/`ncol' {
            if relyr[`r'] == `B'[3, `c'] {
              local m 1
              local ++nhit
              if reldif(`lo'[`r'], `B'[1, `c']) > 1e-6 | reldif(`hi'[`r'], `B'[2, `c']) > 1e-6 local ++bad
            }
          }
          if !`m' {
            local ++nref
            if `lo'[`r'] != 0 | `hi'[`r'] != 0 local ++bad
          }
        }
      }
      else local bad = `bad' + 100
    }
    else local bad = `bad' + 1000
    cap graph drop _all
  }
  return scalar bad = `bad'
  return scalar nref = `nref'
  return scalar nhit = `nhit'
  return scalar ncol = `ncol'
end
foreach b13v in a b {
  local b13opt = cond("`b13v'" == "a", "esfixedbaseperiod", "")
  local b13f : dir . files "esband_b13`b13v'_*"
  foreach f of local b13f {
    cap erase "`f'"
  }
  mkb
  cap noi wooldid y unit t ttre, `SPEC' `b13opt' simultaneousbands makeplots saveplots(esband_b13`b13v') replace
  local b13rc = _rc
  b13chk esband_b13`b13v'
  di "B13`b13v' [`b13opt'] rc=`b13rc' band columns=" r(ncol) " matched rows=" r(nhit) " reference rows=" r(nref) " bad=" r(bad)
  if `b13rc' != 0 | r(ncol) == 0 | r(nhit) != r(ncol) | r(nref) != 1 | r(bad) != 0 {
    di as error "B13`b13v' FAIL: the dashed band does not pinch to 0 at the reference period [`b13opt'] (rc=`b13rc' columns=" r(ncol) " matched=" r(nhit) " reference rows=" r(nref) " bad=" r(bad) ")"
    local ++fails
  }
  else di as result "B13`b13v' PASS: the dashed band carries the stored endpoints and pinches to 0 at the reference [`b13opt']"
  local b13f : dir . files "esband_b13`b13v'_*"
  foreach f of local b13f {
    cap erase "`f'"
  }
}

if `fails' == 0 di as result "ESBAND_ALL_PASS"
else di as error "ESBAND_FAILS=`fails'"
exit, clear
