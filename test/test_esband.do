* test_esband.do - PS-45 (author ruling 2026-09-07): simultaneous (sup-t / uniform)
* confidence bands over an event study's horizons, behind simultaneousbands.
* Run from test/: stata-mp -e do test_esband.do  -> expect ESBAND_ALL_PASS
* B1 nothing is posted without the option (it is defaulted off, and no shipped figure or
*    table changes); B2 the analytic band is b +/- crit*se on the family's own covariance,
*    with crit at or above the pointwise normal quantile - a simultaneous band is never
*    narrower than a pointwise interval; B3 it is reproducible and leaves the caller's
*    random-number stream where it found it; B4 the bootstrap engine builds the band from
*    the replications; B5 permute() refuses the combination; B6 a run that cannot build one
*    says so instead of staying silent; B7 subgroups get their own bands; B8 the plot layer
*    draws without error.
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
    zc = invnormal(1 - (100 - st_numscalar("e(esband_level)"))/200)
    bad = 0
    asym = 0
    narrow = 0
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
      if (cc < zc) narrow = narrow + 1
    }
    st_numscalar("b2bad", bad)
    st_numscalar("b2asym", asym)
    st_numscalar("b2narrow", narrow)
  end
  di "B2 engine=`b2e' families=[`b2f'] level=`b2l' draws=`b2d' crit=" %8.5f `b2c' " z=" %8.5f invnormal(1 - (100 - `b2l')/200) " rows=[`b2rn'] cols=[`b2cn']"
  if "`b2e'" != "analytic-simulated" | "`b2f'" == "" | `b2l' != 95 | `b2d' != 10000 ///
     | "`b2rn'" != "lb ub relyr" | b2bad != 0 | b2asym != 0 | b2narrow != 0 | `b2c' <= invnormal(.975) {
    di as error "B2 FAIL: engine=`b2e' rows=[`b2rn'] mismatched=" b2bad " asymmetric=" b2asym " narrower-than-pointwise=" b2narrow " crit=`b2c'"
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
if "`b9eng'" != "bootstrap-studentized" | "`b9eng2'" != "bootstrap" | `b9seen' < 2 | `b9narrow' != 0 {
  di as error "B9 FAIL: boottstat engine=[`b9eng'] (want bootstrap-studentized), plain engine=[`b9eng2'] (want bootstrap), horizons checked=`b9seen', narrower-than-pointwise=`b9narrow'"
  local ++fails
}
else di as result "B9 PASS: boottstat gives a percentile-t band (`b9seen' horizons, none narrower than its pointwise interval); without boottstat the draw-SD construction stands"

if `fails' == 0 di as result "ESBAND_ALL_PASS"
else di as error "ESBAND_FAILS=`fails'"
exit, clear
