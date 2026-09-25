* test_newcaut.do - the four cautions the author ruled on 2026-09-25 (econ review, section 8),
* each on the design its refuter measured and on a healthy design that must stay silent.
*   N1-N5 CT-2: under lattice(), a treated cohort-period-by-lattice-bin piece with a single
*         nonzero dose borrows another bin's slope. Refuter designs p3b (the bin-3 piece is
*         single-valued) and q1 (the base-bin piece is), each with the count; a healthy
*         lattice stays silent; safety(off) silences it with the rest of the collinearity
*         screens; verbose prints the reworded dropped-term line and not the old one.
*   N6-N7 CT-3: under contreatwithin, periods with no untreated reference rows and a dropped
*         treated cell dummy are named (refuter design h); with a never-treated group it
*         stays silent.
*   N8-N10 ES-3: the analytic simultaneousbands band with cluster() nesting i cautions once
*         the band is built; robust and a cluster() finer than i stay silent; the bootstrap
*         band stays silent.
*   N11-N13 ES-6: a rank-deficient event-study covariance behind the adversarial test
*         cautions with its rank, and e(adversarial_df_att) is the rank minus one; a full-rank
*         covariance stays silent with df = K - 1; under bootstrap() no caution and no df.
* Needles are searched in a named text log after it is closed, so the echoed commands never
* hold them. Run from test/: "<stata>" -e do test_newcaut.do -> expect NEWCAUT_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
local DATA "data"
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop nclog
program define nclog, rclass
  * nclog <logfile> <needle>: 1 if the needle appears on one line of the closed log
  args lf needle
  tempname fh
  local found 0
  file open `fh' using "`lf'", read text
  file read `fh' line
  while r(eof) == 0 {
    if strpos(`"`macval(line)'"', `"`needle'"') != 0 local found = `found' + 1
    file read `fh' line
  }
  file close `fh'
  return scalar found = `found'
end

cap program drop mkct2
program define mkct2
  args mode
  clear
  set seed 4243
  set obs 300
  gen unit = _n
  gen i = ceil(unit/60)
  gen ttre = cond(i<=3, 4, .)
  gen double d = 3*runiform()
  if "`mode'" == "p3b" replace d = 2.5 if i==1 & d>=2
  if "`mode'" == "q1" replace d = 0.5 if i==1 & d<1
  expand 6
  bysort unit: gen t = _n
  gen tr = (t>=ttre & ttre<.)
  gen lat = 1 + (d>=1) + (d>=2)
  gen double eff = cond(lat==1, d, cond(lat==2, 3*d, 5*d))
  gen double y = 0.5*d + t + tr*(1 + eff) + 0.01*rnormal()
end

* the CT-2 caution's opening, split so no echoed line of this file can carry it whole
local ct2a "treated cohort-period-by-lattice-bin"
local ct2b "piece(s) have no within-piece variation in d"

* ---- N1-N4: CT-2 on the refuter designs, a healthy lattice, and under safety(off)
foreach leg in "N1 p3b 3" "N2 q1 3" "N3 healthy 0" "N4 p3b 0 safety(off)" {
  local tag : word 1 of `leg'
  local mode : word 2 of `leg'
  local want : word 3 of `leg'
  local extra : word 4 of `leg'
  mkct2 `mode'
  qui log using newcaut_ct2.txt, text replace name(nc)
  cap noi wooldid y i t ttre, att contreat(d) contreatcontrols(basic) robust lattice(lat) `extra'
  local rc = _rc
  qui log close nc
  nclog newcaut_ct2.txt "`ct2a' `ct2b'"
  local any = r(found)
  nclog newcaut_ct2.txt "Caution: `want' `ct2a' `ct2b'"
  local hit = r(found)
  cap erase newcaut_ct2.txt
  local ok = (`rc' == 0)
  if `want' > 0 local ok = `ok' & `hit' == 1 & `any' == 1
  else local ok = `ok' & `any' == 0
  if `ok' di as result "`tag' PASS: CT-2 [`mode' `extra'] rc=`rc' caution lines=`any' (want count `want')"
  else {
    di as error "`tag' FAIL: CT-2 [`mode' `extra'] rc=`rc' caution lines=`any' with count `want'=`hit'"
    local ++fails
  }
}

* ---- N5: verbose prints the reworded dropped-term line under a lattice, never the old one
mkct2 p3b
qui log using newcaut_v.txt, text replace name(nc)
cap noi wooldid y i t ttre, att contreat(d) contreatcontrols(basic) robust lattice(lat) verbose
local rc = _rc
qui log close nc
local o1 "If using a lattice,"
local o2 "this is likely innocuous"
nclog newcaut_v.txt "`o1' `o2'"
local old = r(found)
local v1 "Under a lattice a dropped piece takes the slope"
local v2 "of another bin in its cell"
nclog newcaut_v.txt "`v1' `v2'"
local new = r(found)
nclog newcaut_v.txt "`ct2a' `ct2b'"
local c2 = r(found)
cap erase newcaut_v.txt
if `rc' == 0 & `old' == 0 & `new' >= 1 & `c2' == 1 di as result "N5 PASS: verbose prints the reworded lattice line (`new'x), not the old one, and CT-2 once"
else {
  di as error "N5 FAIL: rc=`rc' old line=`old' new line=`new' CT-2=`c2'"
  local ++fails
}

* ---- N6-N7: CT-3 on the refuter design h (onsets 3 and 5, no never-treated) and with a
* never-treated group
cap program drop mkct3
program define mkct3
  args mode
  clear
  set seed 7
  set obs 240
  gen unit = _n
  gen i = ceil(unit/60)
  gen ttre = cond(i<=2, 3, 5)
  gen double d = 1 + 2*runiform()
  expand 6
  bysort unit: gen t = _n
  gen tr = (t>=ttre)
  gen double te = i + d
  gen double y = 0.5*d + t + tr*te + 0.01*rnormal()
  if "`mode'" == "healthy" {
    replace i = 5 if unit > 200
    replace ttre = . if unit > 200
    replace tr = 0 if unit > 200
    replace y = 0.5*d + t + 0.01*rnormal() if unit > 200
  }
end
local ct3a "Caution: under contreatwithin,"
local ct3b "period(s) have no untreated reference rows"
foreach leg in "N6 h" "N7 healthy" {
  local tag : word 1 of `leg'
  local mode : word 2 of `leg'
  mkct3 `mode'
  qui log using newcaut_ct3.txt, text replace name(nc)
  cap noi wooldid y i t ttre, att contreat(d) contreatcontrols(basic) contreatwithin robust espre(1) espost(3)
  local rc = _rc
  qui log close nc
  nclog newcaut_ct3.txt "`ct3a' 3 `ct3b' (4, 5, 6)"
  local hit = r(found)
  nclog newcaut_ct3.txt "`ct3b'"
  local any = r(found)
  cap erase newcaut_ct3.txt
  if "`mode'" == "h" local ok = `rc' == 0 & `hit' == 1 & `any' == 1
  else local ok = `rc' == 0 & `any' == 0
  if `ok' di as result "`tag' PASS: CT-3 [`mode'] rc=`rc' caution lines=`any' naming periods 4, 5, 6=`hit'"
  else {
    di as error "`tag' FAIL: CT-3 [`mode'] rc=`rc' caution lines=`any' naming periods 4, 5, 6=`hit'"
    local ++fails
  }
}

* ---- N8-N10: ES-3 on dta_binary (i is the unit)
local es3a "simultaneousbands is built from the analytic covariance"
local es3b "with cluster() at the level of i or above"
use "`DATA'/dta_binary", clear
gen coh = cond(ttre == ., 1, ttre)
tempfile bin
save `bin'
foreach leg in "N8a i cluster(i) 1" "N8b i customvce(cluster_i) 1" "N8c coh cluster(coh) 1" "N9a i robust 0" "N9b coh cluster(i) 0" "N10 i cluster(i)_bootstrap(standard)_reps(19)_rseed(5)_dots(0) 0" {
  local tag : word 1 of `leg'
  local iv : word 2 of `leg'
  local vce = subinstr("`: word 3 of `leg''", "_", " ", .)
  local want : word 4 of `leg'
  use `bin', clear
  qui log using newcaut_es3.txt, text replace name(nc)
  cap noi wooldid y `iv' t ttre, att esprelength(2) espostlength(2) simultaneousbands `vce'
  local rc = _rc
  qui log close nc
  cap confirm matrix e(esband_att)
  local band = _rc == 0
  nclog newcaut_es3.txt "`es3a' `es3b'"
  local hit = r(found)
  cap erase newcaut_es3.txt
  if `rc' == 0 & `band' == 1 & `hit' == `want' di as result "`tag' PASS: ES-3 [i=`iv' `vce'] band=`band' caution=`hit' (want `want')"
  else {
    di as error "`tag' FAIL: ES-3 [i=`iv' `vce'] rc=`rc' band=`band' caution=`hit' (want `want')"
    local ++fails
  }
}

* ---- N11-N13: ES-6 on dta_binary, cluster() on four clusters (rank-deficient), robust
* (full rank), and the bootstrap run (no caution, no df)
local es6a "the event-study covariance for the adversarial test is rank-deficient"
foreach leg in "N11 cluster(c4) 1" "N12 robust 0" "N13 cluster(c4)_bootstrap(standard)_boottstat_reps(19)_rseed(5)_dots(0) 0" {
  local tag : word 1 of `leg'
  local vce = subinstr("`: word 2 of `leg''", "_", " ", .)
  local want : word 3 of `leg'
  use `bin', clear
  gen c4 = mod(i, 4)
  qui log using newcaut_es6.txt, text replace name(nc)
  cap noi wooldid y i t ttre, att esfixedbaseperiod esprelength(4) espostlength(3) adversarial `vce'
  local rc = _rc
  qui log close nc
  nclog newcaut_es6.txt "`es6a'"
  local hit = r(found)
  cap erase newcaut_es6.txt
  local k1 .
  local K .
  cap confirm matrix e(V_es_att)
  if _rc == 0 {
    mata: V = st_matrix("e(V_es_att)"); V = (V + V')/2; L = symeigenvalues(V); st_local("k1", strofreal(sum(L :> max(abs(L))*1e-12))); st_local("K", strofreal(rows(V)))
  }
  local df = e(adversarial_df_att)
  local isboot = strpos("`vce'", "bootstrap") > 0
  if `isboot' local ok = `rc' == 0 & `hit' == 0 & missing(`df')
  else {
    local ok = `rc' == 0 & `hit' == `want' & `df' == `k1' - 1 & !missing(e(adversarial_pval_att))
    if `want' local ok = `ok' & `k1' < `K'
    else local ok = `ok' & `k1' == `K'
  }
  if `ok' di as result "`tag' PASS: ES-6 [`vce'] rank `k1' of `K' caution=`hit' df=`df'"
  else {
    di as error "`tag' FAIL: ES-6 [`vce'] rc=`rc' rank `k1' of `K' caution=`hit' (want `want') df=`df'"
    local ++fails
  }
}

if `fails' == 0 di as result "NEWCAUT_ALL_PASS"
else di as error "NEWCAUT_FAILS=`fails'"
exit, clear
