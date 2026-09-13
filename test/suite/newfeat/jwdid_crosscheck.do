* jwdid_crosscheck.do — newfeat-style external-benchmark cross-check: wooldid vs jwdid
*   jwdid = Rios-Avila's implementation of Wooldridge's extended TWFE (Mundlak) estimator.
* Usage:  stata-mp -e do jwdid_crosscheck.do <path-to-wooldid.ado> <path-to-test/data>
*         (run_full.sh's newfeat component supplies both; see the locals below)
* Sentinels: RES|<leg>|<item>|<value>|<status>   status in PASS FAIL INFO
* Success line at the end: JWDID_XCHECK_ALL_PASS   (skip line: JWDID_XCHECK_SKIPPED)
* mreldif convention throughout: max |wooldid - jwdid| / (|jwdid| + 1).
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 20260824
set sortseed 20260824

* run_full.sh's newfeat component passes TWO positional arguments to every script here:
*   1 = <repo>/wooldidinstall/wooldid.ado      2 = <repo>/test/data
* The suite fixtures this script also needs live one level over, in <repo>/test/suite/data,
* which is why the sdir local is derived from ddir rather than hard-coded.
local wado `1'
local ddir `2'
if `"`wado'"' == "" local wado "../../../wooldidinstall/wooldid.ado"
if `"`ddir'"' == "" local ddir "../../../test/data"
local sdir `"`ddir'/../suite/data"'
qui do `"`wado'"'

* ---- preflight: jwdid is an optional external dependency (as ppmlhdfe is in reghdfe_oracle.do)
cap which jwdid
if _rc {
    dis "JWDID_XCHECK_SKIPPED (jwdid not installed)"
    exit
}
which jwdid
which jwdid_estat

global pf 0
global nres 0
**************************** helpers ****************************
* report one comparison
cap program drop res
program define res
  args leg item value tol vinfo
  * tol == . -> INFO only
  global nres = $nres + 1
  if "`tol'" == "." | "`tol'" == "" | "`vinfo'" == "info" {
    dis "RES|`leg'|`item'|" %12.4e `value' "|INFO"
  }
  else if `value' > `tol' {
    dis "RES|`leg'|`item'|" %12.4e `value' "|FAIL (tol `tol')"
    global pf = $pf + 1
  }
  else dis "RES|`leg'|`item'|" %12.4e `value' "|PASS (tol `tol')"
end

* after wooldid (with gteffects): frames WD (g t wb wse wgt), WDG (g wbg), WDC (t wbc);
* matrices WB WV (ES b/V); scalars watt wse_att; wN
cap program drop wd_stash
program define wd_stash
  tempname GT M
  matrix `GT' = e(gteffects)
  matrix `M' = e(wdidmainresults)
  scalar watt = `M'[1,1]
  scalar wse_att = `M'[1,2]
  scalar wN = e(N)
  cap matrix drop WB
  cap matrix drop WV
  cap matrix WB = e(b_es_att)
  cap matrix WV = e(V_es_att)
  cap frame drop __gt
  frame create __gt
  frame __gt {
    qui svmat double `GT', names(col)
    keep i onset caltime relyr estimate se nobs w_att
    * aggregate unit-level cells to (onset, caltime) with w_att weights (identity when i is the cohort)
    gen double num = w_att*estimate
    gen double senum = w_att*se
    collapse (sum) num senum wgt=w_att nobs (count) ncell=estimate (mean) relyr, by(onset caltime)
    gen double wb = num/wgt
    gen double wse = senum/wgt
    rename onset g
    rename caltime t
    keep g t relyr wb wse wgt nobs ncell
    sort g t
  }
  cap frame drop WD
  frame copy __gt WD
  * post-cell aggregates by cohort and by calendar time
  cap frame drop WDG
  frame copy WD WDG
  frame WDG {
    keep if relyr >= 0
    gen double num = wgt*wb
    collapse (sum) num wgt, by(g)
    gen double wbg = num/wgt
    keep g wbg
  }
  cap frame drop WDC
  frame copy WD WDC
  frame WDC {
    keep if relyr >= 0
    gen double num = wgt*wb
    collapse (sum) num wgt, by(t)
    gen double wbc = num/wgt
    keep t wbc
  }
end

* after a jwdid fit: frame JW (g t jb jse); matrices JB JV JKEY
cap program drop jw_cells
program define jw_cells
  tempname B V
  matrix `B' = e(b)
  matrix `V' = e(V)
  local names : colnames `B'
  local K = 0
  local keep
  local gl
  local tl
  foreach nm of local names {
    if regexm("`nm'", "^([0-9]+)[bn]*\.g#([0-9]+)[bn]*\.t#c\.__tr__$") {
      local K = `K' + 1
      local gl `gl' `=regexs(1)'
      local tl `tl' `=regexs(2)'
      local keep `keep' `nm'
    }
  }
  cap matrix drop JB JV JKEY
  matrix JB = J(1, `K', .)
  matrix JV = J(`K', `K', .)
  matrix JKEY = J(`K', 2, .)
  forval j = 1/`K' {
    local nj : word `j' of `keep'
    local cj = colnumb(`B', "`nj'")
    matrix JB[1, `j'] = `B'[1, `cj']
    matrix JKEY[`j', 1] = `: word `j' of `gl''
    matrix JKEY[`j', 2] = `: word `j' of `tl''
    forval k = 1/`K' {
      local nk : word `k' of `keep'
      local ck = colnumb(`B', "`nk'")
      matrix JV[`j', `k'] = `V'[`cj', `ck']
    }
  }
  cap frame drop JW
  frame create JW double(g t jb jse)
  forval j = 1/`K' {
    frame post JW (JKEY[`j',1]) (JKEY[`j',2]) (JB[1,`j']) (sqrt(JV[`j',`j']))
  }
  frame JW: sort g t
  dis "jw_cells: K=`K' cells"
end

* compare frames WD and JW on (g,t): cell b and se
cap program drop cmp_cells
program define cmp_cells
  args leg tolb tolse seinfo
  cap frame drop __c
  frame copy JW __c
  frame __c {
    qui frlink 1:1 g t, frame(WD)
    qui count if WD == .
    local nomatch = r(N)
    qui frget wb wse wgt ncell, from(WD)
    qui gen double db = abs(wb - jb)/(abs(jb) + 1)
    qui gen double dse = abs(wse - jse)/(abs(jse) + 1)
    qui sum db, meanonly
    local mb = r(max)
    qui sum dse, meanonly
    local ms = r(max)
    qui count
    local nj = r(N)
  }
  frame WD: qui count
  local nw = r(N)
  dis "cmp_cells `leg': jwdid cells=`nj' wooldid (aggregated) cells=`nw' unmatched=`nomatch'"
  if `nomatch' > 0 | `nj' != `nw' {
    dis "RES|`leg'|cell_count|" `nj' - `nw' "|FAIL (jwdid `nj' vs wooldid `nw', unmatched `nomatch')"
    global pf = $pf + 1
  }
  res `leg' cells_b `mb' `tolb'
  res `leg' cells_se `ms' `tolse' `seinfo'
end

* after estat event: align r(b)/r(V) to WB's colnames; rmin shifts index -> rel
cap program drop es_align
program define es_align
  args rmin
  tempname JB JV
  matrix `JB' = r(b)
  matrix `JV' = r(V)
  local jn : colnames `JB'
  local wl : colnames WB
  local K : word count `wl'
  local J : word count `jn'
  cap matrix drop OB OV
  matrix OB = J(1, `K', .)
  matrix OV = J(`K', `K', .)
  local cidx
  forval j = 1/`K' {
    local l : word `j' of `wl'
    if substr("`l'", 1, 3) == "pre" local r = -real(substr("`l'", 4, .))
    else local r = real(substr("`l'", 5, .))
    local cj = 0
    forval c = 1/`J' {
      local nm : word `c' of `jn'
      if regexm("`nm'", "^([0-9]+)[bn]*\.__event__$") {
        if real(regexs(1)) + `rmin' == `r' local cj = `c'
      }
    }
    if `cj' == 0 dis as error "es_align: no jwdid column for rel `r' (label `l')"
    local cidx `cidx' `cj'
  }
  forval j = 1/`K' {
    local cj : word `j' of `cidx'
    if `cj' == 0 continue
    matrix OB[1, `j'] = `JB'[1, `cj']
    forval k = 1/`K' {
      local ck : word `k' of `cidx'
      if `ck' == 0 continue
      matrix OV[`j', `k'] = `JV'[`cj', `ck']
    }
  }
end

* after estat group / calendar: compare r(b) to frame WDG (g wbg) or WDC (t wbc)
cap program drop cmp_over
program define cmp_over
  args leg kind tol
  tempname JB
  matrix `JB' = r(b)
  local jn : colnames `JB'
  local J : word count `jn'
  if "`kind'" == "group" {
    local fr WDG
    local key g
    local wv wbg
    local pat "^([0-9]+)[bn]*\.__group__$"
  }
  else {
    local fr WDC
    local key t
    local wv wbc
    local pat "^([0-9]+)[bn]*\.__calendar__$"
  }
  local mx = 0
  local nmatched = 0
  forval c = 1/`J' {
    local nm : word `c' of `jn'
    local lev = .
    if regexm("`nm'", "`pat'") local lev = real(regexs(1))
    else if "`nm'" == "r2vs1._at" & `J' == 1 {
      frame `fr': qui count
      if r(N) == 1 {
        frame `fr': local lev = `key'[1]
      }
    }
    if `lev' == . {
      dis as error "cmp_over: unparsed colname `nm'"
      continue
    }
    frame `fr' {
      qui count if `key' == `lev'
      if r(N) != 1 {
        dis as error "cmp_over: level `lev' not unique in `fr' (n=" r(N) ")"
        continue
      }
      qui sum `wv' if `key' == `lev', meanonly
      local wval = r(mean)
    }
    local d = abs(`wval' - `JB'[1, `c'])/(abs(`JB'[1, `c']) + 1)
    if `d' > `mx' local mx = `d'
    local nmatched = `nmatched' + 1
  }
  dis "cmp_over `leg' `kind': jwdid levels=`J' matched=`nmatched'"
  res `leg' `kind'_b `mx' `tol'
end

* full comparison after a jwdid fit: cells, simple, event, group, calendar
* args: leg mode(never|notyet) tolb tolv [seinfo] [predopt]
cap program drop jw_compare
program define jw_compare
  args leg mode tolb tolv seinfo predopt
  local jN = e(N)
  if reldif(wN, `jN') > 0 {
    dis "RES|`leg'|N|" wN - `jN' "|FAIL (wooldid N=" wN " jwdid N=`jN')"
    global pf = $pf + 1
  }
  else dis "RES|`leg'|N|" wN "|PASS (equal N)"
  jw_cells
  cmp_cells `leg' `tolb' `tolv' `seinfo'
  * rmin for never mode
  local rmin = 0
  if "`mode'" == "never" {
    tempvar ev
    qui gen double `ev' = t - g if g > 0 & e(sample)
    qui sum `ev', meanonly
    local rmin = r(min)
  }
  timer on 90
  cap noi estat simple, `predopt'
  timer off 90
  if _rc != 0 {
    dis "RES|`leg'|simple_rc|" _rc "|FAIL"
    global pf = $pf + 1
  }
  else {
    tempname sb sv
    matrix `sb' = r(b)
    matrix `sv' = r(V)
    res `leg' att_b `=reldif(watt, `sb'[1,1])' `tolb'
    res `leg' att_se `=reldif(wse_att, sqrt(`sv'[1,1]))' `tolv' `seinfo'
    dis "   att: wooldid " %21.0g watt " (se " %21.0g wse_att ")  jwdid " %21.0g `sb'[1,1] " (se " %21.0g sqrt(`sv'[1,1]) ")"
  }
  timer on 91
  cap noi estat event, `predopt'
  timer off 91
  if _rc != 0 {
    dis "RES|`leg'|event_rc|" _rc "|FAIL"
    global pf = $pf + 1
  }
  else {
    es_align `rmin'
    res `leg' es_b `=mreldif(WB, OB)' `tolb'
    res `leg' es_V `=mreldif(WV, OV)' `tolv' `seinfo'
    * per-horizon se ratio range (diagnostic)
    local K = colsof(WB)
    local rmn = 10
    local rmx = 0
    forval j = 1/`K' {
      local rr = sqrt(WV[`j',`j'])/sqrt(OV[`j',`j'])
      if `rr' < `rmn' local rmn = `rr'
      if `rr' > `rmx' local rmx = `rr'
    }
    dis "   es se ratio wooldid/jwdid range: " %8.5f `rmn' " .. " %8.5f `rmx'
  }
  timer on 92
  cap noi estat group, `predopt'
  timer off 92
  if _rc != 0 {
    dis "RES|`leg'|group_rc|" _rc "|FAIL"
    global pf = $pf + 1
  }
  else cmp_over `leg' group `tolb'
  timer on 93
  cap noi estat calendar, `predopt'
  timer off 93
  if _rc != 0 {
    dis "RES|`leg'|calendar_rc|" _rc "|FAIL"
    global pf = $pf + 1
  }
  else cmp_over `leg' calendar `tolb'
end


* probe: does jwdid's replay line `syntax [anything(everything)], [example *]` reject weights?
cap program drop _wtprobe
program define _wtprobe
  syntax [anything(everything)], [example *]
  dis "wtprobe anything=[`anything'] options=[`options']"
end

* EMULATION of jwdid's linear path with weights (jwdid v1.41 rejects weights at its
* replay-parse line). Requires: an UNWEIGHTED jwdid fit just ran on the current data
* (leaves __tr__ in the data; JKEY from jw_cells gives the exact cell list).
* Runs  reghdfe y <cells> [pw=w], abs(<ivar> t) cluster(i) keepsingletons  (jwdid.ado lines 202-203)
* then re-implements estat simple/event as explicit weight vectors over the cells
* (estat = margins subpop(__etr__==1) at(__tr__=(0 1)) contrast, i.e. the estimation-
* weight-weighted mean of the cell coefficients over treated obs).
cap program drop jw_emul_compare
program define jw_emul_compare
  args leg ivar tolb tolv seinfo
  jw_cells
  local K = rowsof(JKEY)
  local cells
  forval j = 1/`K' {
    local cells `cells' c.__tr__#i`=JKEY[`j',1]'.g#i`=JKEY[`j',2]'.t
  }
  cap noi reghdfe y `cells' [pw=w], abs(`ivar' t) cluster(i) keepsingletons
  local rrc = _rc
  if `rrc' != 0 {
    dis "RES|`leg'|emul_rc|`rrc'|FAIL (emulated weighted reghdfe failed)"
    global pf = $pf + 1
    exit
  }
  local jN = e(N)
  if reldif(wN, `jN') > 0 {
    dis "RES|`leg'|N|" wN - `jN' "|FAIL (wooldid N=" wN " emul N=`jN')"
    global pf = $pf + 1
  }
  else dis "RES|`leg'|N|" wN "|PASS (equal N)"
  jw_cells
  cmp_cells `leg' `tolb' `tolv' `seinfo'
  * cell weight sums
  local K = rowsof(JKEY)
  tempname WS L LES
  matrix `WS' = J(1, `K', 0)
  forval j = 1/`K' {
    qui sum w if g == JKEY[`j',1] & t == JKEY[`j',2] & e(sample)
    matrix `WS'[1, `j'] = r(sum)
  }
  * att: post cells
  matrix `L' = J(1, `K', 0)
  local tot = 0
  forval j = 1/`K' {
    if JKEY[`j',2] >= JKEY[`j',1] {
      matrix `L'[1, `j'] = `WS'[1, `j']
      local tot = `tot' + `WS'[1, `j']
    }
  }
  matrix `L' = `L'/`tot'
  tempname ab av
  matrix `ab' = `L'*JB'
  matrix `av' = `L'*JV*`L''
  res `leg' att_b `=reldif(watt, `ab'[1,1])' `tolb'
  res `leg' att_se `=reldif(wse_att, sqrt(`av'[1,1]))' `tolv' `seinfo'
  dis "   att: wooldid " %21.0g watt " (se " %21.0g wse_att ")  emul " %21.0g `ab'[1,1] " (se " %21.0g sqrt(`av'[1,1]) ")"
  * event-study horizons aligned to WB labels
  local wl : colnames WB
  local M : word count `wl'
  matrix `LES' = J(`M', `K', 0)
  forval m = 1/`M' {
    local l : word `m' of `wl'
    if substr("`l'", 1, 3) == "pre" local r = -real(substr("`l'", 4, .))
    else local r = real(substr("`l'", 5, .))
    local tot = 0
    forval j = 1/`K' {
      if JKEY[`j',2] - JKEY[`j',1] == `r' {
        matrix `LES'[`m', `j'] = `WS'[1, `j']
        local tot = `tot' + `WS'[1, `j']
      }
    }
    if `tot' > 0 {
      forval j = 1/`K' {
        matrix `LES'[`m', `j'] = `LES'[`m', `j']/`tot'
      }
    }
    else dis as error "emul: no cells for rel `r'"
  }
  cap matrix drop OB OV
  matrix OB = (`LES'*JB')'
  matrix OV = `LES'*JV*`LES''
  res `leg' es_b `=mreldif(WB, OB)' `tolb'
  res `leg' es_V `=mreldif(WV, OV)' `tolv' `seinfo'
end


**************************** data: synthetic common-onset panel (40 x 12) ****************************
set obs 40
gen double i = _n
gen double a = rnormal()
gen double ttre = cond(i <= 20, 7, .)
expand 12
bysort i: gen double t = _n
gen double dt = 0.2*t + 0.5*cos(t)
gen double rel = t - ttre
gen double tau = cond(ttre < . & t >= ttre, 0.5 + 0.06*(t - ttre), 0)
gen double y = a + dt + tau + rnormal()*0.4
gen double icoh = cond(ttre < ., 1, i)
gen double g   = cond(ttre < ., ttre, 0)
gen double g2  = cond(ttre < ., 1, 2)
gen double w   = 0.5 + runiform()*2
tempfile sbase
save `sbase'
* unbalanced variant: drop ~15% of obs, never the rel=-1 period (t=6)
gen double u = runiform()
drop if u < 0.15 & t != 6
drop u
tempfile ubase
save `ubase'
timer clear

* ---- J1: syn, icoh, fixed base (-1) 4/4  <->  jwdid never ivar(icoh)
dis "=== J1 ==="
use `sbase', clear
timer on 1
cap noi wooldid y icoh t ttre, att cluster(i) esfixedbaseperiod esprelength(4) espostlength(4) gteffects
timer off 1
dis "J1_W_RC=" _rc
wd_stash
use `sbase', clear
timer on 2
cap noi jwdid y, ivar(icoh) tvar(t) gvar(g) never cluster(i)
timer off 2
local jrc = _rc
dis "J1_J_RC=`jrc'"
if `jrc' == 0 jw_compare J1 never 1e-8 1e-6
else {
  dis "RES|J1|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
  global pf = $pf + 1
}

* ---- J2: syn, icoh, POOLED reference (wooldid default)  <->  jwdid DEFAULT (no never)
dis "=== J2 ==="
use `sbase', clear
timer on 3
cap noi wooldid y icoh t ttre, att cluster(i) espostlength(4) gteffects
timer off 3
dis "J2_W_RC=" _rc
wd_stash
use `sbase', clear
timer on 4
cap noi jwdid y, ivar(icoh) tvar(t) gvar(g) cluster(i)
timer off 4
local jrc = _rc
dis "J2_J_RC=`jrc'"
if `jrc' == 0 jw_compare J2 notyet 1e-8 1e-6
else {
  dis "RES|J2|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
  global pf = $pf + 1
}

* ---- J8b: UNBALANCED syn, g2 cohort-FE fixed base  <->  jwdid never group NOCORR
dis "=== J8b ==="
use `ubase', clear
timer on 5
cap noi wooldid y g2 t ttre, att cluster(i) esfixedbaseperiod esprelength(4) espostlength(4) gteffects
timer off 5
dis "J8b_W_RC=" _rc
wd_stash
use `ubase', clear
timer on 6
cap noi jwdid y, ivar(i) tvar(t) gvar(g) never group nocorr cluster(i)
timer off 6
local jrc = _rc
dis "J8b_J_RC=`jrc'"
if `jrc' == 0 jw_compare J8b never 1e-8 1e-6
else {
  dis "RES|J8b|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
  global pf = $pf + 1
}

**************************** data: dta_binary (staggered 2005 / 2008 / never) ****************************
use `"`ddir'/dta_binary.dta"', clear
gen double icoh = cond(ttre < ., ttre, i)
gen double g    = cond(ttre < ., ttre, 0)
gen double g2   = cond(ttre < ., ttre, 1)
tempfile dbase
save `dbase'

* ---- J3: dta_binary, icoh, fixed base 7/7  <->  jwdid never ivar(icoh)
dis "=== J3 ==="
use `dbase', clear
timer on 7
cap noi wooldid y icoh t ttre, att cluster(i) esfixedbaseperiod esprelength(7) espostlength(7) gteffects
timer off 7
dis "J3_W_RC=" _rc
wd_stash
use `dbase', clear
timer on 8
cap noi jwdid y, ivar(icoh) tvar(t) gvar(g) never cluster(i)
timer off 8
local jrc = _rc
dis "J3_J_RC=`jrc'"
if `jrc' == 0 jw_compare J3 never 1e-8 1e-6
else {
  dis "RES|J3|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
  global pf = $pf + 1
}

* ---- J4: dta_binary, icoh, pooled reference  <->  jwdid default
dis "=== J4 ==="
use `dbase', clear
timer on 9
cap noi wooldid y icoh t ttre, att cluster(i) espostlength(7) gteffects
timer off 9
dis "J4_W_RC=" _rc
wd_stash
use `dbase', clear
timer on 10
cap noi jwdid y, ivar(icoh) tvar(t) gvar(g) cluster(i)
timer off 10
local jrc = _rc
dis "J4_J_RC=`jrc'"
if `jrc' == 0 jw_compare J4 notyet 1e-8 1e-6
else {
  dis "RES|J4|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
  global pf = $pf + 1
}

* ---- J5: dta_binary, wooldid's LITERAL default i = unit id, fixed base  <->  jwdid never ivar(i)
*          b judged at 1e-8 after w_att aggregation; V is INFO only (documented structural difference)
dis "=== J5 ==="
use `dbase', clear
timer on 11
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esprelength(7) espostlength(7) gteffects
timer off 11
dis "J5_W_RC=" _rc
wd_stash
use `dbase', clear
timer on 12
cap noi jwdid y, ivar(i) tvar(t) gvar(g) never cluster(i)
timer off 12
local jrc = _rc
dis "J5_J_RC=`jrc'"
if `jrc' == 0 jw_compare J5 never 1e-8 1e-6 info
else {
  dis "RES|J5|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
  global pf = $pf + 1
}

* ---- J6: dta_binary, icoh fixed base [aw=w]  <->  EMULATION of jwdid's linear path with [pw=w]
*          (jwdid v1.41 rejects weights at its replay-parse line: r(101))
dis "=== J6 ==="
use `dbase', clear
cap noi jwdid y [pw=w], ivar(icoh) tvar(t) gvar(g) never cluster(i)
dis "J6_JWDID_PW_RC=" _rc
use `dbase', clear
timer on 13
cap noi wooldid y icoh t ttre [aw=w], att cluster(i) esfixedbaseperiod esprelength(7) espostlength(7) gteffects
timer off 13
dis "J6_W_RC=" _rc
wd_stash
use `dbase', clear
timer on 14
cap noi jwdid y, ivar(icoh) tvar(t) gvar(g) never cluster(i)
timer off 14
local jrc = _rc
dis "J6_J_RC=`jrc'"
if `jrc' == 0 jw_emul_compare J6 icoh 1e-8 1e-6
else {
  dis "RES|J6|jwdid_rc|`jrc'|FAIL"
  global pf = $pf + 1
}

* ---- J7: dta_binary poisson, g2 cohort-FE fixed base  <->  jwdid never method(poisson), xb scale
dis "=== J7 ==="
use `dbase', clear
timer on 15
cap noi wooldid ycount g2 t ttre, att cluster(i) poisson esfixedbaseperiod esprelength(7) espostlength(7) gteffects
timer off 15
dis "J7_W_RC=" _rc
wd_stash
use `dbase', clear
timer on 16
cap noi jwdid ycount, ivar(i) tvar(t) gvar(g) never cluster(i) method(poisson)
timer off 16
local jrc = _rc
dis "J7_J_RC=`jrc'"
if `jrc' == 0 jw_compare J7 never 1e-6 1e-6 "" predict(xb)
else {
  dis "RES|J7|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
  global pf = $pf + 1
}

* ---- J8: dta_binary, g2 cohort-FE fixed base  <->  jwdid never group ivar(i)
dis "=== J8 ==="
use `dbase', clear
timer on 17
cap noi wooldid y g2 t ttre, att cluster(i) esfixedbaseperiod esprelength(7) espostlength(7) gteffects
timer off 17
dis "J8_W_RC=" _rc
wd_stash
use `dbase', clear
timer on 18
cap noi jwdid y, ivar(i) tvar(t) gvar(g) never group cluster(i)
timer off 18
local jrc = _rc
dis "J8_J_RC=`jrc'"
if `jrc' == 0 jw_compare J8 never 1e-8 1e-6
else {
  dis "RES|J8|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
  global pf = $pf + 1
}

**************************** data: dta_s_admin (realistic administrative panel) ****************************
* ~2,000 FIPS-style counties x 2005-2019, onsets 2009 / 2012 / 2015 + never, unbalanced
* (late entry / early exit), ~2% missing y.  Renamed to the harness's y i t g w names.
capture confirm file `"`sdir'/dta_s_admin.dta"'
if _rc {
    dis "RES|J9|skipped|0|INFO (dta_s_admin.dta absent)"
}
else {
    use `"`sdir'/dta_s_admin.dta"', clear
    rename fips i
    rename year t
    gen double g    = cond(ttre < ., ttre, 0)
    * fips codes overlap the onset years (1001..48999), so shift the control ids
    gen double icoh = cond(ttre < ., ttre, i + 100000)
    gen double w    = pop
    tempfile abase
    save `abase'

    * ---- J9: admin, icoh fixed base 4/4  <->  jwdid never ivar(icoh)
    dis "=== J9 ==="
    use `abase', clear
    timer on 19
    cap noi wooldid y icoh t ttre, att cluster(i) esfixedbaseperiod esprelength(4) espostlength(4) gteffects
    timer off 19
    dis "J9_W_RC=" _rc
    wd_stash
    use `abase', clear
    timer on 20
    cap noi jwdid y, ivar(icoh) tvar(t) gvar(g) never cluster(i)
    timer off 20
    local jrc = _rc
    dis "J9_J_RC=`jrc'"
    if `jrc' == 0 jw_compare J9 never 1e-8 1e-6
    else {
      dis "RES|J9|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
      global pf = $pf + 1
    }

    * ---- J10: admin, icoh pooled reference  <->  jwdid default
    dis "=== J10 ==="
    use `abase', clear
    timer on 21
    cap noi wooldid y icoh t ttre, att cluster(i) espostlength(4) gteffects
    timer off 21
    dis "J10_W_RC=" _rc
    wd_stash
    use `abase', clear
    timer on 22
    cap noi jwdid y, ivar(icoh) tvar(t) gvar(g) cluster(i)
    timer off 22
    local jrc = _rc
    dis "J10_J_RC=`jrc'"
    if `jrc' == 0 jw_compare J10 notyet 1e-8 1e-6
    else {
      dis "RES|J10|jwdid_rc|`jrc'|FAIL (jwdid did not run)"
      global pf = $pf + 1
    }
}

timer list
dis "JWDID_XCHECK_NRES=$nres  NFAIL=$pf"
if $pf == 0 dis "JWDID_XCHECK_ALL_PASS"
else dis "JWDID_XCHECK_FAILED"
