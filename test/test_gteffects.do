* Targeted suite for F3 gteffects: the cohort-by-relative-period building-block table.
* Run from test/:  "<stata>" -e do test_gteffects.do   -> expect GTEFFECTS_ALL_PASS
* Core invariant: the weighted mean of the estimate column under each estimand's weight
* column reproduces that aggregate estimate (the cells ARE the building blocks).
clear all
set more off
set maxvar 20000
set graphics off
set seed 12345
set sortseed 12345
local DATA "data"
do "../wooldidinstall/wooldid.ado"

local fails 0

* ---- G1: basic binary run - shape, names, aggregation identity, sane SEs ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) gteffects
if _rc != 0 {
  dis as error "G1 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G
  matrix `G' = e(gteffects)
  if colsof(`G') != 15 {
    dis as error "G1 FAIL: colsof=" colsof(`G')
    local ++fails
  }
  local cn : colnames `G'
  if "`cn'" != "i onset caltime relyr estimate se p lb_95ci ub_95ci inwindow nobs w_att w_att_it w_att_i w_att_itime" {
    dis as error "G1 FAIL: colnames=[`cn']"
    local ++fails
  }
  local num 0
  local den 0
  local badse 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',4] >= 0 {
      local num = `num' + `G'[`r',12]*`G'[`r',5]
      local den = `den' + `G'[`r',12]
    }
    if `G'[`r',6] != . & `G'[`r',6] <= 0 local ++badse
  }
  if reldif(`num'/`den', _b[main_att]) > 1e-10 {
    dis as error "G1 FAIL: att aggregation identity reldif=" reldif(`num'/`den', _b[main_att])
    local ++fails
  }
  if `badse' != 0 {
    dis as error "G1 FAIL: `badse' nonpositive SEs"
    local ++fails
  }
}

* ---- G2: aw weights - att vs att_it differ, both identities hold, horizon cross-check ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre [aw=w], att att_it cluster(i) esprelength(2) espostlength(2) gteffects
if _rc != 0 {
  dis as error "G2 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G BE
  matrix `G' = e(gteffects)
  matrix `BE' = e(b_es_att)
  local num 0
  local den 0
  local numit 0
  local denit 0
  local h1 0
  local hw 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',4] >= 0 {
      local num = `num' + `G'[`r',12]*`G'[`r',5]
      local den = `den' + `G'[`r',12]
      local numit = `numit' + `G'[`r',13]*`G'[`r',5]
      local denit = `denit' + `G'[`r',13]
    }
    if `G'[`r',4] == 1 & `G'[`r',10] == 1 {
      local h1 = `h1' + `G'[`r',12]*`G'[`r',5]
      local hw = `hw' + `G'[`r',12]
    }
  }
  if reldif(`num'/`den', _b[main_att]) > 1e-10 | reldif(`numit'/`denit', _b[main_att_it]) > 1e-10 {
    dis as error "G2 FAIL: identities att=" reldif(`num'/`den', _b[main_att]) " att_it=" reldif(`numit'/`denit', _b[main_att_it])
    local ++fails
  }
  if reldif(_b[main_att], _b[main_att_it]) <= 1e-12 {
    dis as error "G2 FAIL: att and att_it coincide under aw (test has no power)"
    local ++fails
  }
  local c1 = colnumb(`BE', "post1")
  if reldif(`h1'/`hw', `BE'[1,`c1']) > 1e-10 {
    dis as error "G2 FAIL: post1 horizon cross-check reldif=" reldif(`h1'/`hw', `BE'[1,`c1'])
    local ++fails
  }
}

* ---- G3: esfixedbaseperiod - relyr mapping, caltime identity, window flag ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esfixedbaseperiod esrelativeto(-2) esprelength(3) espostlength(2) gteffects
if _rc != 0 {
  dis as error "G3 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G BE
  matrix `G' = e(gteffects)
  matrix `BE' = e(b_es_att)
  local h0 0
  local hw0 0
  local calbad 0
  local inw 0
  local outw 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',4] == 0 & `G'[`r',10] == 1 {
      local h0 = `h0' + `G'[`r',12]*`G'[`r',5]
      local hw0 = `hw0' + `G'[`r',12]
    }
    if `G'[`r',3] != `G'[`r',2] + `G'[`r',4] local ++calbad
    if `G'[`r',10] == 1 local ++inw
    else local ++outw
  }
  local c0 = colnumb(`BE', "post0")
  if reldif(`h0'/`hw0', `BE'[1,`c0']) > 1e-10 {
    dis as error "G3 FAIL: post0 horizon cross-check reldif=" reldif(`h0'/`hw0', `BE'[1,`c0'])
    local ++fails
  }
  if `calbad' != 0 {
    dis as error "G3 FAIL: caltime != onset + relyr on `calbad' rows"
    local ++fails
  }
  if `inw' == 0 | `outw' == 0 {
    dis as error "G3 FAIL: inwindow split degenerate (inw=`inw' outw=`outw')"
    local ++fails
  }
}

* ---- G4: contreat() is refused at parse ----
use "`DATA'/dta_cont", clear
cap wooldid y i t ttre, contreat(intensity) contreatcontrols(basic) cluster(i) gteffects
if _rc != 198 {
  dis as error "G4 FAIL: rc=" _rc " (expected 198)"
  local ++fails
}

* ---- G5: fastnose synergy - full SEs in the table while e(V) is absent ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) fastnose gteffects
if _rc != 0 {
  dis as error "G5 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G
  matrix `G' = e(gteffects)
  local misse 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',6] == . local ++misse
  }
  cap confirm matrix e(V)
  if `misse' != 0 | _rc == 0 {
    dis as error "G5 FAIL: missing SEs=`misse', e(V) present rc=" _rc
    local ++fails
  }
}

* ---- G6: timetrends - degenerate variances report as missing, never zero ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) timetrends gteffects
if _rc != 0 {
  dis as error "G6 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G
  matrix `G' = e(gteffects)
  local zse 0
  forval r = 1/`=rowsof(`G')' {
    if `G'[`r',6] == 0 local ++zse
  }
  if `zse' != 0 {
    dis as error "G6 FAIL: `zse' exact-zero SEs shipped"
    local ++fails
  }
}

* ---- G7: poisson - runs; the xb-scale note under the table was removed by the author's 2026-09-21 message review ----
use "`DATA'/dta_binary", clear
qui log using gt_t7.txt, text replace name(gt)
cap noi wooldid ycount i t ttre, att cluster(i) poisson gteffects
local rc7 = _rc
qui log close gt
if `rc7' != 0 | strpos(fileread("gt_t7.txt"), "linear-index (xb) scale") != 0 {
  dis as error "G7 FAIL: rc=" `rc7' " xb note printed but was removed"
  local ++fails
}

* ---- G8: subgroups - one full-sample table ----
use "`DATA'/dta_binary", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) subgroups(sg) gteffects
cap confirm matrix e(gteffects)
if _rc != 0 {
  dis as error "G8 FAIL: no table under subgroups"
  local ++fails
}

* ---- G9: 9-digit administrative ids round-trip exactly ----
use "`DATA'/dta_binary", clear
qui gen double bigi = 100000000 + i*7919
cap noi wooldid y bigi t ttre, att cluster(bigi) esprelength(2) espostlength(2) gteffects
if _rc != 0 {
  dis as error "G9 FAIL: rc=" _rc
  local ++fails
}
else {
  tempname G
  matrix `G' = e(gteffects)
  qui glevelsof bigi if ttre != ., clean local(bl)
  local bad 0
  forval r = 1/`=rowsof(`G')' {
    local v = `G'[`r',1]
    if strpos(" `bl' ", " `=string(`v', "%12.0f")' ") == 0 local ++bad
  }
  if `bad' != 0 {
    dis as error "G9 FAIL: `bad' id values do not round-trip"
    local ++fails
  }
}

* ---- G10: replay reprints the table; G11: off-state leaves e() clean ----
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) gteffects
qui log using gt_t10.txt, text replace name(gt)
cap noi wooldid
qui log close gt
if _rc != 0 | strpos(fileread("gt_t10.txt"), "Cohort-by-Relative-Period Effects") == 0 {
  dis as error "G10 FAIL: replay rc=" _rc " or table not reprinted"
  local ++fails
}
use "`DATA'/dta_binary", clear
qui wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2)
cap confirm matrix e(gteffects)
if _rc == 0 {
  dis as error "G11 FAIL: e(gteffects) posted without the option"
  local ++fails
}

* ---- G12 (change audit 2026-09-22, gt-containment-customvce-arm): the containment count and
* its caution under a heteroskedasticity-robust customvce(). every treated cell of
* dta_path_onepercell is one observation, so robust counts all of them; customvce(robust) is the
* same meat and must count the same cells (it used to post no e(gt_ncontained) and print nothing) ----
set linesize 255
use "`DATA'/dta_path_onepercell", clear
cap noi wooldid y i t ttre, att robust esprelength(1) espostlength(1) gteffects
local rc12a = _rc
local nc12a = e(gt_ncontained)
local nr12a = cond(`rc12a' == 0, rowsof(e(gteffects)), .)
use "`DATA'/dta_path_onepercell", clear
qui log using gt_g12.txt, text replace name(gt12)
cap noi wooldid y i t ttre, att customvce(robust) esprelength(1) espostlength(1) gteffects
local rc12b = _rc
qui log close gt12
local nc12b = e(gt_ncontained)
local cau12b = strpos(fileread("gt_g12.txt"), "e(gteffects) cell(s) lie wholly inside a single observation") > 0
cap erase gt_g12.txt
dis "G12: robust rc=`rc12a' ncontained=`nc12a' rows=`nr12a' | customvce(robust) rc=`rc12b' ncontained=`nc12b' caution=`cau12b'"
if `rc12a' != 0 | `rc12b' != 0 | `nc12a' >= . | `nc12a' < 1 | `nc12b' != `nc12a' | `cau12b' != 1 {
  dis as error "G12 FAIL: customvce(robust) must count the one-observation cells as robust does"
  local ++fails
}

* ---- G12c (author ruling 2026-09-25, POIS-2): poisson on the ppmlhdfe path estimates a robust
* variance unasked, so its one-row cells have robust's blind spot and must be counted as robust
* counts them. a strictly positive count outcome keeps every one-row cell out of separation.
* G12d: with cluster() the cluster branch still counts (the arm sits below it). G12e: on the
* customvce(oim) path the variance is model-based and nothing is counted ----
use "`DATA'/dta_path_onepercell", clear
qui su y
qui gen long gtyc = round(y - r(min)) + 1
tempfile gt12pois
qui save `gt12pois'
qui log using gt_g12c.txt, text replace name(gt12c)
cap noi wooldid gtyc i t ttre, att poisson esprelength(1) espostlength(1) gteffects
local rc12c = _rc
qui log close gt12c
local nc12c = e(gt_ncontained)
local vce12c "`e(vce)'"
local cau12c = strpos(fileread("gt_g12c.txt"), "e(gteffects) cell(s) lie wholly inside a single observation") > 0
cap erase gt_g12c.txt
use `gt12pois', clear
qui log using gt_g12d.txt, text replace name(gt12d)
cap noi wooldid gtyc i t ttre, att poisson cluster(i) esprelength(1) espostlength(1) gteffects
local rc12d = _rc
qui log close gt12d
local nc12d = e(gt_ncontained)
local cau12d = strpos(fileread("gt_g12d.txt"), "e(gteffects) cell(s) lie wholly inside a single cluster") > 0
cap erase gt_g12d.txt
use `gt12pois', clear
cap noi wooldid gtyc i t ttre, att poisson customvce(oim) esprelength(1) espostlength(1) gteffects
local rc12e = _rc
local nc12e = e(gt_ncontained)
dis "G12c: poisson rc=`rc12c' ncontained=`nc12c' (robust count `nc12a') caution=`cau12c' vce=`vce12c' | G12d cluster(i) rc=`rc12d' ncontained=`nc12d' caution=`cau12d' | G12e customvce(oim) rc=`rc12e' ncontained=`nc12e'"
if `rc12c' != 0 | `nc12c' != `nc12a' | `cau12c' != 1 | "`vce12c'" != "robust" {
  dis as error "G12c FAIL: poisson on ppmlhdfe must count the one-observation cells as robust does"
  local ++fails
}
if `rc12d' != 0 | `nc12d' >= . | `nc12d' < 1 | `cau12d' != 1 {
  dis as error "G12d FAIL: poisson with cluster(i) must keep the cluster count and its caution"
  local ++fails
}
if `rc12e' != 0 | `nc12e' < . {
  dis as error "G12e FAIL: poisson with customvce(oim) must count nothing"
  local ++fails
}

* ---- G13 (round 2, gteffects-table-format-and-floor; author ruling 2026-09-23): the printed
* gteffects table takes the PS-3 per-column formats of the main, subgroup and pairwise tables
* (wooldid_mlfmt), so p prints to three decimals and relyr as an integer. It used
* format(%9.0g) on every column, so a p printed as .8122894 or 5.64e-34 (and a bootstrap 0 as
* 0). No p-floor line goes under this table (ruled): its cells rest on per-cell draw counts, so
* the "p = 0 above means" line stays under the tables whose floor it states. gt13scan reads
* every printed row of the table back from the log. linesize 255 keeps each row on one line ----
cap program drop gt13scan
program define gt13scan, rclass
  args fn
  local seen 0
  local rows 0
  local badp 0
  local badr 0
  local firstbad
  local flbefore 0
  local flafter 0
  tempname fh
  file open `fh' using `fn', read text
  file read `fh' ln
  while r(eof) == 0 {
    if strpos(`"`macval(ln)'"', "Cohort-by-Relative-Period Effects") > 0 local seen 1
    if strpos(`"`macval(ln)'"', "p = 0 above means") > 0 {
      if `seen' local ++flafter
      else local ++flbefore
    }
    if `seen' & regexm(`"`macval(ln)'"', "^ +r[0-9]+ [|]") {
      local ++rows
      local rest = substr(`"`macval(ln)'"', strpos(`"`macval(ln)'"', "|") + 1, .)
      local pt : word 7 of `rest'
      local rt : word 4 of `rest'
      if !regexm("`pt'", "^[0-9][.][0-9][0-9][0-9]$") & "`pt'" != "." {
        local ++badp
        if "`firstbad'" == "" local firstbad `pt'
      }
      if !regexm("`rt'", "^-?[0-9]+$") local ++badr
    }
    file read `fh' ln
  }
  file close `fh'
  return scalar rows = `rows'
  return scalar badp = `badp'
  return scalar badr = `badr'
  return scalar flbefore = `flbefore'
  return scalar flafter = `flafter'
  return local firstbad `firstbad'
end
set linesize 255
use "`DATA'/dta_simple1", clear
qui log using gt_g13.txt, text replace name(gt13)
cap noi wooldid y i t ttre, att att_it robust esprelength(2) espostlength(2) gteffects
local rc13a = _rc
qui log close gt13
local nr13a = cond(`rc13a' == 0, rowsof(e(gteffects)), .)
gt13scan gt_g13.txt
local rows13a = r(rows)
local badp13a = r(badp)
local badr13a = r(badr)
local fb13a `r(firstbad)'
cap erase gt_g13.txt
dis "G13a: rc=`rc13a' e() rows=`nr13a' printed rows=`rows13a' bad p=`badp13a' (first `fb13a') bad relyr=`badr13a'"
if `rc13a' != 0 | `nr13a' >= . | `rows13a' != `nr13a' | `badp13a' != 0 | `badr13a' != 0 {
  dis as error "G13a FAIL: the gteffects table must print p to three decimals and relyr as an integer"
  local ++fails
}
use "`DATA'/dta_simple1", clear
qui log using gt_g13.txt, text replace name(gt13)
cap noi wooldid y i t ttre, att att_it cluster(i) esprelength(2) espostlength(2) gteffects bootstrap(standard) reps(12) rseed(5) dots(0)
local rc13b = _rc
qui log close gt13
local nr13b = cond(`rc13b' == 0, rowsof(e(gteffects)), .)
gt13scan gt_g13.txt
local rows13b = r(rows)
local badp13b = r(badp)
local fb13b `r(firstbad)'
local flb13b = r(flbefore)
local fla13b = r(flafter)
cap erase gt_g13.txt
dis "G13b: rc=`rc13b' e() rows=`nr13b' printed rows=`rows13b' bad p=`badp13b' (first `fb13b') floor lines before/after the table=`flb13b'/`fla13b'"
if `rc13b' != 0 | `nr13b' >= . | `rows13b' != `nr13b' | `badp13b' != 0 | `flb13b' < 1 | `fla13b' != 0 {
  dis as error "G13b FAIL: under bootstrap() the gteffects p prints to three decimals, with no floor line of its own"
  local ++fails
}

* g5find <file> <needle>: r(found) is 1 when the needle occurs in the file once line breaks,
* spaces and the "> " continuation marks are removed from both, so a message the log wrapped
* still matches. Read only from logs closed before the call
cap program drop g5find
program define g5find, rclass
  args lf needle
  local h = fileread("`lf'")
  local h = subinstr(`"`macval(h)'"', char(13), "", .)
  local h = subinstr(`"`macval(h)'"', char(10), "", .)
  local h = subinstr(`"`macval(h)'"', ">", "", .)
  local h = subinstr(`"`macval(h)'"', " ", "", .)
  local n = subinstr(`"`needle'"', " ", "", .)
  return scalar found = strpos(`"`macval(h)'"', `"`n'"') > 0
end

* ---- G14 (author ruling 2026-09-25, card m055): the containment caution says bootstrap()
* does not repair a contained cell, since such a cell is resampled whole, and it keeps
* printing inside a bootstrap() run, whose observed run is not fastnose. cluster(i) holds
* every cohort-by-period cell of dta_binary inside one cluster, so every cell counts. The
* old last sentence, which named bootstrap() and permute() as remedies, must be gone ----
set linesize 255
local g14new "bootstrap() does not repair this, since such a cell is resampled whole. Only a cluster() strictly finer than the cohort-by-period cells gives valid cell inference. The aggregate rows are a separate question"
local g14old "bootstrap() or permute() gives valid cell inference"
foreach g14m in an bs {
  local g14b
  if "`g14m'" == "bs" local g14b bootstrap(standard) reps(12) rseed(5) dots(0)
  use "`DATA'/dta_binary", clear
  qui log using gt_g14.txt, text replace name(gt14)
  cap noi wooldid y i t ttre, att cluster(i) esprelength(2) espostlength(2) gteffects `g14b'
  local rc14`g14m' = _rc
  qui log close gt14
  local nc14`g14m' = e(gt_ncontained)
  g5find gt_g14.txt "e(gteffects) cell(s) lie wholly inside a single cluster"
  local op14`g14m' = r(found)
  g5find gt_g14.txt "`g14new'"
  local new14`g14m' = r(found)
  g5find gt_g14.txt "`g14old'"
  local old14`g14m' = r(found)
  cap erase gt_g14.txt
}
dis "G14: analytic rc=`rc14an' ncontained=`nc14an' opener=`op14an' new=`new14an' old=`old14an' | bootstrap rc=`rc14bs' ncontained=`nc14bs' opener=`op14bs' new=`new14bs' old=`old14bs'"
if `rc14an' != 0 | `rc14bs' != 0 | `nc14an' >= . | `nc14an' < 1 | `nc14bs' != `nc14an' | `op14an' != 1 | `op14bs' != 1 | `new14an' != 1 | `new14bs' != 1 | `old14an' != 0 | `old14bs' != 0 {
  dis as error "G14 FAIL: the containment caution must carry the 2026-09-25 text, analytic and under bootstrap()"
  local ++fails
}

if `fails' == 0 dis "GTEFFECTS_ALL_PASS"
else dis as error "GTEFFECTS_FAILS=`fails'"
exit, clear
