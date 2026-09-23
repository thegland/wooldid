* test_sgdegen.do - the realized-degeneracy screen on the subgroup pairwise covariance
* (bug sweep BS-15, author ruling 2026-09-07 "blank with a caution").
* Run from test/: stata-mp -e do test_sgdegen.do -> expect SGDEGEN_ALL_PASS
*
* When the clustering is at the level of i or above, the two subgroups' level standard
* errors agree to several digits and the contrast standard error is the cancellation
* residue between them, not an estimated variance. Measured on dta_binary.dta: the
* implied correlation is .99940 on the main row and .99715 on the pre-treatment row under
* cluster(i), against .51141 and .49779 under robust, and the clustered contrast standard
* errors are 1/44 and 1/20 of their robust counterparts.
*
* G1 the screen fires, blanks both pairwise rows and both equality tests, keeps the points
* G2 a healthy robust design is untouched, bit-for-bit
* G3 safety(off) reports instead of withholding
* G4 the BS-15 and jointtests cautions stand side by side, and jointtests still computes
* G5 no false flag on a healthy continuous-treatment subgroup run
* G6 replay reproduces the caution from e()
* G7 display-eqtest-bare-dot: a withheld equality test prints nothing, a reported one prints
*    its value (the string gate printed a bare dot for a scalar posted as missing)
* G8 under suppresssgprimaryeffects the screen fires exactly as in G1: it reads the subgroup
*    level rows of the pairwise passes, posted hidden as e(wdid_sglev)
* G9 the continuous-treatment twin under suppresssgprimaryeffects: the att rows are flagged as
*    on the unsuppressed run, and the ame rows, whose standard errors the cell-containment
*    scrub blanks, are skipped rather than counted
* G10 on a healthy continuous-treatment run under suppresssgprimaryeffects nothing fires, and
*    e(wdid_sglev) carries the level standard errors the unsuppressed run posts in
*    e(wdidsgresults_*), for the att and ame rows alike
clear all
set more off
do "../wooldidinstall/wooldid.ado"
local fails 0

cap program drop scanSG
program define scanSG, rclass
  args lf needle
  * the needles below are short enough to sit inside Stata's 80-column wrap, so a plain
  * line scan is enough and no continuation rejoining is needed - which matters here,
  * because the caution text carries apostrophes that break the rejoining idiom
  tempname fh
  local found 0
  local hits 0
  file open `fh' using "`lf'", read text
  file read `fh' line
  while r(eof) == 0 {
    if strpos(`"`macval(line)'"', `"`needle'"') != 0 {
      local found 1
      local ++hits
    }
    file read `fh' line
  }
  file close `fh'
  return scalar found = `found'
  return scalar hits = `hits'
end

local needle15 Caution: the analytic covariance of the subgroup effects is

* ---- G1 the screen fires, and G6 replay
use data/dta_binary.dta, clear
qui log using sgd_g1.txt, text replace name(sd)
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) esprelength(2) espostlength(2)
local rc1 = _rc
cap noi wooldid
local rc6 = _rc
qui log close sd
local n1 = e(wdid_sgdegen)
local a1 = "`e(wdid_sgdegenact)'"
local r1 = "`e(wdid_sgdegenr)'"
local j1 = e(joint_att_allsgssame)
local jp1 = e(joint_pre_att_allsgssame)
matrix G1P = e(pwsgs)
scanSG sgd_g1.txt "`needle15'"
local c1 = r(found)
* the caution has to appear on the fit AND on the replay
local c1n = r(hits)
* display-eqtest-bare-dot: both equality tests are withheld here, so no "have the same" line
* may print on the fit or on the replay
scanSG sgd_g1.txt "have the same"
local d1 = r(hits)
cap erase sgd_g1.txt
di "G1 BS-15: rc=`rc1' flagged=`n1' act=`a1' rho=`r1' allsgssame=`j1' pre=`jp1' pw se=" G1P[1,2] "/" G1P[2,2] " b=" G1P[1,1] " caution hits=`c1n' (replay rc `rc6')"
if `rc1' != 0 | `n1' != 2 | "`a1'" != "blank" | `j1' < . | `jp1' < . | G1P[1,1] >= . | G1P[1,2] < . | G1P[1,4] < . | G1P[1,6] < . | G1P[2,2] < . | `c1' != 1 {
  di as error "G1 FAIL: the screen did not blank both rows and both equality tests"
  local ++fails
}
if `rc6' != 0 | `c1n' < 2 {
  di as error "G6 FAIL: replay did not reproduce the caution from e() (rc `rc6', `c1n' hit(s))"
  local ++fails
}

* ---- G2 a healthy robust design is untouched
use data/dta_binary.dta, clear
qui log using sgd_g2.txt, text replace name(sd)
cap noi wooldid y i t ttre, att robust subgroups(sg) esprelength(2) espostlength(2)
local rc2 = _rc
qui log close sd
local a2 = "`e(wdid_sgdegenact)'"
local j2 = e(joint_att_allsgssame)
local jp2 = e(joint_pre_att_allsgssame)
matrix G2P = e(pwsgs)
scanSG sgd_g2.txt "`needle15'"
local c2 = r(found)
cap erase sgd_g2.txt
di "G2 BS-15: rc=`rc2' act=[`a2'] caution=`c2' allsgssame=" %10.8f `j2' " pre=" %10.8f `jp2' " pw se=" %10.7f G2P[1,2] "/" %10.7f G2P[2,2]
if `rc2' != 0 | "`a2'" != "" | `c2' != 0 | abs(`j2' - .48310953) > 1e-7 | abs(`jp2' - .43928221) > 1e-7 | abs(G2P[1,2] - .1002352) > 1e-7 | G2P[2,2] >= . {
  di as error "G2 FAIL: the screen fired, or moved a number, on a healthy robust design"
  local ++fails
}

* ---- G3 safety(off) reports rather than withholds
use data/dta_binary.dta, clear
qui log using sgd_g3.txt, text replace name(sd)
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) safety(off) esprelength(2) espostlength(2)
local rc3 = _rc
qui log close sd
local a3 = "`e(wdid_sgdegenact)'"
local j3 = e(joint_att_allsgssame)
matrix G3P = e(pwsgs)
scanSG sgd_g3.txt "`needle15'"
local c3 = r(found)
* display-eqtest-bare-dot control: the reported test still prints its value
scanSG sgd_g3.txt "have the same att: "
local d3 = r(hits)
cap erase sgd_g3.txt
di "G3 BS-15: rc=`rc3' act=[`a3'] caution=`c3' allsgssame=`j3' pw se=" G3P[1,2]
if `rc3' != 0 | "`a3'" != "report" | `c3' != 1 | `j3' >= . | G3P[1,2] >= . {
  di as error "G3 FAIL: safety(off) must report the residue with the caution, not withhold it"
  local ++fails
}

* ---- G4 the two policies side by side
use data/dta_binary.dta, clear
qui log using sgd_g4.txt, text replace name(sd)
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) jointtests esprelength(2) espostlength(2)
local rc4 = _rc
qui log close sd
local jh4 = "`e(joint_horiz)'"
scanSG sgd_g4.txt "`needle15'"
local c4a = r(found)
scanSG sgd_g4.txt "Caution: jointtests was combined with analytic clustered"
local c4b = r(found)
cap erase sgd_g4.txt
di "G4 BS-15: rc=`rc4' bs15 caution=`c4a' jointtests caution=`c4b' joint_horiz=[`jh4']"
if `rc4' != 0 | `c4a' != 1 | `c4b' != 1 | "`jh4'" == "" {
  di as error "G4 FAIL: the BS-15 and jointtests cautions do not stand side by side with the statistic computed"
  local ++fails
}

* ---- G5 no false flag on a healthy continuous-treatment subgroup run
use data/dta_cont.dta, clear
qui gen sg = mod(i, 2) + 1
qui log using sgd_g5.txt, text replace name(sd)
cap noi wooldid y i t ttre, att robust subgroups(sg) contreat(intensity) contreatcontrols(basic) esprelength(2) espostlength(2)
local rc5 = _rc
qui log close sd
* G10 compares its suppressed twin's level rows against these
cap matrix G5S1 = e(wdidsgresults_1)
cap matrix G5S2 = e(wdidsgresults_2)
local a5 = "`e(wdid_sgdegenact)'"
scanSG sgd_g5.txt "`needle15'"
local c5 = r(found)
cap erase sgd_g5.txt
di "G5 BS-15: rc=`rc5' act=[`a5'] caution=`c5'"
if `rc5' != 0 | "`a5'" != "" | `c5' != 0 {
  di as error "G5 FAIL: rc `rc5', act [`a5'], caution `c5' - a false flag (or a broken fixture) on a healthy continuous-treatment subgroup run"
  local ++fails
}

* ---- G7 display-eqtest-bare-dot: the withheld equality tests print nothing (G1, fit and
* replay), the reported ones under safety(off) still print (G3)
di "G7 display-eqtest-bare-dot: withheld have-the-same lines=`d1' (expect 0) reported att lines=`d3' (expect >= 1)"
if `d1' != 0 | `d3' < 1 {
  di as error "G7 FAIL: a withheld equality test printed (`d1' line(s)) or the reported one did not (`d3')"
  local ++fails
}

* ---- G8 the G1 call under suppresssgprimaryeffects: same flag count, same largest implied
* correlation, same blanking, same caution on the fit and on the replay
use data/dta_binary.dta, clear
qui log using sgd_g8.txt, text replace name(sd)
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) esprelength(2) espostlength(2) suppresssgprimaryeffects
local rc8 = _rc
cap noi wooldid
local rc8r = _rc
qui log close sd
local n8 = e(wdid_sgdegen)
local a8 = "`e(wdid_sgdegenact)'"
local r8 = "`e(wdid_sgdegenr)'"
local j8 = e(joint_att_allsgssame)
local jp8 = e(joint_pre_att_allsgssame)
matrix G8P = e(pwsgs)
scanSG sgd_g8.txt "`needle15'"
local c8n = r(hits)
cap erase sgd_g8.txt
di "G8 BS-15 suppress: rc=`rc8' flagged=`n8' act=`a8' rho=`r8' (G1 `r1') allsgssame=`j8' pre=`jp8' pw se=" G8P[1,2] "/" G8P[2,2] " b=" G8P[1,1] " caution hits=`c8n' (replay rc `rc8r')"
if `rc8' != 0 | `rc8r' != 0 | `n8' != 2 | "`a8'" != "blank" | "`r8'" != "`r1'" | `j8' < . | `jp8' < . | G8P[1,1] != G1P[1,1] | G8P[2,1] != G1P[2,1] | G8P[1,2] < . | G8P[1,4] < . | G8P[1,6] < . | G8P[2,2] < . | G8P[2,4] < . | `c8n' < 2 {
  di as error "G8 FAIL: under suppresssgprimaryeffects the screen did not blank both rows and both equality tests as G1 does"
  local ++fails
}

* ---- G9 the continuous-treatment fixture of G5 clustered at i, with and without
* suppresssgprimaryeffects. Clustering at i fires the cell-containment scrub, which blanks the
* ame standard errors, level and contrast alike: the screen skips those rows on both runs and
* flags the att and pre_att rows on both
use data/dta_cont.dta, clear
qui gen sg = mod(i, 2) + 1
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) contreat(intensity) contreatcontrols(basic) esprelength(2) espostlength(2)
local rc9u = _rc
local n9u = e(wdid_sgdegen)
local r9u = "`e(wdid_sgdegenr)'"
qui log using sgd_g9.txt, text replace name(sd)
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) contreat(intensity) contreatcontrols(basic) esprelength(2) espostlength(2) suppresssgprimaryeffects
local rc9 = _rc
qui log close sd
local n9 = e(wdid_sgdegen)
local a9 = "`e(wdid_sgdegenact)'"
local r9 = "`e(wdid_sgdegenr)'"
local j9 = e(joint_att_allsgssame)
matrix G9P = e(pwsgs)
local g9rn : rowfullnames G9P
scanSG sgd_g9.txt "`needle15'"
local c9 = r(found)
cap erase sgd_g9.txt
di "G9 BS-15 contreat suppress: rc=`rc9' (unsuppressed rc `rc9u' flagged `n9u' rho `r9u') flagged=`n9' act=`a9' rho=`r9' allsgssame=`j9' caution=`c9' rows=`g9rn'"
if `rc9u' != 0 | `n9u' != 2 | `rc9' != 0 | `n9' != 2 | "`a9'" != "blank" | "`r9'" != "`r9u'" | `j9' < . | `c9' != 1 | G9P[1,2] < . | G9P[2,2] < . | G9P[3,1] >= . {
  di as error "G9 FAIL: under suppresssgprimaryeffects the continuous-treatment run did not flag exactly the att rows its unsuppressed twin flags"
  local ++fails
}

* ---- G10 the healthy G5 fixture under suppresssgprimaryeffects: no flag, and the hidden level
* rows match the unsuppressed run's level tables. The att rows there come from the same pairwise
* pass; the ame rows from a separate margins call over the same subpopulation, so they are held
* to the 1e-10 relative threshold rather than bit for bit
use data/dta_cont.dta, clear
qui gen sg = mod(i, 2) + 1
qui log using sgd_g10.txt, text replace name(sd)
cap noi wooldid y i t ttre, att robust subgroups(sg) contreat(intensity) contreatcontrols(basic) esprelength(2) espostlength(2) suppresssgprimaryeffects
local rc10 = _rc
qui log close sd
local a10 = "`e(wdid_sgdegenact)'"
scanSG sgd_g10.txt "`needle15'"
local c10 = r(found)
cap erase sgd_g10.txt
local nlev 0
local lvmax 0
local lvbad
cap matrix G10L = e(wdid_sglev)
if _rc == 0 {
  local lvn : rowfullnames G10L
  local kk 0
  foreach wn of local lvn {
    local ++kk
    if !regexm("`wn'", "^([^:]+):([0-9]+)(bn?)?[.]") continue
    local weq = regexs(1)
    local wsg = regexs(2)
    local wit = subinstr("`weq'", "pre_", "", 1)
    local wkey = cond("`wit'" != "`weq'", "pretreatment", "main") + ":`wit'"
    local sn
    cap local sn : rowfullnames G5S`wsg'
    local ss 0
    foreach s of local sn {
      local ++ss
      if "`s'" != "`wkey'" continue
      local ++nlev
      local dd = reldif(G10L[`kk', 2], G5S`wsg'[`ss', 2])
      if `dd' >= . | G10L[`kk', 2] >= . local dd 1
      if `dd' > `lvmax' local lvmax = `dd'
      if `dd' > 1e-10 local lvbad `lvbad' `wn'
    }
  }
}
di "G10 BS-15 contreat suppress robust: rc=`rc10' act=[`a10'] caution=`c10' level rows matched=`nlev' (expect 8) max reldif=`lvmax' off=[`lvbad']"
if `rc10' != 0 | "`a10'" != "" | `c10' != 0 | `nlev' != 8 | "`lvbad'" != "" {
  di as error "G10 FAIL: a flag on a healthy run, or e(wdid_sglev) does not carry the unsuppressed level standard errors for att, pre_att, ame and pre_ame of both subgroups"
  local ++fails
}

if `fails' == 0 di as result "SGDEGEN_ALL_PASS"
else di as error "SGDEGEN_FAILS=`fails'"
