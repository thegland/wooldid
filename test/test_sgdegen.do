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
* G9 the continuous-treatment fixture clustered at i: the att rows are flagged, and the ame
*    rows, whose standard errors the cell-containment scrub blanks, are skipped rather than
*    counted
* (G8 and G10 pinned the suppresssgprimaryeffects route, retired 2026-09-25; numbers kept)
* G11 under permute() with permanalyticse the screen discloses without blanking (the author's
*    ruling D4, 2026-09-23): act=perm, the caution prints on the fit and on the replay, the
*    se_analytic column is the same residue the safety(off) run reports (G3), and the
*    equality tests still post; under the fast-mode default nothing is flagged
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

* ---- G9 the continuous-treatment fixture of G5 clustered at i. Clustering at i fires the
* cell-containment scrub, which blanks the ame standard errors, level and contrast alike: the
* screen skips those rows and flags the att and pre_att rows
use data/dta_cont.dta, clear
qui gen sg = mod(i, 2) + 1
qui log using sgd_g9.txt, text replace name(sd)
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) contreat(intensity) contreatcontrols(basic) esprelength(2) espostlength(2)
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
di "G9 BS-15 contreat: rc=`rc9' flagged=`n9' act=`a9' rho=`r9' allsgssame=`j9' caution=`c9' rows=`g9rn'"
if `rc9' != 0 | `n9' != 2 | "`a9'" != "blank" | `j9' < . | `c9' != 1 | G9P[1,2] < . | G9P[2,2] < . | G9P[3,1] >= . {
  di as error "G9 FAIL: the continuous-treatment run clustered at i did not flag exactly the att rows"
  local ++fails
}

* ---- G11 permute(): disclose, never blank
use data/dta_binary.dta, clear
qui log using sgd_g11.txt, text replace name(sd)
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) esprelength(2) espostlength(2) permute(shuffle) permreps(20) permanalyticse rseed(11) dots(0)
local rc11 = _rc
cap noi wooldid
local rc11r = _rc
qui log close sd
local a11 = "`e(wdid_sgdegenact)'"
local n11 = e(wdid_sgdegen)
local j11 = e(joint_att_allsgssame)
matrix G11P = e(pwsgs)
scanSG sgd_g11.txt "`needle15'"
local c11 = r(hits)
scanSG sgd_g11.txt "safety(off) was specified"
local s11 = r(hits)
cap erase sgd_g11.txt
local d11 = reldif(G11P[1,2], G3P[1,2])
use data/dta_binary.dta, clear
qui log using sgd_g11b.txt, text replace name(sd)
cap noi wooldid y i t ttre, att cluster(i) subgroups(sg) esprelength(2) espostlength(2) permute(shuffle) permreps(20) rseed(11) dots(0)
local rc11b = _rc
qui log close sd
local a11b = "`e(wdid_sgdegenact)'"
scanSG sgd_g11b.txt "`needle15'"
local c11b = r(found)
cap erase sgd_g11b.txt
di "G11 D4 permute: rc=`rc11' act=[`a11'] flagged=`n11' caution hits=`c11' (replay rc `rc11r') safety-text hits=`s11' allsgssame=`j11' se_analytic vs G3 reldif=`d11'; fast default rc=`rc11b' act=[`a11b'] caution=`c11b'"
if `rc11' != 0 | "`a11'" != "perm" | `n11' != 2 | `c11' < 2 | `rc11r' != 0 | `s11' != 0 | `j11' >= . | `d11' > 1e-10 | G11P[1,1] >= . {
  di as error "G11 FAIL: under permute() the screen must disclose on the fit and the replay without blanking or naming safety(off)"
  local ++fails
}
if `rc11b' != 0 | "`a11b'" != "" | `c11b' != 0 {
  di as error "G11 FAIL: the fast-mode default flagged a row it has no standard error for"
  local ++fails
}

if `fails' == 0 di as result "SGDEGEN_ALL_PASS"
else di as error "SGDEGEN_FAILS=`fails'"
