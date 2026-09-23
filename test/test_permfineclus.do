* Targeted suite for T1 (2026-09-16): permute() with cluster() FINER than i, under the
* user's assertion permonsetcohorts that i IS the treatment-onset grouping. Each
* replication then rebuilds i from the onsets it drew, so ttre is single-valued within the
* rebuilt i by construction and a unit can no longer straddle two onsets in one draw.
* Covered here: the parse-time verification of the assertion (both directions of the
* bijection), the subgroups() and stray-companion refusals, that the refusal of finer
* clustering STILL stands without the option, the declared-no-op property at cluster() = i
* (bitwise), the column-set/denominator integrity of the accumulator, and that the exact
* test's observed-is-a-member property survives the rebuild path, and that the assertion reads
* an onset after the last period as never-treated under a finer cluster() (Q7) and as given at
* cluster() = i (Q8). What is NOT here, said so
* the file does not imply coverage it has not got: the check that the onset-to-cohort map is
* right. That one lives in wooldid_permrun (rep 0's rebuilt i must equal the observed i) and
* is measured to fire; Q5 below says why nothing in this file can do it instead.
* Run from test/:  "<stata>" -e do test_permfineclus.do   -> expect PERMFINECLUS_ALL_PASS
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
set seed 13579
set sortseed 13579
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

cap program drop xok
program define xok
  args label cond
  if `cond' != 1 {
    dis as error "`label' FAIL"
    global pf = $pf + 1
  }
  else dis as result "`label' PASS"
end

* fixture: nclus units in 3 onset cohorts (onset 5, onset 7, never-treated), np periods.
*   i     = the unit, the level at which onset is randomized  -> cluster()
*   icoh  = the onset grouping, one level per distinct ttre   -> i
*   icoh2 = a relabelling of icoh (same partition, other codes)
*   icohbad = icoh with the never-treated split in two (breaks the bijection)
*   iunit = i itself, FINER than the onset grouping (breaks the bijection the other way)
*   sg    = a unit-level subgroup
cap program drop mkfc
program define mkfc
  args nclus ntr np eff
  clear
  set seed 4321
  set sortseed 4321
  qui set obs `nclus'
  qui gen i = _n
  qui gen u = runiform()
  sort u
  qui gen ttre = cond(_n <= `ntr', 5, cond(_n <= 2*`ntr', 7, .))
  sort i
  qui gen icoh = cond(ttre == 5, 1, cond(ttre == 7, 2, 3))
  qui gen icoh2 = icoh*7 + 3
  qui gen icohbad = cond(icoh == 3 & mod(i, 2) == 0, 4, icoh)
  qui gen sg = mod(i, 2)
  qui gen a = rnormal()
  qui expand `np'
  qui bysort i: gen t = _n
  qui gen y = a + .1*t + rnormal()
  if `eff' != 0 qui replace y = y + `eff'*(t >= ttre & ttre < .)
  drop a u
  sort i t
end

local ESO esprelength(2) espostlength(2)
local PO permute(shuffle) permreps(99) rseed(4242) dots(0)

*=============================================================== Q1 the assertion is VERIFIED
mkfc 12 3 10 0

* Q1a: i finer than the onset grouping (the unit itself) - onset constant within i holds,
* but several i share an onset, so the map has no unique image
cap noi qui wooldid y i t ttre, att cluster(i) `ESO' `PO' permonsetcohorts
xrc "Q1a i finer than the onset grouping" 459 `=_rc'

* Q1b: the never-treated split across two levels of i - the direction a user who defines
* i as cohort-by-region trips
cap noi qui wooldid y icohbad t ttre, att cluster(icohbad) `ESO' `PO' permonsetcohorts
xrc "Q1b never-treated split across levels of i" 459 `=_rc'

* Q1c: the honest case passes parse
cap noi qui wooldid y icoh t ttre, att cluster(i) `ESO' `PO' permonsetcohorts
xrc "Q1c i IS the onset grouping, cluster finer" 0 `=_rc'

* Q1d: safety(off) does NOT buy past the verification - it is a mechanical precondition
cap noi qui wooldid y i t ttre, att cluster(i) `ESO' `PO' permonsetcohorts safety(off)
xrc "Q1d safety(off) does not bypass the verification" 459 `=_rc'

* Q1e: a design with NO never-treated units, which is where the second direction of the
* verification has no missing onset group to check. wooldid refuses such a design outright -
* a documented design rejection, not a degeneracy (test/suite/data_classes.csv on
* dta_path_nonever: "a design wooldid intentionally rejects (r459) ... not to be run to
* results") - and this option must not change WHICH refusal fires. The screens added here run
* before the data setup that rejects it, so an over-broad one would replace an informative
* message with its own; assert the same return code with the option, without it, and with a
* finer cluster. MEASURED 2026-09-16: all three are r459 carrying the data-setup message.
* This is also why the rebuild's reserved never-treated code (the branch taken when the
* design sample holds no never-treated rows) is unreachable in practice: the observed run
* never reaches the replication loop on such a design.
mkfc 12 6 10 0
cap noi qui wooldid y icoh t ttre, att cluster(icoh) `ESO' `PO'
local rce1 = _rc
cap noi qui wooldid y icoh t ttre, att cluster(icoh) `ESO' `PO' permonsetcohorts
local rce2 = _rc
cap noi qui wooldid y icoh t ttre, att cluster(i) `ESO' `PO' permonsetcohorts
local rce3 = _rc
dis "Q1e no never-treated: rc without=`rce1' with=`rce2' with+finer=`rce3'"
xok "Q1e the option does not pre-empt the no-never-treated design rejection" `=(`rce1' == 459 & `rce2' == 459 & `rce3' == 459)'
mkfc 12 3 10 0

*=============================================================== Q2 the refusals
* Q2a: subgroups() is refused with the option (F13 would fire r411 in every replication)
cap noi qui wooldid y icoh t ttre, att cluster(i) subgroups(sg) `ESO' `PO' permonsetcohorts
xrc "Q2a subgroups() refused under the option" 198 `=_rc'

* Q2b: the option without permute() is a stray companion
cap noi qui wooldid y icoh t ttre, att cluster(i) permonsetcohorts
xrc "Q2b option without permute()" 198 `=_rc'

* Q2c: finer clustering WITHOUT the option is still refused - the relaxation is gated on
* the user's assertion and on nothing else
cap noi qui wooldid y icoh t ttre, att cluster(i) `ESO' `PO'
xrc "Q2c finer clustering still refused without the option" 459 `=_rc'

*=============================================================== Q3 the declared no-op (P3)
* cluster() at i: the straddling problem does not arise, so the rebuild is skipped and the
* run must be the pre-T1 path. Compared BITWISE across the whole reported surface.
tempname B0 P0 C0 D0 K0 B1 P1 C1 D1 K1
cap noi qui wooldid y icoh t ttre, att cluster(icoh) `ESO' `PO'
xrc "Q3 setup (option absent)" 0 `=_rc'
local r3a = _rc
matrix `B0' = e(b)
matrix `P0' = e(perm_p)
matrix `C0' = e(perm_ci)
matrix `K0' = e(perm_kernel)
matrix `D0' = e(permdist_b)
local irb0 = e(perm_irebuilt)

cap noi qui wooldid y icoh t ttre, att cluster(icoh) `ESO' `PO' permonsetcohorts
xrc "Q3 setup (option present)" 0 `=_rc'
local r3b = _rc
matrix `B1' = e(b)
matrix `P1' = e(perm_p)
matrix `C1' = e(perm_ci)
matrix `K1' = e(perm_kernel)
matrix `D1' = e(permdist_b)
local irb1 = e(perm_irebuilt)

if `r3a' == 0 & `r3b' == 0 {
  * every number read out of Mata: a local's decimal round-trip would hide a last-digit move
  local bad 0
  foreach m in B P C K D {
    mata: st_local("eq", strofreal(st_matrix("``m'0'") == st_matrix("``m'1'")))
    if "`eq'" != "1" {
      dis as error "Q3 e(`m') differs between option-absent and option-present"
      local bad 1
    }
  }
  local nmk : colfullnames `P0'
  local nmk1 : colfullnames `P1'
  if "`nmk'" != "`nmk1'" local bad 1
  xok "Q3a P3 negative control: cluster() at i, option present == option absent, bitwise" `=(`bad' == 0)'
  xok "Q3b e(perm_irebuilt) is 0 where the rebuild is a no-op" `=("`irb1'" == "0")'
  xok "Q3c e(perm_irebuilt) is absent without the option" `=("`irb0'" == "." | "`irb0'" == "")'
}
else {
  dis as error "Q3 FAIL (VACUOUS): a setup run failed, so nothing was compared"
  global pf = $pf + 1
}

*=============================================================== Q4 accumulator integrity (P2)
cap noi qui wooldid y icoh t ttre, att cluster(i) `ESO' `PO' permonsetcohorts
xrc "Q4 setup" 0 `=_rc'
if _rc == 0 {
  tempname PN4 PP4
  matrix `PN4' = e(perm_n)
  matrix `PP4' = e(perm_p)
  local K4 = colsof(`PN4')
  local bad 0
  local nlive 0
  forval k = 1/`K4' {
    if `PN4'[1, `k'] != 99 local ++bad
    if `PP4'[1, `k'] < . {
      local ++nlive
      if `PP4'[1, `k'] < 0 | `PP4'[1, `k'] > 1 local ++bad
    }
  }
  dis "Q4 columns=`K4' live p=`nlive' bad=`bad' completed=" e(perm_completed) " failed=" e(perm_failed) " blank=" e(perm_blank)
  xok "Q4a e(perm_n) equals permreps on EVERY accumulated column" `=(`bad' == 0 & `K4' > 0)'
  xok "Q4b every reported p is a live number in [0,1]" `=(`nlive' == `K4')'
  xok "Q4c no replication failed and no draw was blanked" `=(e(perm_completed) == 99 & e(perm_failed) == 0 & e(perm_blank) == 0)'
  xok "Q4d the rebuild actually ran" `=(e(perm_irebuilt) == 1)'
}
else {
  dis as error "Q4 FAIL (VACUOUS): the run failed"
  global pf = $pf + 1
}

*=============================================================== Q5 observed-is-a-member
* survives the rebuild path. permexact enumerates every assignment, so one enumerated column
* IS the observed one; under the rebuild its replication runs on a regrouped i and must still
* reproduce the comparator BITWISE, which is the property the exact test's p-values rest on.
* WHAT THIS DOES NOT CHECK, stated so the file does not imply coverage it has not got: it is
* not the check that catches a wrong onset-to-cohort map. Both sides of the comparison are
* built by the same map from the same onsets, so they agree whatever the map says. The map
* itself is checked in wooldid_permrun, by the assertion that rep 0's rebuilt i equals the
* observed i on the design sample - measured 2026-09-16 to fire (r459) on a deliberately
* reversed map, before any enumeration runs. The rep-0 colfullnames tripwire cannot do it
* either: no accumulated e(b) column name embeds an i code, so a wrong map reproduces the
* column set exactly and changes only the numbers.
mkfc 6 2 10 0
cap noi qui wooldid y icoh t ttre, att cluster(i) `ESO' permute(shuffle) permexact rseed(31) dots(0) permonsetcohorts
xrc "Q5 setup (permexact, cluster finer than i)" 0 `=_rc'
if _rc == 0 {
  tempname D5 CB5
  matrix `D5' = e(permdist_b)
  matrix `CB5' = e(perm_comparator_b)
  mata: st_local("nhit", strofreal(sum(rowsum(st_matrix("`D5'") :== J(rows(st_matrix("`D5'")), 1, 1) # st_matrix("`CB5'")) :== cols(st_matrix("`CB5'")))))
  dis "Q5 enumerated=" e(perm_reps) " completed=" e(perm_completed) " rows matching the comparator bitwise=`nhit'"
  xok "Q5a the observed assignment is a bitwise member of its own reference set" `=(`nhit' >= 1)'
  xok "Q5b the enumeration completed" `=(e(perm_completed) == e(perm_reps) & e(perm_failed) == 0)'
}
else {
  dis as error "Q5 FAIL (VACUOUS): the permexact run failed"
  global pf = $pf + 1
}

*=============================================================== Q6 the companions compose
mkfc 12 3 10 0
qui gen str = 1 + (i > 6)

* Q6a: permstrata() - the generator is constrained, the rebuild is downstream of it
cap noi qui wooldid y icoh t ttre, att cluster(i) `ESO' permute(shuffle) permreps(99) ///
    permstrata(str) rseed(77) dots(0) permonsetcohorts
xrc "Q6a permstrata() composes" 0 `=_rc'
if _rc == 0 {
  tempname PN6
  matrix `PN6' = e(perm_n)
  local bad 0
  forval k = 1/`=colsof(`PN6')' {
    if `PN6'[1, `k'] != 99 local ++bad
  }
  xok "Q6a2 permstrata(): every column keeps its full denominator" `=(`bad' == 0)'
}

* Q6b: permtstat - the side benefit. At cluster() = i the help cautions that the
* studentized statistic is unprotected because the cluster-robust se cannot see the treated
* cells' variance; a finer cluster is what makes it informative, so it must at least run
cap noi qui wooldid y icoh t ttre, att cluster(i) `ESO' `PO' permonsetcohorts permtstat
xrc "Q6b permtstat with cluster() below i" 0 `=_rc'
if _rc == 0 xok "Q6b2 permtstat completed every replication" `=(e(perm_completed) == 99 & e(perm_failed) == 0)'

* Q6c: permute(probabilities) - the one method whose drawn onset LEVEL SET can be a strict
* subset of the observed one, so the rebuilt i can have fewer levels in a replication. The
* accumulator must still keep every column: the level set of i moves the per-replication
* ESTIMAND (att_i averages over the cohorts the draw produced), not the column names
local PP (ttre 5 7 .)
forval c = 1/12 {
  local PP `PP' (`c' .25 .25 .5)
}
cap noi qui wooldid y icoh t ttre, att att_i cluster(i) `ESO' permute(probabilities) ///
    permreps(100) permprobs(`PP') rseed(9) dots(0) permonsetcohorts
xrc "Q6c permute(probabilities) composes" 0 `=_rc'
if _rc == 0 {
  tempname PN6c
  matrix `PN6c' = e(perm_n)
  local bad 0
  forval k = 1/`=colsof(`PN6c')' {
    if `PN6c'[1, `k'] != 100 local ++bad
  }
  dis "Q6c columns=" colsof(`PN6c') " short=`bad' completed=" e(perm_completed) " blank=" e(perm_blank)
  xok "Q6c2 probabilities: a shrinking level set of i costs no accumulator column" `=(`bad' == 0)'
}

*=============================================================== Q7 onsets after the last period
* The permutation design reads an onset after the last period as never-treated
* (test_permute.do P13), and under a cluster() finer than i, where the rebuild keys on that
* design onset, the assertion is checked on it, so i must pool the late adopters with the
* never-treated. Fixture: 16 units x t = 1..10, units 1-4 first
* treated at 6, 5-8 at 12, 9-16 never treated; cluster() is the unit, finer than i.
*   icoha = 1 for onset 6, 2 for the late adopters and the never-treated together
*   icohb = 1 for onset 6, 2 for the late adopters, 3 for the never-treated
cap program drop mkfclate
program define mkfclate
  clear
  set seed 2468
  set sortseed 2468
  qui set obs 16
  qui gen i = _n
  qui gen double ttre = cond(_n <= 4, 6, cond(_n <= 8, 12, .))
  qui gen icoha = cond(ttre == 6, 1, 2)
  qui gen icohb = cond(ttre == 6, 1, cond(ttre == 12, 2, 3))
  qui expand 10
  qui bysort i: gen t = _n
  qui gen y = 1.5*(t >= ttre & ttre < .) + 0.1*i + rnormal()
end

* Q7a: icoha is the design's onset grouping, so it runs, the rebuild engages, and the count is
* the design's 16!/(4!12!) = 1,820
mkfclate
cap noi qui wooldid y icoha t ttre, att cluster(i) permute(shuffle) permreps(20) rseed(71) dots(0) permonsetcohorts
local rc7 = _rc
xrc "Q7a late adopters pooled with the never-treated in i" 0 `rc7'
if `rc7' == 0 xok "Q7a2 rebuild engaged over the 1,820 design assignments" `=(e(perm_irebuilt) == 1 & e(perm_npossible) == 1820 & e(perm_completed) == 20)'

* Q7b: icohb splits the design's never-treated level across two levels of i: refused
mkfclate
cap noi qui wooldid y icohb t ttre, att cluster(i) permute(shuffle) permreps(20) rseed(71) dots(0) permonsetcohorts
xrc "Q7b late adopters in their own level of i" 459 `=_rc'

*=============================================================== Q8 onsets after the last period at cluster() = i
* At cluster() = i nothing is rebuilt and the option changes no number, so the assertion is
* checked on the onsets as given there: a late adopter is a level of its own, as at any other
* onset. Q7's fixture with cluster() at i, so each cluster is one cohort.
tempname B8a P8a D8a B8b P8b D8b
* Q8a: icohb gives the late adopters their own level of i: it runs, nothing is rebuilt, and
* the run is bitwise the run without the option
mkfclate
cap noi qui wooldid y icohb t ttre, att cluster(icohb) permute(shuffle) permreps(20) rseed(71) dots(0)
local rc8a = _rc
xrc "Q8a setup (option absent)" 0 `rc8a'
if `rc8a' == 0 {
  matrix `B8a' = e(b)
  matrix `P8a' = e(perm_p)
  matrix `D8a' = e(permdist_b)
}
mkfclate
cap noi qui wooldid y icohb t ttre, att cluster(icohb) permute(shuffle) permreps(20) rseed(71) dots(0) permonsetcohorts
local rc8b = _rc
xrc "Q8a late adopters in their own level of i at cluster() = i" 0 `rc8b'
if `rc8a' == 0 & `rc8b' == 0 {
  matrix `B8b' = e(b)
  matrix `P8b' = e(perm_p)
  matrix `D8b' = e(permdist_b)
  local bad 0
  foreach m in B P D {
    mata: st_local("eq", strofreal(st_matrix("``m'8a'") == st_matrix("``m'8b'")))
    if "`eq'" != "1" local bad 1
  }
  xok "Q8a2 nothing rebuilt, and option present == option absent, bitwise" `=(`bad' == 0 & e(perm_irebuilt) == 0)'
}

* Q8b: icoha pools the late adopters with the never-treated in one level of i, which the onsets
* as given split in two: refused at cluster() = i by the assertion's check, while the same
* run without the option runs
mkfclate
cap noi qui wooldid y icoha t ttre, att cluster(icoha) permute(shuffle) permreps(20) rseed(71) dots(0)
xrc "Q8b setup (option absent)" 0 `=_rc'
mkfclate
qui log using pfc_q8b.txt, text replace name(q8)
cap noi wooldid y icoha t ttre, att cluster(icoha) permute(shuffle) permreps(20) rseed(71) dots(0) permonsetcohorts
local rc8 = _rc
qui log close q8
xrc "Q8b late adopters pooled with the never-treated in i at cluster() = i" 459 `rc8'
xok "Q8b2 the refusal is the assertion's" `=(strpos(fileread("pfc_q8b.txt"), "Error: permonsetcohorts requires that icoha represent treatment-onset groups") > 0)'
cap erase pfc_q8b.txt

if $pf == 0 dis as result "PERMFINECLUS_ALL_PASS"
else dis as error "PERMFINECLUS_FAILS=$pf"
