* test_stored.do - the stored-result contract added by the 2026-09-08 rulings
* (batch 3: CS-60, PS-1, CS-54, CS-57, CS-4).
* Run from test/: stata-mp -e do test_stored.do -> expect STORED_ALL_PASS
*
* S1  CS-60: the four subgroup joint p-values are numeric SCALARS, not string macros,
*     and they carry the same values the string macros used to
* S2  PS-1: a clustered run posts e(N_clust), e(clustvar), e(vce)=cluster,
*     e(vcetype)=Robust alongside the existing e(Nclus)
* S3  PS-1: an unclustered OLS run posts e(vce)=ols and no e(clustvar); a poisson run
*     posts e(vce)=robust, e(vcetype)=Robust (ruling 2026-09-25, POIS-1: ppmlhdfe's
*     variance is robust unasked), and so does a pweighted linear run (SEB-05); the
*     header's e(wdid_vcespec) says vce(robust) on both; customvce(oim) still posts oim
* S4  PS-1: a bootstrap run relabels the variance as bootstrap/Bootstrap
* S5  CS-54: e(boot_percentile_mainresults) exists and is element-identical to the
*     legacy e(boot_pctl_mainresults)
* S6  CS-57: e(perm_bandlo_*)/e(perm_bandhi_*) exist for every permuted estimand and
*     equal the e(perm_nullband) rows by name; the deprecated aliases e(perm_lb_*)/
*     e(perm_ub_*) are no longer posted (ruling 2026-09-25)
* S8  PS-62: the within-family event-study covariance blocks are spliced into e(V), every
*     cross-family cell stays zero, and e(wdid_vblocks)/e(wdid_vblockn) name what landed
* S7  CS-4: e(wdid_vcvstate) reports what e(V) actually is - analytic-blocks once an event
*     study exists (PS-62), analytic-diagonal on the
*     analytic path, bootstrap under bootstrap(), absent when no e(V) was posted
* S10 e(histogramestimates): one row per treated cohort i, columns att_i (att) and att_itime
*     (att_it), and the att_i column reaggregates to the att with post-period row counts
clear all
set more off
do "../wooldidinstall/wooldid.ado"
local fails 0
clear
set obs 16
gen unit = _n
gen s = cond(unit<=8, 1, 2)
gen ttre = cond(mod(unit,4)==1, 4, cond(mod(unit,4)==2, 6, .))
expand 8
bysort unit: gen t = _n
set seed 31
gen y = 2*(t>=ttre & ttre<.) + 0.1*unit + 0.05*t + rnormal()
gen byte ycount = max(0, round(y + 5))
tempfile panel
save `panel'

* ---- S1 CS-60: subgroup joint p-values are scalars, not string macros
* the four names are joint_horiz_s*, joint_horiz_i_s*, joint_Choriz_s* and
* joint_Choriz_i_s*; the two Choriz twins need an identified contreat event study and
* are not reachable from a compact fixture, so this pins the two horiz families and the
* absence of ANY joint_*_s* macro, which is what would catch a missed site
clear
set obs 48
gen unit = _n
gen s = cond(mod(unit,2)==0, 1, 2)
gen ttre = cond(mod(unit,4)==0, 5, cond(mod(unit,4)==1, 7, .))
expand 12
bysort unit: gen t = _n
set seed 31
gen dose = runiform() + 0.2*(t>=ttre & ttre<.)
gen y = 2*(t>=ttre & ttre<.) + 0.5*dose + 0.1*unit + 0.05*t + rnormal()
tempfile wide
save `wide'
qui wooldid y unit t ttre, att att_it cluster(unit) subgroups(s) contreat(dose) contreatcontrols(basic) jointtests esprelength(2) espostlength(2)
local sc1 : e(scalars)
local mc1 : e(macros)
local s1bad
local s1seen 0
foreach nm of local mc1 {
  if regexm("`nm'", "^joint_C?horiz(_i)?_s[0-9]+$") local s1bad `s1bad' macro:`nm'
}
foreach nm of local sc1 {
  if regexm("`nm'", "^joint_C?horiz(_i)?_s[0-9]+$") local ++s1seen
}
if `s1seen' < 4 | "`s1bad'" != "" {
  dis as error "S1 FAIL: scalars seen=`s1seen' (want 4 or more); still stored as macros:`s1bad'"
  local ++fails
}
else dis as result "S1 PASS: `s1seen' subgroup joint p-value(s) stored as scalars, none as macros"

* ---- S2 PS-1 clustered
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2)
local nclus_new = e(N_clust)
local nclus_old = e(Nclus)
local cvar = "`e(clustvar)'"
local vce  = "`e(vce)'"
local vct  = "`e(vcetype)'"
if `nclus_new' != `nclus_old' | `nclus_new' >= . | "`cvar'" != "unit" | "`vce'" != "cluster" | "`vct'" != "Robust" {
  dis as error "S2 FAIL: N_clust=`nclus_new' Nclus=`nclus_old' clustvar=`cvar' vce=`vce' vcetype=`vct'"
  local ++fails
}
else dis as result "S2 PASS: clustered run posts N_clust=`nclus_new' clustvar=`cvar' vce=`vce' vcetype=`vct'"

* ---- S3 PS-1 unclustered and poisson
use `panel', clear
qui wooldid y unit t ttre, att esprelength(2) espostlength(2)
local vce3 = "`e(vce)'"
local cvar3 = "`e(clustvar)'"
use `panel', clear
cap qui wooldid ycount unit t ttre, att poisson esprelength(2) espostlength(2)
local rc3p = _rc
local vce3p = "`e(vce)'|`e(vcetype)'|`e(wdid_vcespec)'"
tempname b3p v3p
if `rc3p' == 0 {
  matrix `b3p' = e(b)
  matrix `v3p' = e(V)
}
* the same poisson run with robust typed: the label was the only difference, so the numbers
* must agree bit for bit
use `panel', clear
cap qui wooldid ycount unit t ttre, att poisson robust esprelength(2) espostlength(2)
local rc3r = _rc
local d3r = .
if `rc3p' == 0 & `rc3r' == 0 local d3r = mreldif(`b3p', e(b)) + mreldif(`v3p', e(V))
* pweights force a robust variance on the linear path
use `panel', clear
gen double pw = 1 + mod(unit, 3)
qui wooldid y unit t ttre [pw=pw], att esprelength(2) espostlength(2)
local vce3w = "`e(vce)'|`e(vcetype)'|`e(wdid_vcespec)'"
* customvce(oim) under poisson keeps oim, on the customvce arm
use `panel', clear
cap qui wooldid ycount unit t ttre, att poisson customvce(oim) esprelength(2) espostlength(2)
local rc3o = _rc
local vce3o = "`e(vce)'|`e(vcetype)'"
if "`vce3'" != "ols" | "`cvar3'" != "" | `rc3p' != 0 | "`vce3p'" != "robust|Robust|vce(robust)" | `rc3r' != 0 | `d3r' != 0 | "`vce3w'" != "robust|Robust|vce(robust)" | `rc3o' != 0 | "`vce3o'" != "oim|" {
  dis as error "S3 FAIL: ols run vce=`vce3' clustvar=`cvar3'; poisson rc=`rc3p' vce=`vce3p' (robust-typed rc=`rc3r' reldif=`d3r'); pweight vce=`vce3w'; customvce(oim) rc=`rc3o' vce=`vce3o'"
  local ++fails
}
else dis as result "S3 PASS: unclustered vce=`vce3' (no clustvar); poisson vce=`vce3p' (same numbers as robust typed); pweight vce=`vce3w'; customvce(oim) vce=`vce3o'"

* ---- S4/S5/S7b bootstrap
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) bootstrap(standard) reps(12) rseed(4) esprelength(2) espostlength(2)
local vce4 = "`e(vce)'"
local vct4 = "`e(vcetype)'"
local st4  = "`e(wdid_vcvstate)'"
if "`vce4'" != "bootstrap" | "`vct4'" != "Bootstrap" {
  dis as error "S4 FAIL: bootstrap run vce=`vce4' vcetype=`vct4'"
  local ++fails
}
else dis as result "S4 PASS: bootstrap run vce=`vce4' vcetype=`vct4'"
cap matrix PA = e(boot_percentile_mainresults)
local rc5 = _rc
cap matrix PL = e(boot_pctl_mainresults)
local rc5b = _rc
local d5 = .
if `rc5' == 0 & `rc5b' == 0 local d5 = mreldif(PA, PL)
if `rc5' != 0 | `rc5b' != 0 | `d5' != 0 {
  dis as error "S5 FAIL: alias rc `rc5' legacy rc `rc5b' mreldif `d5'"
  local ++fails
}
else dis as result "S5 PASS: e(boot_percentile_mainresults) matches the legacy name exactly"
if "`st4'" != "bootstrap" {
  dis as error "S7b FAIL: e(wdid_vcvstate) = `st4' under bootstrap()"
  local ++fails
}
else dis as result "S7b PASS: e(wdid_vcvstate) = `st4' under bootstrap()"

* ---- S6 CS-57 permutation band names
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(4) esprelength(2) espostlength(2)
local s6bad
local s6n 0
tempname S6NB
matrix `S6NB' = e(perm_nullband)
foreach e in `e(wdid_permest)' {
  local ++s6n
  local s6c = colnumb(`S6NB', "main_`e'")
  if `s6c' >= . local s6bad `s6bad' nocol_`e'
  else {
    if e(perm_bandlo_`e') != `S6NB'[1, `s6c'] local s6bad `s6bad' bandlo_`e'
    if e(perm_bandhi_`e') != `S6NB'[2, `s6c'] local s6bad `s6bad' bandhi_`e'
  }
  if e(perm_bandlo_`e') >= . local s6bad `s6bad' missing_`e'
  if e(perm_lb_`e') < . | e(perm_ub_`e') < . local s6bad `s6bad' alias_`e'
}
if `s6n' == 0 | "`s6bad'" != "" {
  dis as error "S6 FAIL: estimands=`s6n' mismatches:`s6bad'"
  local ++fails
}
else dis as result "S6 PASS: `s6n' permuted estimand(s) carry band scalars equal to e(perm_nullband), and no deprecated alias"
local st6 = "`e(wdid_vcvstate)'"
local vce6a = "`e(vce)'`e(vcetype)'"

* ---- S7 CS-4 analytic and absent
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2)
local st7 = "`e(wdid_vcvstate)'"
cap confirm matrix e(V)
local hasv7 = _rc == 0
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) fastnose esprelength(2) espostlength(2)
local st7b = "`e(wdid_vcvstate)'"
cap confirm matrix e(V)
local hasv7b = _rc == 0
* a run with no e(V) leaves e(vce)/e(vcetype) empty - they name what that matrix is -
* while e(clustvar) names the design and stands either way
local lab7b = "`e(vce)'`e(vcetype)'"
local cvar7b = "`e(clustvar)'"
* RE-PINNED 2026-09-10 (PS-62): this run has an event study, so its e(V) carries the
* family's genuine covariance blocks and the state reads analytic-blocks. The
* analytic-DIAGONAL state is still exercised - it is what a run with no event study
* leaves - so both words are checked rather than one replacing the other
if "`st7'" != "analytic-blocks" | `hasv7' != 1 | "`st7b'" != "absent" | `hasv7b' != 0 | "`lab7b'" != "" | "`cvar7b'" != "unit" {
  dis as error "S7 FAIL: analytic state=`st7' hasV=`hasv7'; fastnose state=`st7b' hasV=`hasv7b' vce+vcetype=[`lab7b'] clustvar=`cvar7b'"
  local ++fails
}
else dis as result "S7 PASS: analytic=`st7' (V present), fastnose=`st7b' (no V, vce cleared, clustvar=`cvar7b')"
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit)
local st7c = "`e(wdid_vcvstate)'"
local vb7c = "`e(wdid_vblocks)'"
if "`st7c'" != "analytic-diagonal" | "`vb7c'" != "" {
  dis as error "S7d FAIL: a run with no event study must stay analytic-diagonal: state=`st7c' blocks=[`vb7c']"
  local ++fails
}
else dis as result "S7d PASS: no event study -> `st7c', no blocks"
* the permute fast-mode default posts no analytic variance; permanalyticse restores it,
* and the state has to follow e(V) in both directions or the wire is only half-checked
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(4) permanalyticse esprelength(2) espostlength(2)
local st6b = "`e(wdid_vcvstate)'"
cap confirm matrix e(V)
local hasv6b = _rc == 0
local vce6 = "`e(vce)'"
* RE-PINNED 2026-09-10 (PS-62): permute() reposts no covariance of its own, so the
* matrix under permanalyticse is still the worker's spliced one and the label must keep
* saying so rather than flattening back to analytic-diagonal
if "`st6'" != "absent" | "`vce6a'" != "" | "`st6b'" != "analytic-blocks" | `hasv6b' != 1 | "`vce6'" != "cluster" {
  dis as error "S7c FAIL: fast-mode state=`st6' vce+vcetype=[`vce6a']; permanalyticse state=`st6b' hasV=`hasv6b' vce=`vce6'"
  local ++fails
}
else dis as result "S7c PASS: permute fast-mode=`st6' (no V, vce cleared), permanalyticse=`st6b' (V present, vce=`vce6')"

* ---- S8 PS-62: the spliced blocks are the family covariance, and nothing else moved ----
* Every off-diagonal wooldid writes must be the corresponding e(V_es_*) cell exactly; every
* cell BETWEEN two families must still be zero, because those covariances are not estimated
* and a zero there means "absent", not "independent".
use `panel', clear
qui wooldid y unit t ttre, att att_it cluster(unit) esprelength(2) espostlength(2)
local s8st = "`e(wdid_vcvstate)'"
local s8fam = "`e(wdid_vblocks)'"
local s8n = e(wdid_vblockn)
mata:
  V = st_matrix("e(V)")
  vs = st_matrixcolstripe("e(V)")
  fams = tokens(st_local("s8fam"))
  bad = 0
  hit = 0
  for (f = 1; f <= cols(fams); f++) {
    E = st_matrix("e(V_es_" + fams[f] + ")")
    es = st_matrixcolstripe("e(V_es_" + fams[f] + ")")
    px = "ES_" + fams[f] + "_"
    for (j = 1; j <= rows(es); j++) {
      for (i = j + 1; i <= rows(es); i++) {
        a = 0
        b = 0
        for (r = 1; r <= rows(vs); r++) {
          if (vs[r, 2] == px + es[j, 2]) a = r
          if (vs[r, 2] == px + es[i, 2]) b = r
        }
        if (a == 0 | b == 0) continue
        if (V[a, a] <= 0 | V[a, a] >= . | V[b, b] <= 0 | V[b, b] >= .) continue
        if (E[j, i] >= .) continue
        hit = hit + 1
        if (V[a, b] != E[j, i] | V[b, a] != E[j, i]) bad = bad + 1
      }
    }
  }
  ka = 0
  kb = 0
  for (r = 1; r <= rows(vs); r++) {
    if (vs[r, 2] == "ES_att_post0") ka = r
    if (vs[r, 2] == "ES_att_it_post0") kb = r
  }
  st_numscalar("s8hit", hit)
  st_numscalar("s8bad", bad)
  st_numscalar("s8cross", (ka > 0 & kb > 0 ? V[ka, kb] : .))
  st_numscalar("s8found", (ka > 0 & kb > 0))
end
di "S8 state=`s8st' families=[`s8fam'] pairs=`s8n' checked=" s8hit " mismatched=" s8bad " cross-family cell=" s8cross
if "`s8st'" != "analytic-blocks" | "`s8fam'" == "" | s8hit < 1 | s8bad != 0 | `s8n' != s8hit | s8found != 1 | s8cross != 0 {
  dis as error "S8 FAIL: state=`s8st' families=[`s8fam'] posted pairs=`s8n' checked=" s8hit " mismatched=" s8bad " cross=" s8cross
  local ++fails
}
else dis as result "S8 PASS: `s8n' off-diagonal pair(s) across [`s8fam'] equal e(V_es_*) exactly; the att-vs-att_it cell is still 0"

* ---- S8b PS-62 consumer side: what the splice is FOR - a cross-horizon lincom on e(V) ----
* the se of post0 - post1 from lincom (which reads e(V) and e(df_r)) must equal the family
* covariance's own sqrt(v00 + v11 - 2 v01); before PS-62 lincom silently used the zero
* off-diagonal and reported sqrt(v00 + v11)
use `panel', clear
qui wooldid y unit t ttre, att cluster(unit) esprelength(2) espostlength(2)
cap noi lincom ES_att_post0 - ES_att_post1
local s8brc = _rc
local s8bse = r(se)
local s8best = r(estimate)
mata: VE = st_matrix("e(V_es_att)"); nm = st_matrixcolstripe("e(V_es_att)")[., 2]; a = selectindex(nm :== "post0"); b = selectindex(nm :== "post1"); BE = st_matrix("e(b_es_att)"); st_numscalar("s8b_se", sqrt(VE[a,a] + VE[b,b] - 2*VE[a,b])); st_numscalar("s8b_est", BE[1,a] - BE[1,b]); st_numscalar("s8b_diag", sqrt(VE[a,a] + VE[b,b]))
if `s8brc' != 0 | reldif(`s8bse', s8b_se) > 1e-12 | reldif(`s8best', s8b_est) > 1e-12 | reldif(s8b_se, s8b_diag) < 1e-6 {
  dis as error "S8b FAIL: lincom rc=`s8brc' se=`s8bse' vs family sqrt(v00+v11-2v01)=" s8b_se " (diagonal-only would give " s8b_diag ")"
  local ++fails
}
else dis as result "S8b PASS: lincom post0 - post1 on the spliced e(V) reproduces the family covariance (se " %10.8f `s8bse' ", diagonal-only would be " %10.8f s8b_diag ")"

* ---- S10 (change audit 2026-09-22, cleanmatrices-cert-gap): histogramcohorteffects' stored
* result, which no cert suite read. i = unit here, so the treated cohorts are the 8 treated units
* (onset 4 or 6 on t = 1..8). att alone gives one column, att_i; att att_it gives att_i and
* att_itime. the rows come in ascending i, so the att_i column, weighted by each unit's
* post-period row count (5 for onset 4, 3 for onset 6), reaggregates to e(b)'s att ----
use `panel', clear
cap noi wooldid y unit t ttre, att histogramcohorteffects
local rc10a = _rc
local r10a = .
local c10a
local id10 = .
if `rc10a' == 0 {
  tempname H10 B10
  matrix `H10' = e(histogramestimates)
  local r10a = rowsof(`H10')
  local c10a : colnames `H10'
  matrix `B10' = e(b)
  local att10 = `B10'[1,1]
  * treated units in ascending order: unit mod 4 == 1 has onset 4 (5 post rows), == 2 onset 6 (3)
  local num 0
  local den 0
  local k 0
  forval u = 1/16 {
    if mod(`u', 4) == 1 | mod(`u', 4) == 2 {
      local ++k
      local n = cond(mod(`u', 4) == 1, 5, 3)
      local num = `num' + `n'*`H10'[`k', 1]
      local den = `den' + `n'
    }
  }
  local id10 = reldif(`num'/`den', `att10')
}
use `panel', clear
cap noi wooldid y unit t ttre, att att_it histogramcohorteffects
local rc10b = _rc
local r10b = .
local c10b
if `rc10b' == 0 {
  tempname H10b
  matrix `H10b' = e(histogramestimates)
  local r10b = rowsof(`H10b')
  local c10b : colnames `H10b'
}
dis "S10: att rc=`rc10a' rows=`r10a' cols=[`c10a'] reaggregation reldif=`id10' | att att_it rc=`rc10b' rows=`r10b' cols=[`c10b']"
if `rc10a' != 0 | `r10a' != 8 | "`c10a'" != "att_i" | `id10' > 1e-10 | `rc10b' != 0 | `r10b' != 8 | "`c10b'" != "att_i att_itime" {
  dis as error "S10 FAIL: e(histogramestimates) shape, columns or reaggregation"
  local ++fails
}
else dis as result "S10 PASS: e(histogramestimates) has one row per treated cohort, the requested columns, and reaggregates to att"

* ---- S9 PS-62 follow-through: cleanmatrices must not take the splice down with it ----
* cleanmatrices runs `clear mata' inside the worker, which drops wooldid's OWN Mata kernels
* for the rest of the session whenever the ado was loaded with `do <path>' - which is how
* this file loads it, and how the suite's dual-build shards load theirs. (An adopath-resolved
* ado recovers on its own, because Stata re-reads the file when a Mata function is missing;
* that asymmetry is why the refusal this combination once carried was removed - SI-17/SI-18.)
* PS-62 put a Mata call on the mainline ANALYTIC path, where none had been, so before the
* guard this exact call died r(3499) instead of estimating - the 501-config assessment caught
* it, no cert leg did. The contract: the run completes and lands on the pre-PS-62 surface,
* an e(V) with no blocks, rather than a half-spliced one or an error.
* KEEP THIS LAST: it clears wooldid's Mata for everything after it.
use `panel', clear
cap qui wooldid y unit t ttre, att cluster(unit) cleanmatrices esprelength(2) espostlength(2)
local s9rc = _rc
local s9st = "`e(wdid_vcvstate)'"
local s9vb = "`e(wdid_vblocks)'"
cap confirm matrix e(V)
local s9hasv = _rc == 0
if `s9rc' != 0 | `s9hasv' != 1 | "`s9st'" != "analytic-diagonal" | "`s9vb'" != "" {
  dis as error "S9 FAIL: cleanmatrices rc=`s9rc' hasV=`s9hasv' state=`s9st' blocks=[`s9vb']"
  local ++fails
}
else dis as result "S9 PASS: cleanmatrices completes (rc `s9rc'), e(V) present, state=`s9st', no blocks"

if `fails' == 0 dis as result "STORED_ALL_PASS"
else dis as error "STORED_FAILS=`fails'"
