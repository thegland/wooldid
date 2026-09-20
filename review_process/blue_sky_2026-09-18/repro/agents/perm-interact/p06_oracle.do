*==============================================================================
* p06: EXACT-ENUMERATION ORACLE for the permute() replication model.
* 6 clusters, onset multiset {2003,2003,2005,.,.,.} -> exactly 60 distinct
* shuffle assignments.  wooldid's own e(permdist_b) under permexact is compared
* against a null distribution I enumerate and fit BY HAND with the same options.
* If an option does not reach the replications, the hand-built distribution with
* that option will NOT match wooldid's.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 4242
set sortseed 4242

clear
set obs 6
gen long unit = _n
gen double ttre = .
replace ttre = 2003 if inlist(unit,1,2)
replace ttre = 2005 if unit==3
expand 7
bysort unit: gen int t = 1999 + _n
gen double a_i = rnormal()
gen double g_t = 0.3*(t-2002)
gen double zbig = rnormal() + 0.5*unit
gen double y = a_i + g_t + 1.0*zbig + 2.0*(ttre<. & t>=ttre) + rnormal()*0.25
save "p06_data.dta", replace

* ---- enumerate the 60 distinct assignments ----
mata:
  base = (1 \ 1 \ 2 \ 0 \ 0 \ 0)
  info = cvpermutesetup(base)
  out = J(0,1,"")
  while (1) {
    p = cvpermute(info)
    if (rows(p)==0) break
    s = ""
    for (j=1; j<=6; j++) {
      if (p[j]==1) s = s + "2003|"
      else if (p[j]==2) s = s + "2005|"
      else s = s + ".|"
    }
    out = out \ s
  }
  out = uniqrows(out)
  st_local("nasg", strofreal(rows(out)))
  st_local("asglist", invtokens(out'))
end
dis "P06_NASG=" `nasg'

cap prog drop oracle
prog def oracle
  syntax anything(name=tag) , ASGLIST(string asis) [OPTS(string asis)]
  local vals
  local nfail 0
  foreach s of local asglist {
    preserve
    qui use "p06_data.dta", clear
    qui gen double ttre2 = .
    local u 0
    foreach tok in `=subinstr("`s'","|"," ",.)' {
      local ++u
      if "`tok'" != "." qui replace ttre2 = `tok' if unit == `u'
    }
    cap qui wooldid y unit t ttre2 , att cluster(unit) fastnose `opts'
    if _rc == 0 local vals `vals' `=string(_b[main_att], "%21x")'
    else local ++nfail
    restore
  }
  dis "ORACLE|`tag'|N=" wordcount("`vals'") "|nfail=`nfail'"
  c_local oraclevals `vals'
end

mata:
void wd_cmpdist(string scalar tag, real scalar c, string scalar ovals)
{
    real colvector W, O
    string rowvector tk
    real scalar k, mx
    W = st_matrix("PDSNAP")[., c]
    tk = tokens(ovals)
    O = J(cols(tk), 1, .)
    for (k=1; k<=cols(tk); k++) O[k] = strtoreal(tk[k])
    W = select(W, W :< .)
    O = select(O, O :< .)
    printf("WDCMP|%s|nwooldid=%f|noracle=%f\n", tag, rows(W), rows(O))
    if (rows(W) == rows(O) & rows(W) > 0) {
        W = sort(W, 1)
        O = sort(O, 1)
        mx = max(abs(W - O) :/ (abs(O) :+ 1e-30))
        printf("WDCMP|%s|MAXRELDIFF=%21x\n", tag, mx)
    }
}
end
cap prog drop wdcmp
prog def wdcmp
  syntax anything(name=tag) , VALS(string asis)
  local c = colnumb(PDSNAP, "main_att")
  if `c' >= . {
    dis "WDCMP|`tag'|NO main_att COLUMN; colnames = `: colnames PDSNAP'"
    exit
  }
  mata: wd_cmpdist("`tag'", `c', "`vals'")
end

* ---------- CASE 1: no extra options ----------
use "p06_data.dta", clear
cap noi wooldid y unit t ttre , att cluster(unit) permute(shuffle) permexact dots(0)
dis "C1_RC=" _rc "|npossible=" e(perm_npossible) "|done=" e(perm_completed) "|fail=" e(perm_failed)
cap matrix drop PDSNAP
matrix PDSNAP = e(permdist_b)
oracle plain , asglist(`asglist')
local ov_plain `oraclevals'
wdcmp plain , vals(`ov_plain')

* ---------- CASE 2: controls(zbig) ----------
use "p06_data.dta", clear
cap noi wooldid y unit t ttre , att cluster(unit) permute(shuffle) permexact dots(0) controls(zbig)
dis "C2_RC=" _rc "|npossible=" e(perm_npossible) "|done=" e(perm_completed) "|fail=" e(perm_failed)
cap matrix drop PDSNAP
matrix PDSNAP = e(permdist_b)
oracle ctrl , asglist(`asglist') opts(controls(zbig))
local ov_ctrl `oraclevals'
wdcmp ctrl , vals(`ov_ctrl')
* NEGATIVE CONTROL: must NOT match the no-controls oracle
wdcmp ctrl_vs_plainoracle , vals(`ov_plain')

* ---------- CASE 3: gteffects (not on the allow-list; must be a true no-op) ----------
use "p06_data.dta", clear
cap noi wooldid y unit t ttre , att cluster(unit) permute(shuffle) permexact dots(0) gteffects
dis "C3_RC=" _rc "|gt_blank=" "`e(perm_gt_blank)'"
cap matrix drop PDSNAP
matrix PDSNAP = e(permdist_b)
wdcmp gtef , vals(`ov_plain')

* ---------- CASE 4: timetrends ----------
use "p06_data.dta", clear
cap noi wooldid y unit t ttre , att cluster(unit) permute(shuffle) permexact dots(0) timetrends
dis "C4_RC=" _rc
cap matrix drop PDSNAP
matrix PDSNAP = e(permdist_b)
oracle tt , asglist(`asglist') opts(timetrends)
wdcmp tt , vals(`oraclevals')
dis "P06_DONE"
