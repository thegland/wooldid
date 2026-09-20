*==========================================================================
* p3_perm_and_mats.do -- (1) do the permute()/permtstat tables address a
*   coefficient by the e(b) name?  (2) full e() matrix inventory so the
*   THREE subgroup naming schemes can be tabulated.  SELF-CONTAINED.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

cap prog drop dumpn
prog def dumpn
  args tag
  local nn : colnames e(b)
  local k = 0
  foreach n of local nn {
    local k = `k' + 1
    di "NAME|`tag'|`k'|`n'"
  }
end
cap prog drop dumpmat
prog def dumpmat
  args mname tag
  cap confirm matrix e(`mname')
  if _rc != 0 {
    di "MATMISS|`tag'|`mname'"
    exit
  }
  tempname M
  matrix `M' = e(`mname')
  local cn : colnames `M'
  local k = 0
  foreach c of local cn {
    local k = `k' + 1
    di "MCOL|`tag'|`mname'|`k'|`c'"
  }
end
cap prog drop dumpinv
prog def dumpinv
  args tag
  local M : e(matrices)
  foreach m of local M {
    di "INV|`tag'|`m'"
  }
end

* ---------------- DGP ----------------
set obs 80
gen long unit = _n
gen double u_fe = rnormal()
gen byte sgA = 1 + mod(_n, 2)
gen byte sgB = cond(mod(_n,2)==0, 3, 17)
gen byte cohcode = mod(_n, 4)
gen int ttre = .
replace ttre = 2003 if cohcode == 1
replace ttre = 2005 if cohcode == 2
replace ttre = 2007 if cohcode == 3
expand 10
bysort unit: gen int year = 1999 + _n
gen double t_fe = 0.2*(year - 2000)
gen double post = (year >= ttre) & ttre < .
gen double tau  = 1.0 + 0.3*(year - ttre) if post
replace tau = 0 if !post
gen double y = u_fe + t_fe + tau + 0.30*rnormal()
gen double dose = 1 + abs(u_fe) + 0.2*(year-2000)
save "p3_data.dta", replace

* ============ N: subgroup values 3 / 17 -- full e() inventory ============
di "=== RUN N_sgB ==="
cap noi wooldid y unit year ttre, att att_i subgroups(sgB) esprelength(2) espostlength(2) robust
di "N_RC=" _rc
if _rc == 0 {
  dumpinv "N_sgB"
  foreach cand in b_es_att_s3 b_es_att_s17 b_es_att_s1 b_es_att_s2 ///
                  wdidsgresults_3 wdidsgresults_17 wdidsgresults_1 {
    cap confirm matrix e(`cand')
    di "CONFIRM|N|`cand'|rc=" _rc
  }
  * scalar naming for subgroups
  foreach sc in adversarial_pval_att_s3 adversarial_pval_att_SG3 joint_att_allsgssame {
    cap confirm scalar e(`sc')
    di "SCALAR|N|`sc'|rc=" _rc
  }
}

* ============ O: contreat + subgroups -- inventory (ame subgroup ES?) ====
di "=== RUN O_amesg ==="
use "p3_data.dta", clear
cap noi wooldid y unit year ttre, att att_it contreat(dose) contreatcontrols(basic) ///
    subgroups(sgA) esprelength(2) espostlength(2) robust
di "O_RC=" _rc
if _rc == 0 dumpinv "O_amesg"

* ============ M: permute + permtstat + permsave ============
di "=== RUN M_perm ==="
use "p3_data.dta", clear
cap erase p3perm.dta
cap noi wooldid y unit year ttre, att subgroups(sgA) esprelength(1) espostlength(1) ///
    cluster(unit) permute(shuffle) permreps(20) rseed(777) permtstat permsave(p3perm) replace
di "M_RC=" _rc
if _rc == 0 {
  dumpn "M_perm"
  dumpinv "M_perm"
  foreach mm in perm_stat perm_n perm_p_lo perm_comparator_b {
    dumpmat `mm' "M_perm"
  }
}
cap confirm file "p3perm.dta"
if _rc == 0 {
  use "p3perm.dta", clear
  qui ds
  foreach v of local r(varlist) {
    di "SAVEVAR|p3perm|`v'"
  }
}
else di "SAVEMISS|p3perm"

di "PROBE3_DONE"
