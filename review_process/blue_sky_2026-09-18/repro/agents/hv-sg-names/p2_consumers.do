*==========================================================================
* p2_consumers.do -- do the DOWNSTREAM tables the help names (permtstat,
*   bootsave(), bootdist, perm kernels) address a coefficient by the e(b)
*   name, and what is the ame-family naming?  Plus: negative / non-integer
*   subgroup values.  SELF-CONTAINED.
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
  di "NCOL|`tag'|`k'"
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
  local rn : rownames `M'
  di "MROWS|`tag'|`mname'|`rn'"
end

* ---------------- DGP ----------------
set obs 80
gen long unit = _n
gen double u_fe = rnormal()
gen byte sgA = 1 + mod(_n, 2)
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
save "p2_data.dta", replace

* ============ G: continuous treatment (ame family) + subgroups ============
di "=== RUN G_amesg ==="
cap noi wooldid y unit year ttre, att att_it contreat(dose) contreatcontrols(basic) ///
    subgroups(sgA) esprelength(2) espostlength(2) robust
di "G_RC=" _rc
if _rc == 0 dumpn "G_amesg"

* ============ H: bootstrap + boottstat + bootsave ============
di "=== RUN H_boot ==="
cap erase p2boot.dta
cap noi wooldid y unit year ttre, att att_i subgroups(sgA) esprelength(2) espostlength(2) ///
    robust bootstrap(standard) reps(25) rseed(777) boottstat bootsave(p2boot) replace
di "H_RC=" _rc
if _rc == 0 {
  dumpn "H_boot"
  dumpmat bootdist       "H_boot"
  dumpmat boot_tstud     "H_boot"
  dumpmat boot_ptails    "H_boot"
  dumpmat boot_wald_mainresults "H_boot"
}

* ============ I: permute + permtstat + permsave ============
di "=== RUN I_perm ==="
cap erase p2perm.dta
use "p2_data.dta", clear
cap noi wooldid y unit year ttre, att att_i subgroups(sgA) esprelength(2) espostlength(2) ///
    robust permute(shuffle) permreps(25) rseed(777) permtstat permsave(p2perm) replace
di "I_RC=" _rc
if _rc == 0 {
  dumpn "I_perm"
  dumpmat perm_stat  "I_perm"
  dumpmat perm_n     "I_perm"
  dumpmat perm_p_lo  "I_perm"
}

* ============ inspect the two saved replication datasets ============
foreach f in p2boot p2perm {
  cap confirm file "`f'.dta"
  if _rc == 0 {
    use "`f'.dta", clear
    qui ds
    local vl `r(varlist)'
    foreach v of local vl {
      di "SAVEVAR|`f'|`v'"
    }
  }
  else di "SAVEMISS|`f'"
}

* ============ J: negative subgroup value ============
di "=== RUN J_negsg ==="
use "p2_data.dta", clear
gen int sgN = cond(mod(unit,2)==0, -1, 1)
cap noi wooldid y unit year ttre, att subgroups(sgN) robust
di "J_RC=" _rc
if _rc == 0 dumpn "J_negsg"

* ============ K: non-integer subgroup value ============
di "=== RUN K_fracsg ==="
use "p2_data.dta", clear
gen double sgF = cond(mod(unit,2)==0, 1.5, 2)
cap noi wooldid y unit year ttre, att subgroups(sgF) robust
di "K_RC=" _rc
if _rc == 0 dumpn "K_fracsg"

* ============ L: zero-valued subgroup (legal: non-negative integer) ======
di "=== RUN L_zerosg ==="
use "p2_data.dta", clear
gen byte sgZ = mod(unit,2)
cap noi wooldid y unit year ttre, att subgroups(sgZ) esprelength(2) espostlength(2) robust
di "L_RC=" _rc
if _rc == 0 dumpn "L_zerosg"

di "PROBE2_DONE"
