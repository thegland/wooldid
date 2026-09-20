* ===================================================================
* probe2_surface.do -- SELF-CONTAINED
* Consequence surface of the fastnose event-study posting:
*  E1  A vs B bit comparison of e(b_es_att) (analytic vs fastnose)
*  E2  fastnose + subgroups(): does any SE leak into the sg tables?
*  E3  fastnose + contreat(): ame family behaviour
*  E4  permute() fast mode: e(b_es_att) must NOT post (carve-out)
*  E5  fastnose + timetrends + honestdid-signature checks
* ===================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

* ---------------- panel (same DGP as probe1, plus sg + dose) ----------
clear
set obs 48
gen long i = _n
gen byte coh = mod(_n,4)
gen int ttre = .
replace ttre = 6  if coh==1
replace ttre = 9  if coh==2
replace ttre = 12 if coh==3
gen byte sg = mod(_n,2)
gen double ui = rnormal()
gen double dose_u = 0.5 + runiform()
expand 16
bysort i: gen int t = _n
gen double et = rnormal()*0.5
gen byte D = (ttre<. & t>=ttre)
gen double rel = t - ttre
gen double dose = cond(D==1, dose_u, 0)
gen double y = 1 + 0.30*ui + 0.05*t + cond(D==1, 0.8 + 0.15*rel, 0) + 0.4*sg + et
gen double yc = 1 + 0.30*ui + 0.05*t + 0.9*dose + et
assert !missing(y) & !missing(yc)
compress
save "panel2.dta", replace

* =================== E1: A vs B bit comparison =====================
dis "=== E1a analytic ==="
use "panel2.dta", clear
qui wooldid y i t ttre, att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod
matrix E1A = e(b_es_att)
dis "E1A_NUMPRE=" e(numpre_att)
dis "=== E1b fastnose ==="
use "panel2.dta", clear
qui wooldid y i t ttre, att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod fastnose
matrix E1B = e(b_es_att)
dis "E1B_NUMPRE=" e(numpre_att)
local nA : colnames E1A
local nB : colnames E1B
dis "E1_NAMES_MATCH=" ("`nA'"=="`nB'")
dis "E1_NAMES=" "`nA'"
local k = colsof(E1A)
local maxad 0
forval j = 1/`k' {
  local d = abs(E1A[1,`j']-E1B[1,`j'])
  if `d' > `maxad' local maxad = `d'
  dis "E1_CELL " `j' " A=" %21x E1A[1,`j'] " B=" %21x E1B[1,`j'] " same=" (E1A[1,`j']==E1B[1,`j'])
}
dis "E1_MAXABSDIFF=" `maxad'

* =================== E2: fastnose + subgroups ======================
dis "=== E2 fastnose + subgroups ==="
use "panel2.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod subgroups(sg) fastnose
dis "E2_RC=" _rc
cap confirm matrix e(V)
dis "E2_V_RC=" _rc
foreach s in 0 1 {
  cap confirm matrix e(b_es_att_s`s')
  dis "E2_BES_s`s'_RC=" _rc
  cap confirm matrix e(V_es_att_s`s')
  dis "E2_VES_s`s'_RC=" _rc
  cap confirm matrix e(wdidsgresults_`s')
  dis "E2_SGTAB_s`s'_RC=" _rc
  cap mat list e(wdidsgresults_`s')
}
cap confirm matrix e(pwsgs)
dis "E2_PWSGS_RC=" _rc
cap mat list e(pwsgs)
cap mat list e(wdidmainresults)
* explicit leak test: any non-missing se/t/p/lb/ub anywhere in the sg tables?
foreach mm in wdidmainresults wdidsgresults_0 wdidsgresults_1 pwsgs {
  cap confirm matrix e(`mm')
  if _rc == 0 {
    tempname M
    matrix `M' = e(`mm')
    local nr = rowsof(`M')
    local nc = colsof(`M')
    local cn : colnames `M'
    local leak 0
    forval rr = 1/`nr' {
      local cc 0
      foreach c of local cn {
        local ++cc
        if inlist("`c'","se","t","z","p","lb_95ci","ub_95ci") {
          if !missing(`M'[`rr',`cc']) local ++leak
        }
      }
    }
    dis "E2_LEAK_`mm'=" `leak'
  }
}

* =================== E3: fastnose + contreat =======================
dis "=== E3 fastnose + contreat ==="
use "panel2.dta", clear
cap noi wooldid yc i t ttre, att cluster(i) contreat(dose) esprelength(3) espostlength(3) esfixedbaseperiod fastnose
dis "E3_RC=" _rc
foreach f in ame att {
  cap confirm matrix e(b_es_`f')
  dis "E3_BES_`f'_RC=" _rc
  cap confirm matrix e(V_es_`f')
  dis "E3_VES_`f'_RC=" _rc
  dis "E3_NUMPRE_`f'=" e(numpre_`f')
}
dis "=== E3b contreat analytic twin ==="
use "panel2.dta", clear
cap noi wooldid yc i t ttre, att cluster(i) contreat(dose) esprelength(3) espostlength(3) esfixedbaseperiod
dis "E3b_RC=" _rc
foreach f in ame att {
  cap confirm matrix e(b_es_`f')
  dis "E3b_BES_`f'_RC=" _rc
  cap confirm matrix e(V_es_`f')
  dis "E3b_VES_`f'_RC=" _rc
  dis "E3b_NUMPRE_`f'=" e(numpre_`f')
}

* =================== E4: permute fast mode =========================
dis "=== E4 permute fast mode ==="
use "panel2.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(3) espostlength(3) esfixedbaseperiod permute(ttre) permreps(20) rseed(7)
dis "E4_RC=" _rc
cap confirm matrix e(b_es_att)
dis "E4_BES_RC=" _rc
cap confirm matrix e(V_es_att)
dis "E4_VES_RC=" _rc
dis "E4_NUMPRE=" e(numpre_att)

dis "=== PROBE2 DONE ==="
