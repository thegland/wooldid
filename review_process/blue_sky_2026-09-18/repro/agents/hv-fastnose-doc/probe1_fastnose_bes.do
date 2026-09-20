* ===================================================================
* probe1_fastnose_bes.do  -- SELF-CONTAINED
* Question: on the pinned build, does `fastnose` post e(b_es_att) and
* e(numpre_att)?  And what do the analytic / timetrends / impvalues
* twins post?  The sthlp sentence at :1325 says the NO-identified-
* variance families (fastnose, typically timetrends, all-zero VCE)
* post b without V.
* ===================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

* ---------- a plain, well-behaved staggered panel ----------
* 48 units, 16 periods, three treated cohorts + never-treated
clear
set obs 48
gen long i = _n
gen byte coh = mod(_n,4)            // 0 = never treated
gen int ttre = .
replace ttre = 6  if coh==1
replace ttre = 9  if coh==2
replace ttre = 12 if coh==3
gen double ui = rnormal()
expand 16
bysort i: gen int t = _n
gen double et = rnormal()*0.5
gen byte D = (ttre<. & t>=ttre)
gen double rel = t - ttre
gen double y = 1 + 0.30*ui + 0.05*t + cond(D==1, 0.8 + 0.15*rel, 0) + et
assert !missing(y)
compress
save "panel1.dta", replace

* ---------- helper macro text repeated inline for self-containment ----
* ---------------- A: analytic twin ----------------
dis "=== RUN A: analytic ==="
use "panel1.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod
dis "A_RC=" _rc
cap confirm matrix e(b_es_att)
dis "A_BES_RC=" _rc
cap confirm matrix e(V_es_att)
dis "A_VES_RC=" _rc
dis "A_NUMPRE=" e(numpre_att)
dis "A_VCVSTATE=" "`e(wdid_vcvstate)'"
cap mat list e(b_es_att), format(%21x)
cap local anA : colnames e(b_es_att)
dis "A_BESNAMES=" "`anA'"
cap local bn : colnames e(b)
dis "A_BNAMES=" "`bn'"
mat A_bes = e(b_es_att)

* ---------------- B: fastnose ----------------
dis "=== RUN B: fastnose ==="
use "panel1.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod fastnose
dis "B_RC=" _rc
cap confirm matrix e(b_es_att)
dis "B_BES_RC=" _rc
cap confirm matrix e(V_es_att)
dis "B_VES_RC=" _rc
dis "B_NUMPRE=" e(numpre_att)
dis "B_VCVSTATE=" "`e(wdid_vcvstate)'"
cap confirm matrix e(V)
dis "B_V_RC=" _rc
cap mat list e(b_es_att), format(%21x)
cap local anB : colnames e(b_es_att)
dis "B_BESNAMES=" "`anB'"
cap local bnB : colnames e(b)
dis "B_BNAMES=" "`bnB'"
cap mat B_bes = e(b_es_att)

* bit-comparison A vs B
cap {
  local kA = colsof(A_bes)
  local kB = colsof(B_bes)
  dis "AB_COLS=" `kA' " " `kB'
  if `kA'==`kB' {
    local maxad 0
    forval j=1/`kA' {
      local d = abs(A_bes[1,`j']-B_bes[1,`j'])
      if `d' > `maxad' local maxad = `d'
      dis "AB_CELL " `j' " A=" %21x A_bes[1,`j'] " B=" %21x B_bes[1,`j']
    }
    dis "AB_MAXABSDIFF=" %21x `maxad'
  }
}
dis "AB_CMP_RC=" _rc

* ---------------- C: timetrends (analytic) ----------------
dis "=== RUN C: timetrends ==="
use "panel1.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod timetrends
dis "C_RC=" _rc
cap confirm matrix e(b_es_att)
dis "C_BES_RC=" _rc
cap confirm matrix e(V_es_att)
dis "C_VES_RC=" _rc
dis "C_NUMPRE=" e(numpre_att)
dis "C_VCVSTATE=" "`e(wdid_vcvstate)'"

* ---------------- D: impvalues ----------------
dis "=== RUN D: impvalues ==="
use "panel1.dta", clear
cap noi wooldid y i t ttre, att cluster(i) esprelength(5) espostlength(6) esfixedbaseperiod impvalues
dis "D_RC=" _rc
cap confirm matrix e(b_es_att)
dis "D_BES_RC=" _rc
cap confirm matrix e(V_es_att)
dis "D_VES_RC=" _rc
dis "D_NUMPRE=" e(numpre_att)
dis "D_VCVSTATE=" "`e(wdid_vcvstate)'"

dis "=== PROBE1 DONE ==="
