*==========================================================================
* r1.do  — SKEPTIC's independent test of NFC-01.
*   Own data (state x individual x year micro panel, i = state), own generator,
*   own seeds, own metric.  No repo test-suite file, no repo .dta, no reuse of
*   the prober's script.
*
* QUESTION 1: does default-vs-nofastcollapse move e(V) above rtol 1e-10 when a
*             continuous control varies inside the collapse cell?
* QUESTION 2: is the headline a DIAGONAL (printed) move or an off-diagonal artifact?
* QUESTION 3: is the collapse special, or is margins' variance already
*             storage-order dependent at the same scale (shuffle arms)?
* QUESTION 4: how does it scale with the control's magnitude (z ~ N(0,1) vs
*             age ~ 35 vs income ~ 40000)?
* QUESTION 5: structural control — one row per cell (aggregated panel) must sit
*             at the floor; if it does not, the difference is not the collapse.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 71113
set sortseed 71113

*------------------------------------------------------------------ data
* 30 states, 25 individuals each, 8 years.  Treatment is at the STATE level and
* staggered; the microdata gives 25 rows per (state, year, targetflag) cell, which
* is what the collapse collapses.  This is an ordinary administrative shape.
local NS 30
local NI 25
local T  8
set obs `NS'
gen int state = _n
gen double u_s = rnormal()
gen int ttre = .
replace ttre = 2014 if state <= 10
replace ttre = 2016 if state > 10 & state <= 20
expand `NI'
bysort state: gen int pid = _n
gen double u_i = rnormal()
expand `T'
bysort state pid: gen int year = 2010 + _n
gen double age = 25 + 20*runiform()            // within-cell varying, scale ~35
gen double z   = rnormal()                     // within-cell varying, scale ~1
gen double inc = 30000 + 20000*runiform()      // within-cell varying, scale ~4e4
bysort state year (pid): gen double agecell = age[1]   // constant inside the cell
gen byte D = (ttre < . & year >= ttre)
gen double y = 2 + u_s + u_i + 0.05*(year-2010) + 0.4*D ///
             + 0.01*age + 0.10*z + 0.000004*inc + 0.5*rnormal()
label var age "age in years"
compress
save "_r1.dta", replace

* row-permuted twin: same data, different storage order
set seed 909090
gen double _s = runiform()
sort _s
drop _s
save "_r1_shuf.dta", replace

* structural control: ONE row per (state, year) — cells cannot vary internally
use "_r1.dta", clear
gcollapse (mean) y age z inc agecell, by(state year ttre)
save "_r1_agg.dta", replace

*------------------------------------------------------------------ comparator
mata:
void cmpm(string scalar an, string scalar bn) {
  real matrix A, B, S, Dm
  A = st_matrix(an); B = st_matrix(bn)
  if (rows(A)!=rows(B) | cols(A)!=cols(B)) {
    st_numscalar("mxall", .); st_numscalar("mxdia", .); return
  }
  S  = abs(A); S = S :* (S:>=abs(B)) + abs(B) :* (S:<abs(B))
  Dm = abs(A-B) :/ (S :+ (S:==0))
  st_numscalar("mxall", max(Dm))
  if (rows(A)==cols(A)) st_numscalar("mxdia", max(diagonal(Dm)))
  else st_numscalar("mxdia", .)
}
end

*------------------------------------------------------------------ arms
local ES att att_i att_it esprelength(2) espostlength(2)

dis _n as text "{hline 96}"
dis as text "arm                                      fcolA/B |  maxrel V(all) | maxrel V(diag) | maxrel b"
dis as text "{hline 96}"

foreach arm in NOISE NOCTRL Z AGE INC CELLCONST ONEROW SHUFOFF SHUFON {
  local dtaA "_r1.dta"
  local dtaB "_r1.dta"
  local ctrl controls(age)
  local optA
  local optB nofastcollapse
  local lbl
  if "`arm'"=="NOISE" {
    local lbl "NOISE   same arm twice, controls(age)"
    local optB
  }
  if "`arm'"=="NOCTRL" {
    local lbl "NOCTRL  no controls"
    local ctrl
  }
  if "`arm'"=="Z" {
    local lbl "Z       controls(z)      scale ~1"
    local ctrl controls(z)
  }
  if "`arm'"=="AGE" {
    local lbl "AGE     controls(age)    scale ~35"
  }
  if "`arm'"=="INC" {
    local lbl "INC     controls(inc)    scale ~4e4"
    local ctrl controls(inc)
  }
  if "`arm'"=="CELLCONST" {
    local lbl "CELLC   controls(agecell) cell-constant"
    local ctrl controls(agecell)
  }
  if "`arm'"=="ONEROW" {
    local lbl "ONEROW  aggregated panel, controls(age)"
    local dtaA "_r1_agg.dta"
    local dtaB "_r1_agg.dta"
  }
  if "`arm'"=="SHUFOFF" {
    local lbl "SHUFOFF row shuffle, collapse OFF both"
    local optA nofastcollapse
    local dtaB "_r1_shuf.dta"
  }
  if "`arm'"=="SHUFON" {
    local lbl "SHUFON  row shuffle, collapse ON both"
    local optB
    local dtaB "_r1_shuf.dta"
  }

  use "`dtaA'", clear
  cap noi qui wooldid y state year ttre, `ES' `ctrl' `optA'
  local rcA = _rc
  if `rcA'==0 {
    matrix VA = e(V)
    matrix BA = e(b)
    matrix MA = e(wdidmainresults)
    local fA `e(wdid_fcol)'
    local ncA `e(wdid_fcncell)'
  }
  use "`dtaB'", clear
  cap noi qui wooldid y state year ttre, `ES' `ctrl' `optB'
  local rcB = _rc
  if `rcB'==0 {
    matrix VB = e(V)
    matrix BB = e(b)
    matrix MB = e(wdidmainresults)
    local fB `e(wdid_fcol)'
    local ncB `e(wdid_fcncell)'
  }
  if `rcA'!=0 | `rcB'!=0 {
    dis as error "`lbl' : rcA=`rcA' rcB=`rcB'"
    continue
  }
  mata: cmpm("VA","VB")
  scalar sVall = mxall
  scalar sVdia = mxdia
  mata: cmpm("BA","BB")
  scalar sBall = mxall
  mata: cmpm("MA","MB")
  scalar sMall = mxall
  dis as result %-40s "`lbl'" "   `fA'/`fB'   |  " %9.2e sVall "     |  " %9.2e sVdia "     | " %9.2e sBall
  dis as text  "        (cells `ncA'/`ncB';  max rel on e(wdidmainresults) = " %9.2e sMall ")"
}
dis as text "{hline 96}"

*------------------------------------------------------------------ the printed digits
* hex + wide-format dump of the FIRST diagonal SE, both arms, controls(age) and controls(inc)
foreach c in age inc {
  foreach a in ON OFF {
    local o
    if "`a'"=="OFF" local o nofastcollapse
    use "_r1.dta", clear
    qui wooldid y state year ttre, att att_i att_it esprelength(2) espostlength(2) controls(`c') `o'
    matrix VV = e(V)
    scalar s1 = sqrt(VV[1,1])
    dis as text "SE[1] controls(`c') collapse-`a' : " %21x s1 "   " %20.18g s1
  }
}

*------------------------------------------------------------------ printed table diff
dis "<<<<< PRINTED_A collapse ON  controls(inc age) >>>>>"
use "_r1.dta", clear
cap noi wooldid y state year ttre, att att_it esprelength(2) espostlength(2) controls(inc age)
dis "<<<<< PRINTED_B nofastcollapse controls(inc age) >>>>>"
use "_r1.dta", clear
cap noi wooldid y state year ttre, att att_it esprelength(2) espostlength(2) controls(inc age) nofastcollapse
dis "<<<<< PRINTED_END >>>>>"

dis "R1_DONE"
