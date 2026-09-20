*==========================================================================
* r3.do — SKEPTIC, second attempt at NFC-01, on the SHAPE the prober used but
*         with my own generator: 150 units in 3 treated cohorts + never-treated,
*         10 periods, i = COHORT, unit effects present, control on a percentage
*         scale (50 + 8*x1) exactly as the claim describes.
*
* r1.do could not reproduce the claim on a state x individual x year micro panel
* (max rel on e(V) = 6.8e-16 with controls(age) and controls(inc)).  Two things
* differ from the prober's cell: (a) i = cohort rather than a many-valued unit id,
* (b) the control's contribution to the LINEAR PREDICTION is ~3 rather than ~0.3.
*
* The LADDER below is the decisive experiment: the same data, the same fitted
* coefficient, the control shifted by a pure constant c so that only |xb| moves.
* If the default-vs-nofastcollapse gap is a rounding quantum it must scale
* smoothly and stay at machine precision relative to the estimate.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 515151
set sortseed 515151

*------------------------------------------------------------------ data (own)
set obs 150
gen long unit = _n
gen double _u = runiform()
sort _u unit
gen long rk = _n
gen long cohort = 4
gen double ttre = .
replace cohort = 1 in 1/33
replace ttre = 2004 in 1/33
replace cohort = 2 in 34/66
replace ttre = 2006 in 34/66
replace cohort = 3 in 67/98
replace ttre = 2008 in 67/98
gen double a_i = rnormal()
gen long region = mod(unit,5)+1
drop _u rk
expand 10
bysort unit: gen int t = 1999 + _n
sort t unit
by t: gen double _f = rnormal() if _n==1
by t: replace _f = _f[1]
gen double ft = 0.5*_f
drop _f
gen double x1 = rnormal()
gen double x2 = 0.5*a_i + 0.5*rnormal()
gen byte post = (ttre < . & t >= ttre)
gen double tau = 0
replace tau = 0.50 if cohort==1 & post
replace tau = 0.80 if cohort==2 & post
replace tau = 0.30 if cohort==3 & post
gen double y = a_i + ft + tau + 0.5*x1 - 0.25*x2 + 0.30*rnormal()
* the percentage-scale control the claim names: 50 + 8*x1
gen double xpct = 50 + 8*x1
compress
save "_r3.dta", replace

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
dis _n as text "{hline 100}"
dis as text "arm                                        fcolA/B | maxrel V(all) | maxrel V(diag) | maxrel main | maxrel b"
dis as text "{hline 100}"

foreach arm in NOISE PCT PCTNOJ X1 NOCTRL {
  local ES att att_i att_it esprelength(3) espostlength(3) jointtests
  local ctrl controls(xpct)
  local optA
  local optB nofastcollapse
  local lbl
  if "`arm'"=="NOISE" {
    local lbl "NOISE  same arm twice, controls(xpct)"
    local optB
  }
  if "`arm'"=="PCT" {
    local lbl "PCT    controls(xpct), jointtests  <- claim cell"
  }
  if "`arm'"=="PCTNOJ" {
    local lbl "PCTNOJ controls(xpct), NO jointtests"
    local ES att att_i att_it esprelength(3) espostlength(3)
  }
  if "`arm'"=="X1" {
    local lbl "X1     controls(x1), scale ~1"
    local ctrl controls(x1)
  }
  if "`arm'"=="NOCTRL" {
    local lbl "NOCTRL no controls"
    local ctrl
  }
  use "_r3.dta", clear
  cap noi qui wooldid y cohort t ttre, `ES' `ctrl' `optA'
  local rcA = _rc
  if `rcA'==0 {
    matrix VA = e(V)
    matrix BA = e(b)
    matrix MA = e(wdidmainresults)
    local fA `e(wdid_fcol)'
  }
  use "_r3.dta", clear
  cap noi qui wooldid y cohort t ttre, `ES' `ctrl' `optB'
  local rcB = _rc
  if `rcB'==0 {
    matrix VB = e(V)
    matrix BB = e(b)
    matrix MB = e(wdidmainresults)
    local fB `e(wdid_fcol)'
  }
  if `rcA'!=0 | `rcB'!=0 {
    dis as error "`lbl' : rcA=`rcA' rcB=`rcB'"
    continue
  }
  mata: cmpm("VA","VB")
  scalar sVall = mxall
  scalar sVdia = mxdia
  mata: cmpm("MA","MB")
  scalar sMall = mxall
  mata: cmpm("BA","BB")
  scalar sBall = mxall
  dis as result %-42s "`lbl'" "  `fA'/`fB'  | " %9.2e sVall "   | " %9.2e sVdia "    | " %9.2e sMall "  | " %9.2e sBall
}
dis as text "{hline 100}"

*------------------------------------------------------------------ THE LADDER
* control = c + 8*x1.  The fitted coefficient (0.0625) is the same at every c;
* only the level of the linear prediction moves.
dis _n as text "{hline 100}"
dis as text "LADDER: controls(c + 8*x1).  Only |xb| moves."
dis as text "  c        mean|xb|     maxrel V(all)   maxrel V(diag)   maxrel main   maxrel b"
dis as text "{hline 100}"
foreach c in 0 5 50 500 5000 {
  use "_r3.dta", clear
  gen double xc = `c' + 8*x1
  qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) jointtests controls(xc)
  matrix VA = e(V)
  matrix BA = e(b)
  matrix MA = e(wdidmainresults)
  use "_r3.dta", clear
  gen double xc = `c' + 8*x1
  qui wooldid y cohort t ttre, att att_i att_it esprelength(3) espostlength(3) jointtests controls(xc) nofastcollapse
  matrix VB = e(V)
  matrix BB = e(b)
  matrix MB = e(wdidmainresults)
  * a proxy for |xb|: the control's own contribution, 0.0625*c
  local xbb = 0.0625*`c'
  mata: cmpm("VA","VB")
  scalar sVall = mxall
  scalar sVdia = mxdia
  mata: cmpm("MA","MB")
  scalar sMall = mxall
  mata: cmpm("BA","BB")
  scalar sBall = mxall
  dis as result "  c=" %-8.0g `c' "  " %9.3g `xbb' "     " %9.2e sVall "       " %9.2e sVdia "      " %9.2e sMall "    " %9.2e sBall
}
dis as text "{hline 100}"

dis "R3_DONE"
