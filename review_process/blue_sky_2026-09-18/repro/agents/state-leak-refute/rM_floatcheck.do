* rM_floatcheck.do -- SKEPTIC: is the row-order sensitivity a surviving B21-class FLOAT store
* inside wooldid, or is it reghdfe's?
* Suspects: ado:1619-1621, three UNTYPED gegen totals that build the cell weight denominators
*   qui gegen `w_per_it' = total(`esamp' * `w')  , by(`i' `t')
*   qui gegen `w_per_i'  = total(`w'*`esamp')    , by(`posttreatment' `i')
*   qui gegen `t_per_i'  = total(`tag')          , by(`posttreatment' `i')
* If untyped gegen defaults to float here, w_att_it/w_att_i/w_att_itime divide by a float total.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid
dis "DEFAULT_TYPE=" "`c(type)'"

*=========== M0: what storage type does an UNTYPED gegen total() get here?
clear
set obs 100
gen int g = 1 + mod(_n,7)
gen double x = runiform()
qui gegen tot_untyped = total(x), by(g)
qui gegen double tot_double = total(x), by(g)
qui egen  tot_egen = total(x), by(g)
local ty : type tot_untyped
local td : type tot_double
local te : type tot_egen
dis "M0|gegen_untyped_type|`ty'"
dis "M0|gegen_double_type|`td'"
dis "M0|egen_untyped_type|`te'"
qui gen double dd = abs(tot_untyped - tot_double)
qui su dd, meanonly
dis "M0|maxabs_untyped_vs_double|" %21x r(max)

*=========== data
clear
set obs 30
gen int st = _n
gen int onset = .
replace onset = 2008 if st <= 8
replace onset = 2011 if st > 8 & st <= 16
replace onset = 2013 if st > 16 & st <= 22
gen double sfe = rnormal()
expand 10
bys st: gen int yr = 2004 + _n
gen double y = sfe + 0.02*(yr-2004) + 0.35*((yr>=onset) & onset<.) + rnormal()*0.5
gen double wt  = 1 + runiform()
gen double iwt = 1 + mod(st,5)
gen byte post = (yr >= onset) & onset < .
sort st yr
save "rM_a.dta", replace
gen double shuf = runiform()
sort shuf
drop shuf
save "rM_b.dta", replace

cap prog drop cellcmp
prog def cellcmp
  syntax , a(string) b(string) tag(string)
  local R = rowsof(`a')
  local C = colsof(`a')
  scalar MAXR = 0
  forvalues i = 1/`R' {
    forvalues j = 1/`C' {
      local va = `a'[`i',`j']
      local vb = `b'[`i',`j']
      if `va' < . & `vb' < . {
        local r = abs(`va'-`vb')/max(abs(`va'),1e-300)
        if `r' > MAXR scalar MAXR = `r'
      }
    }
  }
  dis "MAXREL|`tag'|" MAXR
end

cap prog drop runpair
prog def runpair
  syntax , tag(string) opts(string) [wexp(string) settype(string)]
  cap {
    if "`settype'" != "" set type `settype'
    use "rM_a.dta", clear
    qui wooldid y st yr onset `wexp', `opts'
    matrix MA = e(wdidmainresults)
    use "rM_b.dta", clear
    qui wooldid y st yr onset `wexp', `opts'
    matrix MB = e(wdidmainresults)
    if "`settype'" != "" set type float
  }
  if _rc != 0 {
    dis "MAXREL|`tag'|FAILED_RC_" _rc
    if "`settype'" != "" cap set type float
    exit
  }
  cellcmp , a(MA) b(MB) tag(`tag')
end

local ES esprelength(3) espostlength(3)

* --- the firing arm, and the same arm under `set type double'
runpair , tag(M1_att_es_aw_typefloat) opts(cluster(st) att `ES') wexp([aw=wt])
runpair , tag(M2_att_es_aw_typedouble) opts(cluster(st) att `ES') wexp([aw=wt]) settype(double)

* --- the cell-weighted estimands, which DO build the untyped gegen totals (ado:1619-1621)
runpair , tag(M3_cellwtd_aw_typefloat)  opts(cluster(st) att_it att_i att_itime `ES') wexp([aw=wt])
runpair , tag(M4_cellwtd_aw_typedouble) opts(cluster(st) att_it att_i att_itime `ES') wexp([aw=wt]) settype(double)
runpair , tag(M5_cellwtd_unweighted)    opts(cluster(st) att_it att_i att_itime `ES')
runpair , tag(M6_cellwtd_intw)          opts(cluster(st) att_it att_i att_itime `ES') wexp([aw=iwt])

*=========== M7: reghdfe ALONE on the two orderings -- attribution
cap prog drop rhpair
prog def rhpair
  syntax , tag(string) [wexp(string)]
  use "rM_a.dta", clear
  qui reghdfe y post i.st#c.post `wexp', absorb(st yr) cluster(st) tol(1e-9)
  matrix RA = e(V)
  matrix RBb = e(b)
  use "rM_b.dta", clear
  qui reghdfe y post i.st#c.post `wexp', absorb(st yr) cluster(st) tol(1e-9)
  matrix RC = e(V)
  matrix RDb = e(b)
  cellcmp , a(RA) b(RC) tag(`tag'_V)
  cellcmp , a(RBb) b(RDb) tag(`tag'_b)
end
rhpair , tag(M7_reghdfe_alone_aw) wexp([aw=wt])
rhpair , tag(M8_reghdfe_alone_intw) wexp([aw=iwt])
rhpair , tag(M9_reghdfe_alone_now)

dis "ALLDONE"
