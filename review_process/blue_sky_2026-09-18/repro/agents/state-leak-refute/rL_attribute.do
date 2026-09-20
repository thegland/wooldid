* rL_attribute.do -- SKEPTIC: attribute the row-order sensitivity isolated in rK.
* rK: att+esprelength(3)+espostlength(3)+[aw=wt] moves 3.73e-08 relative between two row
* orders of the SAME rows; the moving cells are the pre3 horizon's se/t/p/lb/ub.
* Hypothesis H_regtol: it is reghdfe's iterative-absorption convergence tolerance
*   (wooldid default regtol(9) = tol(1e-9), sthlp:1117), not a wooldid store.
* Test: tighten/loosen regtol() and bypass reghdfe entirely (unconditionalse -> regress).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

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
gen double one = 1
sort st yr
save "rL_a.dta", replace
gen double shuf = runiform()
sort shuf
drop shuf
save "rL_b.dta", replace

cap prog drop cellcmp
prog def cellcmp
  syntax , a(string) b(string) tag(string)
  local R = rowsof(`a')
  local C = colsof(`a')
  local rn : rownames `a'
  local cn : colnames `a'
  scalar MAXR = 0
  forvalues i = 1/`R' {
    forvalues j = 1/`C' {
      local va = `a'[`i',`j']
      local vb = `b'[`i',`j']
      if `va' < . & `vb' < . {
        local d = abs(`va'-`vb')
        local den = max(abs(`va'),1e-300)
        local r = `d'/`den'
        if `r' > MAXR scalar MAXR = `r'
      }
    }
  }
  dis "MAXREL|`tag'|" MAXR
end

cap prog drop runpair
prog def runpair
  syntax , tag(string) opts(string) [wexp(string)]
  cap {
    use "rL_a.dta", clear
    qui wooldid y st yr onset `wexp', `opts'
    matrix LA = e(wdidmainresults)
    use "rL_b.dta", clear
    qui wooldid y st yr onset `wexp', `opts'
    matrix LB = e(wdidmainresults)
  }
  if _rc != 0 {
    dis "MAXREL|`tag'|FAILED_RC_" _rc
    exit
  }
  cellcmp , a(LA) b(LB) tag(`tag')
end

local ES esprelength(3) espostlength(3)

runpair , tag(L1_regtol9_default)  opts(cluster(st) att `ES')                  wexp([aw=wt])
runpair , tag(L2_regtol12)         opts(cluster(st) att `ES' regtol(12))       wexp([aw=wt])
runpair , tag(L3_regtol14)         opts(cluster(st) att `ES' regtol(14))       wexp([aw=wt])
runpair , tag(L4_regtol6)          opts(cluster(st) att `ES' regtol(6))        wexp([aw=wt])
runpair , tag(L5_uncondse_noHDFE)  opts(cluster(st) att `ES' unconditionalse)  wexp([aw=wt])
runpair , tag(L6_customvce_noHDFE) opts(att `ES' customvce(cluster st))        wexp([aw=wt])
runpair , tag(L7_pweights)         opts(cluster(st) att `ES')                  wexp([pw=wt])
runpair , tag(L8_int_aweights)     opts(cluster(st) att `ES')                  wexp([aw=iwt])
runpair , tag(L9_constant_aweight) opts(cluster(st) att `ES')                  wexp([aw=one])
runpair , tag(L10_unweighted)      opts(cluster(st) att `ES')
runpair , tag(L11_robust_aw)       opts(robust att `ES')                       wexp([aw=wt])
runpair , tag(L12_fe_aw)           opts(cluster(st) att `ES' fe(st))           wexp([aw=wt])

dis "ALLDONE"
