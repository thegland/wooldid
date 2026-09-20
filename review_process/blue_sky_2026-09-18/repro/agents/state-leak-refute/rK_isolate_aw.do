* rK_isolate_aw.do -- SKEPTIC: isolate the 3.7e-08 row-order sensitivity found in rJ arm J2.
* J1 (unweighted, att/att_i/att_it/att_itime) moved only 3.6e-15 under a row shuffle.
* J2 ([aw=wt] + esprelength/espostlength + jointtests) moved 3.7e-08 -- 370x the project's
* 1e-10 indifference threshold, and about 0.3 float-ulp, the signature of a single-precision store.
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
gen double wt = 1 + runiform()
sort st yr
save "rK_a.dta", replace
gen double shuf = runiform()
sort shuf
drop shuf
save "rK_b.dta", replace

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
        if `r' > 1e-10 {
          local rr : word `i' of `rn'
          local cc : word `j' of `cn'
          dis "CELL|`tag'|`rr'|`cc'|a=" %21x `va' "|b=" %21x `vb' "|rel=" `r'
        }
        if `r' > MAXR scalar MAXR = `r'
      }
    }
  }
  dis "MAXREL|`tag'|" MAXR
end

cap prog drop runpair
prog def runpair
  syntax , tag(string) opts(string) [wexp(string)]
  use "rK_a.dta", clear
  qui wooldid y st yr onset `wexp', `opts'
  matrix KA = e(wdidmainresults)
  use "rK_b.dta", clear
  qui wooldid y st yr onset `wexp', `opts'
  matrix KB = e(wdidmainresults)
  cellcmp , a(KA) b(KB) tag(`tag')
end

runpair , tag(K1_att_only_noW)            opts(cluster(st) att)
runpair , tag(K2_att_only_aw)             opts(cluster(st) att) wexp([aw=wt])
runpair , tag(K3_att_es_noW)              opts(cluster(st) att esprelength(3) espostlength(3))
runpair , tag(K4_att_es_aw)               opts(cluster(st) att esprelength(3) espostlength(3)) wexp([aw=wt])
runpair , tag(K5_att_es_jt_noW)           opts(cluster(st) att esprelength(3) espostlength(3) jointtests)
runpair , tag(K6_att_es_jt_aw)            opts(cluster(st) att esprelength(3) espostlength(3) jointtests) wexp([aw=wt])
runpair , tag(K7_att_only_pw)             opts(cluster(st) att) wexp([pw=wt])
runpair , tag(K8_att_only_fw)             opts(cluster(st) att) wexp([fw=wt])
runpair , tag(K9_atti_aw)                 opts(cluster(st) att_i att_it att_itime) wexp([aw=wt])
runpair , tag(K10_att_aw_robust)          opts(robust att) wexp([aw=wt])

* K11: is the weight variable's STORAGE the driver? recast wt to float and repeat K2.
use "rK_a.dta", clear
recast float wt
save "rK_af.dta", replace
use "rK_b.dta", clear
recast float wt
save "rK_bf.dta", replace
use "rK_af.dta", clear
qui wooldid y st yr onset [aw=wt], cluster(st) att
matrix KA = e(wdidmainresults)
use "rK_bf.dta", clear
qui wooldid y st yr onset [aw=wt], cluster(st) att
matrix KB = e(wdidmainresults)
cellcmp , a(KA) b(KB) tag(K11_att_aw_floatwt)

* K12: integer weights (exactly representable) -- does the sensitivity vanish?
use "rK_a.dta", clear
gen double iwt = 1 + mod(_n,5)
save "rK_ai.dta", replace
use "rK_b.dta", clear
gen double iwt = 1 + mod(st,5)
qui replace iwt = 1 + mod(st,5)
save "rK_bi.dta", replace
use "rK_ai.dta", clear
qui replace iwt = 1 + mod(st,5)
save "rK_ai.dta", replace
use "rK_ai.dta", clear
qui wooldid y st yr onset [aw=iwt], cluster(st) att
matrix KA = e(wdidmainresults)
use "rK_bi.dta", clear
qui wooldid y st yr onset [aw=iwt], cluster(st) att
matrix KB = e(wdidmainresults)
cellcmp , a(KA) b(KB) tag(K12_att_aw_intweights)

dis "ALLDONE"
