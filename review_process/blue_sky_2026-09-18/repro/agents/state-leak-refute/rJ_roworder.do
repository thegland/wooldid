* rJ_roworder.do -- SKEPTIC new-probe: is wooldid's answer invariant to the caller's ROW ORDER?
* Same rows, different physical order. A user who merges or sorts differently must get the
* same number. Project indifference threshold: rtol 1e-10.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

*--------------------------------------------------- unit-level panel
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
save "rJ_unit.dta", replace
gen double shuf = runiform()
sort shuf
drop shuf
save "rJ_unit_shuf.dta", replace
use "rJ_unit.dta", clear
gsort -st -yr
save "rJ_unit_rev.dta", replace

*--------------------------------------------------- microdata (collapse path)
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
expand 40
gen double y = sfe + 0.02*(yr-2004) + 0.35*((yr>=onset) & onset<.) + rnormal()
gen double wt = 1 + runiform()
sort st yr
save "rJ_micro.dta", replace
gen double shuf = runiform()
sort shuf
drop shuf
save "rJ_micro_shuf.dta", replace

cap prog drop cmpmat
prog def cmpmat
  syntax , a(string) b(string) tag(string)
  mata: st_numscalar("MABS", max(abs(st_matrix("`a'"):-st_matrix("`b'"))))
  mata: st_numscalar("MREL", max(abs((st_matrix("`a'"):-st_matrix("`b'")):/editmissing(abs(st_matrix("`a'")),1))))
  dis "CMP|`tag'|maxabs|" %21x MABS
  dis "CMP|`tag'|maxrel|" %21x MREL
  dis "CMP|`tag'|maxrel_dec|" MREL
  dis "CMP|`tag'|BIT|" (MABS==0)
end

*=========== J1 unit-level, plain att
use "rJ_unit.dta", clear
qui wooldid y st yr onset, cluster(st) att att_i att_it att_itime
matrix J1a = e(wdidmainresults)
use "rJ_unit_shuf.dta", clear
qui wooldid y st yr onset, cluster(st) att att_i att_it att_itime
matrix J1b = e(wdidmainresults)
use "rJ_unit_rev.dta", clear
qui wooldid y st yr onset, cluster(st) att att_i att_it att_itime
matrix J1c = e(wdidmainresults)
cmpmat , a(J1a) b(J1b) tag(J1_unit_shuffled)
cmpmat , a(J1a) b(J1c) tag(J1_unit_reversed)

*=========== J2 unit-level, event study + weights
use "rJ_unit.dta", clear
qui wooldid y st yr onset [aw=wt], cluster(st) att esprelength(3) espostlength(3) jointtests
matrix J2a = e(wdidmainresults)
use "rJ_unit_shuf.dta", clear
qui wooldid y st yr onset [aw=wt], cluster(st) att esprelength(3) espostlength(3) jointtests
matrix J2b = e(wdidmainresults)
cmpmat , a(J2a) b(J2b) tag(J2_unit_es_weighted)

*=========== J3 microdata, default (fast collapse)
use "rJ_micro.dta", clear
qui wooldid y st yr onset, cluster(st) att att_i att_it att_itime
matrix J3a = e(wdidmainresults)
use "rJ_micro_shuf.dta", clear
qui wooldid y st yr onset, cluster(st) att att_i att_it att_itime
matrix J3b = e(wdidmainresults)
cmpmat , a(J3a) b(J3b) tag(J3_micro_fastcollapse)

*=========== J4 microdata, nofastcollapse
use "rJ_micro.dta", clear
qui wooldid y st yr onset, cluster(st) att att_i att_it att_itime nofastcollapse
matrix J4a = e(wdidmainresults)
use "rJ_micro_shuf.dta", clear
qui wooldid y st yr onset, cluster(st) att att_i att_it att_itime nofastcollapse
matrix J4b = e(wdidmainresults)
cmpmat , a(J4a) b(J4b) tag(J4_micro_nofast)

*=========== J5 bootstrap with rseed() -- same rows, different order
use "rJ_unit.dta", clear
set seed 55
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(50) rseed(4242)
matrix J5a = e(wdidmainresults)
use "rJ_unit_shuf.dta", clear
set seed 55
qui wooldid y st yr onset, cluster(st) att bootstrap(standard) reps(50) rseed(4242)
matrix J5b = e(wdidmainresults)
cmpmat , a(J5a) b(J5b) tag(J5_boot_rseed)

*=========== J6 permute(shuffle) -- same rows, different order, same seed
use "rJ_unit.dta", clear
set seed 66
qui wooldid y st yr onset, cluster(st) att permute(shuffle) permreps(50)
matrix J6a = e(wdidmainresults)
use "rJ_unit_shuf.dta", clear
set seed 66
qui wooldid y st yr onset, cluster(st) att permute(shuffle) permreps(50)
matrix J6b = e(wdidmainresults)
cmpmat , a(J6a) b(J6b) tag(J6_permute_shuffle)

*=========== J7 same order, different sortseed -- does set sortseed move anything?
use "rJ_unit.dta", clear
set sortseed 111
qui wooldid y st yr onset, cluster(st) att att_i att_it att_itime
matrix J7a = e(wdidmainresults)
use "rJ_unit.dta", clear
set sortseed 999999
qui wooldid y st yr onset, cluster(st) att att_i att_it att_itime
matrix J7b = e(wdidmainresults)
cmpmat , a(J7a) b(J7b) tag(J7_sortseed)

dis "ALLDONE"
