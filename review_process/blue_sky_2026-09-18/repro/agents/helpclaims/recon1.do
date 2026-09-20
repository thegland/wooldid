clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

di "SENTINEL_VERSION"
di "stata_version = " c(stata_version)
cap which hdidregress
di "hdidregress rc = " _rc
cap which csdid
di "csdid rc = " _rc
cap which did_imputation
di "did_imputation rc = " _rc

* ---- tiny staggered panel: i = onset cohort (5 levels), 12 periods ----
clear
set obs 5
gen i = _n
gen ttre = .
replace ttre = 5 if i==1
replace ttre = 7 if i==2
replace ttre = 6 if i==3
replace ttre = 9 if i==4
* i==5 never treated
gen sg = cond(i<=2, 1, 2)
expand 12
bysort i: gen t = _n
expand 6
gen double u = rnormal()
gen double y = 1 + 0.3*t + 0.5*i + u + 2*( !missing(ttre) & t>=ttre )
gen double wt = 1 + mod(_n,3)
gen long fwv = 1 + mod(_n,3)
gen byte cl = mod(_n,8)+1

di "PRESTATE level=" c(level) " emptycells=" c(emptycells)
local rng0 = c(rngstate)
di "PRESTATE sortedby = `: sortedby'"
qui datasignature
di "PRESTATE datasig = " r(datasignature)

* ================= RUN A =================
di "SENTINEL_A_START"
cap noi wooldid y i t ttre, att att_i subgroups(sg) esprelength(2) espostlength(2) esfixedbaseperiod gteffects
di "A rc = " _rc
di "SENTINEL_A_EB"
cap noi matrix list e(b), format(%12.6f)
di "SENTINEL_A_ESAMP"
qui count if e(sample)
di "count if e(sample) = " r(N) "   e(N) = " e(N)
di "e(properties) = " e(properties)
di "e(wdid_vcvstate) = " e(wdid_vcvstate)
di "e(cmdline) = " e(cmdline)
di "e(vce)=" e(vce) "  e(vcetype)=" e(vcetype)
di "SENTINEL_A_POSTSTATE"
di "POSTSTATE level=" c(level) " emptycells=" c(emptycells)
di "rngstate unchanged = " ("`rng0'" == "`c(rngstate)'")
di "POSTSTATE sortedby = `: sortedby'"
qui datasignature
di "POSTSTATE datasig = " r(datasignature)
di "SENTINEL_A_ERETURN"
cap noi ereturn list, all

* ================= RUN B: identical rerun, hex compare =================
di "SENTINEL_B_START"
tempname bA
matrix `bA' = e(b)
cap noi wooldid y i t ttre, att att_i subgroups(sg) esprelength(2) espostlength(2) esfixedbaseperiod gteffects
di "B rc = " _rc
tempname bB
matrix `bB' = e(b)
cap noi mata: st_local("same", strofreal(all(st_matrix("`bA'") :== st_matrix("`bB'"))))
di "RERUN BITWISE IDENTICAL e(b) = `same'"

* ================= RUN C: att_i only, ES names =================
di "SENTINEL_C_START"
cap noi wooldid y i t ttre, att_i esprelength(2) espostlength(2) esfixedbaseperiod
di "C rc = " _rc
di "C e(b) colnames:"
cap noi di "`: colnames e(b)'"

* ================= RUN D: fweight =================
di "SENTINEL_D_START"
cap noi wooldid y i t ttre [fw=fwv], att
di "D rc = " _rc
di "SENTINEL_D_IW"
cap noi wooldid y i t ttre [iw=wt], att
di "D_iw rc = " _rc

* ================= RUN E: aweight, main/pretreat names only =================
di "SENTINEL_E_START"
cap noi wooldid y i t ttre [aw=wt], att att_it att_i att_itime
di "E rc = " _rc
cap noi di "E e(b) colnames: `: colnames e(b)'"

di "SENTINEL_DONE"
