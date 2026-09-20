*==============================================================================
* REPRO (skeptic re-derivation of EC-6 and EC-7).
*  EC-6: no e(wtype)/e(wexp) on any weighted run, and e(vce)="ols" under [pw=]
*        although the standard errors are BIT-IDENTICAL to the robust ones.
*  EC-7: the help says subgroup coefficients take the suffix _s#; the code builds
*        the infix _SG<value>_ , so the documented name errors r(111).
*
* Run:  /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do repro_R67_wt_sgnames.do
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 5150
set sortseed 5150
which wooldid

clear
set obs 120
gen long id = _n
gen byte clu = ceil(id/5)
gen byte sg  = mod(id,2)
gen double ui = rnormal()
gen int ttre = .
replace ttre = 2004 if mod(id,4)==0
replace ttre = 2007 if mod(id,4)==1
replace ttre = 2010 if mod(id,4)==2
expand 13
bys id: gen int t = 2000 + _n - 1
gen double y = 1 + ui + 0.05*(t-2000) + rnormal()
replace y = y + 1.5 + 0.4*(t-ttre) + 0.8*sg if t>=ttre & ttre<.
* a unit-constant survey-style sampling weight
bys id (t): gen double pw = exp(0.5*ui[1] + 0.3)*20
by id: replace pw = pw[1]

local SPEC "att esfixedbaseperiod esprelength(3) espostlength(3)"

*---------------------------- EC-6 ----------------------------
dis "--- wooldid ---"
wooldid y id t ttre, `SPEC'
mat V = e(V)
dis "unweighted   se = " %21x sqrt(V[1,1]) "  e(vce)=" e(vce)
wooldid y id t ttre [aw=pw], `SPEC'
mat V = e(V)
dis "[aw=pw]      se = " %21x sqrt(V[1,1]) "  e(vce)=" e(vce) "  e(wtype)=" e(wtype)
wooldid y id t ttre [pw=pw], `SPEC'
mat V = e(V)
dis "[pw=pw]      se = " %21x sqrt(V[1,1]) "  e(vce)=" e(vce) "  e(wtype)=" e(wtype) "  e(vcetype)=" e(vcetype)
wooldid y id t ttre [aw=pw], `SPEC' robust
mat V = e(V)
dis "[aw=pw] rob  se = " %21x sqrt(V[1,1]) "  e(vce)=" e(vce)
dis " ^^ the [pw=] and the [aw=]+robust standard errors are bit-identical, yet"
dis "    e(vce) reads 'ols' on the pweighted run and carries no e(vcetype)."
dis "--- the contract every other estimation command honours ---"
qui reg y i.t [pw=pw]
dis "reg [pw=pw]  e(vce)=" e(vce) "  e(wtype)=" e(wtype) "  e(wexp)=" e(wexp)
qui reg y i.t [aw=pw]
dis "reg [aw=pw]  e(vce)=" e(vce) "  e(wtype)=" e(wtype) "  e(wexp)=" e(wexp)

*---------------------------- EC-7 ----------------------------
wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sg) cluster(clu)
dis "e(b) colnames: `: colnames e(b)'"
foreach nm in main_att_s1 main_SG1_att pretreatment_att_s1 pretreatment_SG1_att ES_att_post0_s1 ES_att_SG1_post0 {
  cap lincom `nm'
  dis "lincom `nm'  rc=" _rc
}
dis "the _s# suffix IS correct for every other subgroup-keyed e() name:"
cap confirm matrix e(b_es_att_s1)
dis "  e(b_es_att_s1)   rc=" _rc
cap confirm matrix e(wdidsgresults_1)
dis "  e(wdidsgresults_1) rc=" _rc
dis "REPRO_R67_DONE"
