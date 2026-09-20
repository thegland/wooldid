* SKEPTIC probe of EC-6 (weights in e()) and EC-7 (subgroup coefficient names)
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
gen byte sg = mod(id,2)
gen double ui = rnormal()
gen int ttre = .
replace ttre = 2004 if mod(id,4)==0
replace ttre = 2007 if mod(id,4)==1
replace ttre = 2010 if mod(id,4)==2
expand 13
bys id: gen int t = 2000 + _n - 1
gen double y = 1 + ui + 0.05*(t-2000) + rnormal()
replace y = y + 1.5 + 0.4*(t-ttre) + 0.8*sg if t>=ttre & ttre<.
* a realistic survey-style sampling weight, unit-constant
bys id (t): gen double pw = exp(0.5*ui[1] + 0.3) * 20
by id: replace pw = pw[1]

local SPEC "att esfixedbaseperiod esprelength(3) espostlength(3)"

*------------- EC-6: what does e() record about weights? -------------
foreach W in "none" "aw" "pw" "awrob" "pwrob" "awclu" "pwclu" {
  if "`W'" == "none"  cap noisily wooldid y id t ttre, `SPEC'
  if "`W'" == "aw"    cap noisily wooldid y id t ttre [aw=pw], `SPEC'
  if "`W'" == "pw"    cap noisily wooldid y id t ttre [pw=pw], `SPEC'
  if "`W'" == "awrob" cap noisily wooldid y id t ttre [aw=pw], `SPEC' robust
  if "`W'" == "pwrob" cap noisily wooldid y id t ttre [pw=pw], `SPEC' robust
  if "`W'" == "awclu" cap noisily wooldid y id t ttre [aw=pw], `SPEC' cluster(clu)
  if "`W'" == "pwclu" cap noisily wooldid y id t ttre [pw=pw], `SPEC' cluster(clu)
  local rc = _rc
  if `rc'==0 {
    mat b = e(b)
    mat V = e(V)
    local bb = b[1,1]
    local ss = sqrt(V[1,1])
    dis "WT_`W' rc=`rc' b=" %21x `bb' "  se=" %21x `ss'
    dis "WT_`W' vce=|" e(vce) "| vcetype=|" e(vcetype) "| wtype=|" e(wtype) "| wexp=|" e(wexp) "| N=" e(N)
    dis "WT_`W' cmdline=|" e(cmdline) "|"
  }
  else dis "WT_`W' rc=`rc'"
}

* does an ordinary Stata command post wtype/wexp on the same data? (the contract check)
cap noisily reg y i.t [pw=pw] if id<=40
dis "REG_pw vce=|" e(vce) "| wtype=|" e(wtype) "| wexp=|" e(wexp) "|"
cap noisily reg y i.t [aw=pw] if id<=40
dis "REG_aw vce=|" e(vce) "| wtype=|" e(wtype) "| wexp=|" e(wexp) "|"

* does esttab / estimates table notice?  and can a .ster tell?
cap noisily wooldid y id t ttre [pw=pw], `SPEC'
cap estimates store WPW
cap noisily wooldid y id t ttre, `SPEC'
cap estimates store WNONE
cap noisily estimates table WPW WNONE, b se stats(N)
dis "ESTTAB_RC=" _rc

*------------- EC-7: subgroup coefficient names -------------
cap noisily wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sg) cluster(clu)
dis "SG_RC=" _rc
dis "SG colnames: `: colnames e(b)'"
foreach nm in main_att_s1 main_att_s0 main_SG1_att main_SG0_att pretreatment_att_s1 pretreatment_SG1_att ES_att_post0_s1 ES_att_SG1_post0 {
  cap lincom `nm'
  dis "LINCOM `nm' rc=" _rc
}
* what does the help literally say, and what does the code build?
dis "SG e(b_es_att_s1) present? "
cap confirm matrix e(b_es_att_s1)
dis "  rc=" _rc
cap confirm scalar e(numpre_att_s1)
dis "  numpre_att_s1 rc=" _rc
dis "ALLDONE_P3"
