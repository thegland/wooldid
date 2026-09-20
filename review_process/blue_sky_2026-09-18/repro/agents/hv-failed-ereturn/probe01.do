*=====================================================================
* probe01.do -- HV-FAILED-ERETURN
* Q1: what does ORDINARY Stata do to e() when an estimation fails?
* Q2: does a failed wooldid leave reghdfe's internal model in e()?
* Q3: which wooldid failure classes do / do not?
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

*---------------------------------------------------------------
* e() fingerprint helper
*---------------------------------------------------------------
prog def efp
  syntax , tag(string)
  local kb = .
  cap local kb = colsof(e(b))
  local nm ""
  cap local nm : colnames e(b)
  local nm3 = substr("`nm'",1,70)
  local mac : e(macros)
  local haswdid = 0
  if strpos("`mac'","wdid") > 0 local haswdid 1
  if strpos("`mac'","estimand") > 0 local haswdid 1
  local ns = .
  cap count if e(sample)
  cap local ns = r(N)
  dis "FP|`tag'|cmd=" "`e(cmd)'" "|kb=`kb'|N=" e(N) "|nsamp=`ns'|wdidmac=`haswdid'"
  dis "FP|`tag'|cmdline=" substr(`"`e(cmdline)'"',1,90)
  dis "FP|`tag'|names=" "`nm3'"
end

*===============================================================
* PART 1 -- baseline Stata semantics
*===============================================================
dis "=== PART1 BUILTIN ==="
sysuse auto, clear
regress price mpg
efp, tag(P1good)
cap noisily regress price nosuchvarhere
dis "P1RC_builtin=" _rc
efp, tag(P1builtinfail)

dis "=== PART1 ADO WRAPPER ==="
prog def zzfail, eclass
  regress price weight foreign
  error 322
end
regress price mpg
efp, tag(P1good2)
cap noisily zzfail
dis "P1RC_ado=" _rc
efp, tag(P1adofail)

dis "=== PART1 ADO WRAPPER WITH _estimates hold ==="
prog def zzfail2, eclass
  tempname h
  _estimates hold `h', restore nullok
  regress price weight foreign
  error 322
end
regress price mpg
cap noisily zzfail2
dis "P1RC_ado2=" _rc
efp, tag(P1adoheld)

*===============================================================
* PART 2 -- my own staggered panel
*===============================================================
dis "=== PART2 DATA ==="
clear
set obs 48
gen int unit = _n
gen byte grp  = mod(_n,4)
gen int onset = .
replace onset = 2004 if grp==1
replace onset = 2006 if grp==2
replace onset = 2008 if grp==3
gen double a_i = rnormal()
gen double dose = 0
replace dose = 0.5 + runiform() if !missing(onset)
gen byte state = 1 + mod(_n,6)
expand 10
bys unit: gen int year = 1999 + _n
gen byte post = !missing(onset) & year >= onset
gen double tau = 0.8*post*(year - onset + 1)
gen double x1 = rnormal()
gen double y = a_i + 0.10*(year-2000) + tau + 0.4*x1 + 0.30*rnormal()
* dose variable that is zero outside treated periods (untreatedZero shape)
gen double dosep = dose*post
label var y "outcome"
sum y

*---------------------------------------------------------------
* GOOD RUN
*---------------------------------------------------------------
dis "=== PART2 GOOD WOOLDID ==="
cap noisily wooldid y unit year onset, att
dis "P2RC_good=" _rc
efp, tag(P2good)
scalar goodatt = _b[att]
dis "P2_goodatt=" goodatt
mat GB = e(b)
mat list GB

*---------------------------------------------------------------
* FAILURE A: parse-stage (misspelled option)
*---------------------------------------------------------------
dis "=== FAILURE A parse ==="
cap noisily wooldid y unit year onset, attt
dis "RC_A=" _rc
efp, tag(A_parse)

* restore a good run before each subsequent failure test
qui wooldid y unit year onset, att

*---------------------------------------------------------------
* FAILURE B: recon's suggested repro (contreat, i = unit)
*---------------------------------------------------------------
dis "=== FAILURE B contreat untreatedZero i=unit ==="
cap noisily wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)
dis "RC_B=" _rc
efp, tag(B_contreat)

qui wooldid y unit year onset, att

*---------------------------------------------------------------
* FAILURE C: ccc_absorb with unit-level i (documented late error)
*---------------------------------------------------------------
dis "=== FAILURE C ccc_absorb ==="
cap noisily wooldid y unit year onset, att ccc_absorb(x1)
dis "RC_C=" _rc
efp, tag(C_cccabsorb)

qui wooldid y unit year onset, att

*---------------------------------------------------------------
* FAILURE D: fe() that is collinear with a treated cell
*   fe(cellid) where cellid uniquely identifies one treated cohort-period
*---------------------------------------------------------------
dis "=== FAILURE D fe collinear with a treated cell ==="
gen int cellid = 1
replace cellid = 2 if grp==1 & year==2004
cap noisily wooldid y unit year onset, att fe(cellid)
dis "RC_D=" _rc
efp, tag(D_fecollin)

qui wooldid y unit year onset, att

*---------------------------------------------------------------
* FAILURE E: control variable collinear with a treated cell dummy
*---------------------------------------------------------------
dis "=== FAILURE E control collinear ==="
gen byte ccol = (grp==2 & year==2006)
cap noisily wooldid y unit year onset, att controls(ccol)
dis "RC_E=" _rc
efp, tag(E_ctrlcollin)

qui wooldid y unit year onset, att

*---------------------------------------------------------------
* FAILURE F: poisson path (ppmlhdfe) with a negative outcome -> refused early?
*   and poisson with ccc_absorb -> late error
*---------------------------------------------------------------
dis "=== FAILURE F poisson ccc_absorb ==="
gen double ypos = exp(y/4)
cap noisily wooldid ypos unit year onset, att poisson ccc_absorb(x1)
dis "RC_F=" _rc
efp, tag(F_poisson)

qui wooldid y unit year onset, att

*---------------------------------------------------------------
* FAILURE G: unconditionalse path (plain regress) with ccc_absorb
*---------------------------------------------------------------
dis "=== FAILURE G unconditionalse ==="
cap noisily wooldid y unit year onset, att unconditionalse ccc_absorb(x1)
dis "RC_G=" _rc
efp, tag(G_uncond)

dis "=== PROBE01 DONE ==="
