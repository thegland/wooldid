*=====================================================================
* probe02.do -- HV-FAILED-ERETURN
* Map which wooldid failure classes leave reghdfe's internal model in e().
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
which wooldid

prog def efp
  syntax , tag(string)
  local kb = .
  cap local kb = colsof(e(b))
  local nm ""
  cap local nm : colnames e(b)
  local nm3 = substr("`nm'",1,60)
  local ns = .
  cap count if e(sample)
  cap local ns = r(N)
  dis "FP|`tag'|cmd=" "`e(cmd)'" "|kb=`kb'|N=" e(N) "|nsamp=`ns'"
  dis "FP|`tag'|cmdline=" substr(`"`e(cmdline)'"',1,80)
  dis "FP|`tag'|names=" "`nm3'"
end

*--------------------- DATA: plain staggered panel --------------------
clear
set obs 60
gen int unit = _n
gen byte grp  = mod(_n,5)
gen double onset = .
replace onset = 2004 if grp==1
replace onset = 2006 if grp==2
replace onset = 2008 if grp==3
gen double a_i = rnormal()
gen byte state = 1 + mod(_n,10)
gen double dose0 = 0.5 + runiform()
expand 10
bys unit: gen int year = 1999 + _n
gen byte post = !missing(onset) & year >= onset
gen double x1 = rnormal()
gen double tau = 0
replace tau = 0.8*(year-onset+1) if post
gen double y = a_i + 0.10*(year-2000) + tau + 0.4*x1 + 0.30*rnormal()
gen double dosep = dose0*post
egen long sy = group(state year)
gen byte sg = mod(unit,2)
count

*--------------------- GOOD RUN ---------------------------------------
dis "=== GOOD ==="
wooldid y unit year onset, att
efp, tag(GOOD)
dis "GOOD_att=" _b[main_att]
dis "GOOD_wdidmac=" "`: e(macros)'"

*--------------------- helper: run a failure then fingerprint ---------
* (written out longhand so each block is visible in the log)

dis "=== A: parse-stage typo ==="
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, attt
dis "RC|A=" _rc
efp, tag(A_parse)

dis "=== B: recon repro contreat untreatedZero i=unit ==="
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, att contreat(dosep) contreatcontrols(untreatedZero)
dis "RC|B=" _rc
efp, tag(B_contreat)

dis "=== C: ccc_absorb with unit-level i ==="
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, att ccc_absorb(x1)
dis "RC|C=" _rc
efp, tag(C_ccc)

dis "=== D: fe(state#year) saturating ==="
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, att fe(sy)
dis "RC|D=" _rc
efp, tag(D_fe)

dis "=== E: control collinear with one treated cell ==="
qui wooldid y unit year onset, att
gen byte ccol = (grp==2 & year==2006)
cap noisily wooldid y unit year onset, att controls(ccol)
dis "RC|E=" _rc
efp, tag(E_ctrl)

dis "=== F: poisson + ccc_absorb ==="
qui wooldid y unit year onset, att
gen double ypos = exp(y/4)
cap noisily wooldid ypos unit year onset, att poisson ccc_absorb(x1)
dis "RC|F=" _rc
efp, tag(F_poi)

dis "=== G: unconditionalse (plain regress path) + ccc_absorb ==="
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, att unconditionalse ccc_absorb(x1)
dis "RC|G=" _rc
efp, tag(G_unc)

dis "=== H: subgroups with an empty treated cell ==="
qui wooldid y unit year onset, att
gen byte sgbad = (unit==1)
cap noisily wooldid y unit year onset, att subgroups(sgbad)
dis "RC|H=" _rc
efp, tag(H_sg)

dis "=== I: interactivecontrols collinear ==="
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, att interactivecontrols(ccol)
dis "RC|I=" _rc
efp, tag(I_ic)

dis "=== J: cluster() on a var not constant within i ==="
qui wooldid y unit year onset, att
cap noisily wooldid y unit year onset, att cluster(year)
dis "RC|J=" _rc
efp, tag(J_clus)

dis "PROBE02 DONE"
