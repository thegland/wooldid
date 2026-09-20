*=======================================================================
* r1b_consumers.do -- follow-up to r1_sgnames.do.
*   (1) coefplot: is drop(*_s1*) really a SILENT no-op?  (count properly;
*       coefplot's generate(pfx) makes pfx+"b", not pfx+"_b")
*   (2) e(wdidsgresults_1) row names (the table the user is steered to)
*   (3) do the BOOTSTRAP tables really carry the _SG infix too?  The
*       prober's "the fix belongs in the help, not the code" rests on
*       the _SG infix being used CONSISTENTLY everywhere.
* SELF-CONTAINED.
*=======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

prog def dumprows
  args mname tag
  local rn : rownames `mname'
  dis "DUMPROWS_`tag'_NROW=" rowsof(`mname')
  forval j = 1/`=rowsof(`mname')' {
    local nj : word `j' of `rn'
    dis "DR_`tag' `j' nm=<`nj'>"
  }
end
prog def dumpcolsp
  args mname tag
  local cn : colnames `mname'
  dis "DUMPCOLS_`tag'_NCOL=" colsof(`mname')
  forval j = 1/`=colsof(`mname')' {
    local nj : word `j' of `cn'
    dis "DC_`tag' `j' nm=<`nj'>"
  }
end

clear
set obs 90
gen long unit  = _n
gen byte sg    = 1 + mod(_n, 2)
gen byte cls   = mod(_n, 5)
gen int  ttre  = .
replace  ttre  = 2004 if cls == 1
replace  ttre  = 2006 if cls == 2
replace  ttre  = 2008 if cls == 3
gen double ai  = rnormal()
expand 10
bysort unit: gen int year = 2000 + _n - 1
gen byte post  = (ttre < .) & (year >= ttre)
gen double tau = 1.0 + 0.5*(sg == 2)
gen double y   = ai + 0.10*(year-2000) + tau*post + 0.30*rnormal()
compress

qui wooldid y unit year ttre, att att_it subgroups(sg) ///
    esprelength(2) espostlength(2) robust
dis "RUN_RC=" _rc

tempname Bm
matrix `Bm' = e(b)
dumpcolsp `Bm' EB

dis "--- e(wdidsgresults_1) ---"
tempname SG1
matrix `SG1' = e(wdidsgresults_1)
dumprows `SG1' SGRES1

*=================== coefplot: counts, done right ======================
cap drop cpA*
qui coefplot, generate(cpA) nodraw
qui count if cpAb < .
dis "CP_PLAIN_NCOEF=" r(N)
cap drop cpB*
noi coefplot, drop(*_s1*) generate(cpB) nodraw
dis "CP_DROPHELP_RC=" _rc
qui count if cpBb < .
dis "CP_DROPHELP_NCOEF=" r(N)
cap drop cpC*
qui coefplot, drop(*SG1*) generate(cpC) nodraw
qui count if cpCb < .
dis "CP_DROPTRUE_NCOEF=" r(N)
cap drop cpD*
noi coefplot, keep(*_s1*) generate(cpD) nodraw
dis "CP_KEEPHELP_RC=" _rc
cap drop cpE*
qui coefplot, keep(*SG1*) generate(cpE) nodraw
qui count if cpEb < .
dis "CP_KEEPTRUE_NCOEF=" r(N)

*=================== bootstrap tables: _SG or _s ? =====================
cap noi wooldid y unit year ttre, att subgroups(sg) esprelength(2) ///
    espostlength(2) robust bootstrap(standard) reps(8) rseed(4321)
dis "BOOT_RC=" _rc
foreach mm in boot_wald_mainresults bootdist V_boot b_boot {
  tempname M`=strlen("`mm'")'
  cap matrix BTMP = e(`mm')
  if _rc == 0 {
    dis "BOOTMAT `mm' :"
    cap dumpcolsp BTMP BC_`mm'
    cap dumprows BTMP BR_`mm'
    matrix drop BTMP
  }
  else dis "BOOTMAT `mm' ABSENT rc=" _rc
}
foreach en in boot_wald_sgresults_1 boot_wald_sgresults_s1 ///
              boot_pctl_sgresults_1 {
  cap confirm matrix e(`en')
  dis "BOOTSG `en' RC=" _rc
}

dis "R1B_DONE"
