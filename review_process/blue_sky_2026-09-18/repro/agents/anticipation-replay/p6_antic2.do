*==============================================================================
* p6_antic2.do -- anticipation() follow-ups
*  (a) subgroups(ttre): sthlp:826 says labels stay the RECORDED onset, which
*      contradicts the same paragraph's "exactly equivalent to replacing ttre
*      with ttre - #".  Measure which one the code does.
*  (b) semantic check: with anticipation(k), relative period 0 must be the
*      calendar period ttre-k.  Read it off e(gteffects).
*  (c) k=0 no-op and k<0 refusal (also in p2, repeated self-contained).
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/anticipation-replay"
which wooldid
local dN 100
local dT 14
local donsets "2005 2008 2010"
include "`W'/mkdata.do"
save "`W'/d6.dta", replace
include "`W'/dumpe.do"

* (a) subgroups(ttre)
qui use "`W'/d6.dta", clear
qui gen double ttrelab = cond(ttre < ., ttre, 9999)
cap noi wooldid y i t ttre, att subgroups(ttrelab) robust espre(2) espost(2) anticipation(2)
dis "SGA_RC=" _rc
if _rc == 0 dumpe , f("`W'/out_c28_A.txt")
if _rc == 0 dis "SGA_sglist=[" e(wdid_sglist) "]"
qui use "`W'/d6.dta", clear
qui replace ttre = ttre - 2
qui gen double ttrelab = cond(ttre < ., ttre, 9999)
cap noi wooldid y i t ttre, att subgroups(ttrelab) robust espre(2) espost(2)
dis "SGB_RC=" _rc
if _rc == 0 dumpe , f("`W'/out_c28_B.txt")
if _rc == 0 dis "SGB_sglist=[" e(wdid_sglist) "]"

* subgroups() pointing at the ONSET VARIABLE ITSELF (what sthlp:826 talks about)
qui use "`W'/d6.dta", clear
cap noi wooldid y i t ttre, att subgroups(ttre) robust espre(2) espost(2) anticipation(2)
dis "SGT_RC=" _rc
if _rc == 0 dis "SGT_sglist=[" e(wdid_sglist) "]"

* (b) semantic: relyr==0 must sit at caltime == recorded onset - k
qui use "`W'/d6.dta", clear
cap noi wooldid y i t ttre, att robust espre(3) espost(3) gteffects anticipation(2)
dis "SEM_RC=" _rc
cap noi {
  tempname G
  mat `G' = e(gteffects)
  local bad = 0
  local n = rowsof(`G')
  forvalues r = 1/`n' {
    if `G'[`r',4] == 0 {
      * relyr == 0: columns are i onset caltime relyr ...
      if `G'[`r',3] != `G'[`r',2] local ++bad
    }
  }
  dis "SEM_rows=" `n' "  relyr0_caltime_ne_onset=" `bad'
  * show a few rows
  forvalues r = 1/6 {
    dis "row `r': i=" `G'[`r',1] " onset=" `G'[`r',2] " caltime=" `G'[`r',3] " relyr=" `G'[`r',4]
  }
}
dis "SEMCHK_RC=" _rc

* the same run with no anticipation, for the recorded-onset reference
qui use "`W'/d6.dta", clear
cap noi wooldid y i t ttre, att robust espre(3) espost(3) gteffects
cap noi {
  tempname G0
  mat `G0' = e(gteffects)
  forvalues r = 1/6 {
    dis "noantic row `r': i=" `G0'[`r',1] " onset=" `G0'[`r',2] " caltime=" `G0'[`r',3] " relyr=" `G0'[`r',4]
  }
}

* (c) k=0 / k<0
qui use "`W'/d6.dta", clear
cap noi wooldid y i t ttre, att robust espre(2) espost(2) anticipation(0)
dis "K0_RC=" _rc
if _rc == 0 dumpe , f("`W'/out_k0_A.txt")
if _rc == 0 dis "K0_e_anticipation=" e(anticipation)
qui use "`W'/d6.dta", clear
cap noi wooldid y i t ttre, att robust espre(2) espost(2)
dis "K0B_RC=" _rc
if _rc == 0 dumpe , f("`W'/out_k0_B.txt")
qui use "`W'/d6.dta", clear
cap noi wooldid y i t ttre, att robust anticipation(-1)
dis "KNEG_RC=" _rc
cap noi wooldid y i t ttre, att robust anticipation(2.5)
dis "KFRAC_RC=" _rc

dis "P6_DONE"
