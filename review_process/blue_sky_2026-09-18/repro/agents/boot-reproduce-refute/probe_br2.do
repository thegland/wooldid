* =================================================================
* probe_br2.do -- SKEPTIC repro of BR-2 (reps() resolution caution is
* unreachable without cluster()) AND of SK-1 (the comma refusal's own
* suggested remedy writes a silently mis-named file).
* Self-contained. Run:
*   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -e do probe_br2.do
* =================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 424242
set sortseed 424242

set obs 40
gen long unit = _n
gen int  ttre = .
replace  ttre = 2004 if unit <= 13
replace  ttre = 2006 if unit > 13 & unit <= 26
gen double ufe = rnormal()
expand 10
bysort unit: gen int year = 1999 + _n
gen double te  = cond(missing(ttre), 0, cond(year >= ttre, 0.5, 0))
gen double y   = 1 + ufe + 0.12*(year-1999) + te + rnormal()*0.7

dis "===== A: CLUSTERED, reps(12): caution expected ====="
cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(12) rseed(11) dots(0)
dis "B2_A_RC=" _rc

dis "===== B: UNCLUSTERED, reps(12): same coarse grid, claim says silent ====="
cap noi wooldid y unit year ttre, att bootstrap(standard) reps(12) rseed(11) dots(0)
dis "B2_B_RC=" _rc

dis "===== C: UNCLUSTERED + boottstat reps(4): floor 2/5=.4, nothing can reject ====="
cap noi wooldid y unit year ttre, att bootstrap(standard) boottstat reps(4) rseed(11) dots(0)
dis "B2_C_RC=" _rc

dis "===== D: UNCLUSTERED + bootdisplay(wald) reps(6) ====="
cap noi wooldid y unit year ttre, att bootstrap(standard) bootdisplay(wald) reps(6) rseed(11) dots(0)
dis "B2_D_RC=" _rc

dis "===== E: permute() sibling caution with permreps(20), NO cluster-gating claim ====="
cap noi wooldid y unit year ttre, att cluster(unit) permute(shuffle) permreps(20) rseed(11) dots(0)
dis "B2_E_RC=" _rc

dis "===== SK1: follow the comma refusal's own printed remedy ====="
cap erase "skb2.dta"
cap erase "skb2  replace.dta"
cap erase "skb2 replace.dta"
cap noi wooldid y unit year ttre, att cluster(unit) bootstrap(standard) reps(20) rseed(11) bootsave(skb2  replace) replace dots(0)
dis "SK1_RC=" _rc
cap confirm file "skb2.dta"
dis "SK1_PLAIN_RC=" _rc
cap confirm file "skb2 replace.dta"
dis "SK1_SPACE1_RC=" _rc
cap confirm file "skb2  replace.dta"
dis "SK1_SPACE2_RC=" _rc

dis "PROBE_BR2_DONE"
