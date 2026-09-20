*! r04b — sharpen/limit FCW-04.  The reps() resolution caution sits inside
*! `if "`clus'" != "" {' (ado:6459-6501).  But the DISCLOSURE FOOTER at ado:11784-11785
*! is NOT cluster-gated, so on the percentile and boottstat legs the resolution is still
*! stated.  ado:11785 is `else if "`wdbdsp'" != "wald"', so under bootdisplay(wald) the
*! footer says nothing either.  Does ANY line tell an unclustered wald run that its
*! joint/equality/adversarial p-values cannot fall below 1/(reps+1)?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 5551
set sortseed 5551
set obs 12
gen long cty = _n
gen onset = .
replace onset = 4 if cty <= 4
replace onset = 6 if cty > 4 & cty <= 8
expand 4
sort cty
gen long unit = _n
gen double uf = rnormal(0,1)
expand 7
bysort unit (cty): gen int t = _n
gen byte D = (onset < . & t >= onset)
gen double y = 2 + uf + .1*t + 1*D + rnormal(0,1)
keep y unit cty t onset

dis "===== W1: jointtests + bootdisplay(wald) + reps(10), NO cluster()"
cap noi wooldid y unit t onset, att att_i jointtests bootstrap(standard) reps(10) rseed(9) bootdisplay(wald) dots(0)
dis "W1_RC=" _rc
dis "W1 joint p = " e(boot_p_joint_att)

dis "===== W2: identical plus cluster(cty)  (the caution SHOULD print here)"
cap noi wooldid y unit t onset, att att_i jointtests cluster(cty) bootstrap(standard) reps(10) rseed(9) bootdisplay(wald) dots(0)
dis "W2_RC=" _rc
dis "===== DONE"
