* SELF-CONTAINED REPRO — fw-weights/FW-01
* A weight written as an EXPRESSION whose top-level operator is + or - is pasted verbatim
* into wooldid's estimand-weight arithmetic and binds WRONG:
*   ado:1619  qui gegen `w_per_it' = total(`esamp' * `w') , by(`i' `t')
*   ado:1630  qui gen double `w_att_i'     = `w' `customweighttext'/`w_per_i'
*   ado:1637  qui gen double `w_att_itime' = `w' `customweighttext' /(`w_per_it'*`t_per_i')
*   ado:1644  qui gen double `w_att_it'    = `w' `customweighttext'/`w_per_it'
*   ado:1650  qui gen double `w_att'       = `w' `customweighttext'
* With [aw = wf + wm] line 1644 becomes   wf + wm/W_it   instead of  (wf+wm)/W_it.
* `att' has no denominator so it is EXACT; att_it, att_i and att_itime are silently wrong.
*
* On the DEFAULT estimation path reghdfe refuses the expression first
* (reghdfe.Estimate runs `unab exp : wf+wm, min(1) max(1)' -> "+ invalid name", r(198)),
* so the bug is masked.  customvce() and unconditionalse estimate with plain `regress',
* which accepts any weight expression, and there the wrong numbers are reported with no
* error and no caution.  customvce() needs no safety(off).
*
* wf+wm is constructed to equal wtot EXACTLY, so the two runs below must be identical.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 606
set sortseed 606

clear
set obs 24
gen i = _n
gen ttre = .
replace ttre = 4 if i <= 6
replace ttre = 5 if i > 6  & i <= 12
replace ttre = 6 if i > 12 & i <= 18
gen double ufe = rnormal()
expand 8
bysort i: gen t = _n
expand 3
gen byte post = (t >= ttre) & !missing(ttre)
gen double y = ufe + 0.15*t + 1.25*post + rnormal()*0.5
* the weight the researcher wants: female headcount + male headcount
gen double wf = 10 + int(runiform()*40)
gen double wm = 10 + int(runiform()*40)
gen double wtot = wf + wm
assert wtot == wf + wm

dis "=== 0. the default (reghdfe) path masks the bug with reghdfe's own error ==="
cap noi wooldid y i t ttre [aw=wf+wm], att att_i att_it att_itime
dis "DEFAULTPATH_RC=" _rc
dis "  (reghdfe.Estimate: unab exp : wf+wm, min(1) max(1))"
cap noi qui regress y post i.i i.t [aw=wf+wm]
dis "PLAIN_REGRESS_ACCEPTS_IT_RC=" _rc

dis "=== 1. customvce(robust): [aw=wtot], a single variable ==="
wooldid y i t ttre [aw=wtot], att att_i att_it att_itime customvce(robust)
matrix B1 = e(b)

dis "=== 2. customvce(robust): [aw=wf+wm], the IDENTICAL weight ==="
wooldid y i t ttre [aw=wf+wm], att att_i att_it att_itime customvce(robust)
matrix B2 = e(b)

dis "=== 3. the two must be identical ==="
mata:
A = st_matrix("B1"); B = st_matrix("B2")
nm = st_matrixcolstripe("B1")[.,2]'
for (j=1; j<=cols(A); j++) {
  printf("%-18s  [aw=wtot] %18.12g   [aw=wf+wm] %18.12g   reldiff %10.3e\n",
     nm[j], A[1,j], B[1,j], abs(A[1,j]-B[1,j])/abs(A[1,j]))
}
end

dis "=== 4. the same on the unconditionalse path ==="
qui wooldid y i t ttre [aw=wtot], att att_i att_it att_itime unconditionalse robust safety(off)
matrix C1 = e(b)
qui wooldid y i t ttre [aw=wf+wm], att att_i att_it att_itime unconditionalse robust safety(off)
matrix C2 = e(b)
mata:
A = st_matrix("C1"); B = st_matrix("C2")
nm = st_matrixcolstripe("C1")[.,2]'
for (j=1; j<=cols(A); j++) {
  printf("UNCOND %-18s [aw=wtot] %18.12g   [aw=wf+wm] %18.12g   reldiff %10.3e\n",
     nm[j], A[1,j], B[1,j], abs(A[1,j]-B[1,j])/abs(A[1,j]))
}
end
dis "REPRO_DONE"
