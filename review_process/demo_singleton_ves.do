* demo_singleton_ves.do — SI-12(i): when does the "singleton-cell V_es" phenomenon
* arise, and what does it look like?
*
* WHEN IT ARISES: whenever `i' is a UNIT-level identifier, so each treated
* cohort-by-period cell contains exactly one observation — which is wooldid's most
* ordinary binary usage. In the saturated regression each singleton treated cell is
* fit EXACTLY by its own dummy, so treated residuals are identically zero and the
* cluster-robust V for the event-study block is built from CONTROL-side residual
* variation only. The result is a structurally different variance estimator from the
* pooled event study's: same point estimates (b and att match a pooled-cohort
* formulation at ~1e-15), but ES standard errors that differ by roughly 0.84-1.07x
* on measured designs. It matters mostly for consumers of the full ES covariance —
* honestdid and pretrendpower handoffs — where "which V estimator is this?" affects
* the downstream analysis. (This is the oracle O16 finding; it is NOT the Group A
* refusal case — nothing here is invalid, just structurally different, hence a
* documentation caveat rather than a screen.)
*
* The demo runs the SAME common-onset data twice: (1) unit-level i (singleton
* cells), (2) the pooled single-cohort formulation (icoh = 1 for all treated, the
* reghdfe-oracle trick), both clustered on the unit. Compare: identical b, visibly
* different ES standard errors.
*
* Run from the repo root:   "<stata>" -e do review_process/demo_singleton_ves.do
* (or pass the wooldid.ado path as argument 1). Deterministic; ~20 seconds.

clear all
set more off
set maxvar 20000
set graphics off
set linesize 255

local ADO `1'
if `"`ADO'"' == "" local ADO wooldidinstall/wooldid.ado
if substr(`"`ADO'"', -4, .) != ".ado" local ADO `ADO'/wooldid.ado
do `"`ADO'"'

* ---- common-onset panel: 40 units, onset 7 for units 1..20, t = 1..12 ----
clear
set seed 20260824
qui set obs 40
gen double i = _n
gen double a = rnormal()
gen double ttre = cond(i <= 20, 7, .)
qui expand 12
bysort i: gen double t = _n
gen double y = a + 0.2*t + 0.5*cos(t) + cond(ttre < . & t >= ttre, 0.5 + 0.06*(t - ttre), 0) + rnormal()*0.4
gen double icoh = cond(ttre < ., 1, i)
sort i t
tempfile d
save `d'

* ---- Run 1: unit-level i (every treated cell a singleton) ----
qui wooldid y i t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod
tempname B1 V1
matrix `B1' = e(b_es_att)
matrix `V1' = e(V_es_att)
scalar att1 = _b[main_att]

* ---- Run 2: pooled single treated cohort (multi-unit cells), same clustering ----
use `d', clear
qui wooldid y icoh t ttre, att cluster(i) esprelength(4) espostlength(4) esfixedbaseperiod
tempname B2 V2
matrix `B2' = e(b_es_att)
matrix `V2' = e(V_es_att)
scalar att2 = _b[main_att]

dis _n as result "---- singleton-cell vs pooled-cohort event study, identical data ----"
dis as text "att (unit-level i) = " %12.8f att1 "   att (pooled cohort) = " %12.8f att2 "   reldif = " %8.1e reldif(att1, att2)
dis as text ""
dis as text "  horizon    b (unit i)    b (pooled)      SE (unit i)   SE (pooled)   SE ratio"
local lbl : colnames `B1'
forval j = 1/8 {
  local l : word `j' of `lbl'
  local s1 = sqrt(`V1'[`j', `j'])
  local s2 = sqrt(`V2'[`j', `j'])
  dis as text "  " %-8s "`l'" %12.6f `B1'[1,`j'] %14.6f `B2'[1,`j'] %15.6f `s1' %14.6f `s2' %11.3f `s1'/`s2'
}
dis as text ""
dis as text "Point estimates agree to machine precision; the SEs differ because the"
dis as text "unit-level formulation's treated residuals are identically zero (each"
dis as text "singleton cell is fit exactly), so its clustered V_es is built from"
dis as text "control-side variation only. Neither is 'wrong' for the headline att; the"
dis as text "caveat matters when e(b_es_att)/e(V_es_att) feed honestdid or"
dis as text "pretrendpower, where the two formulations hand those tools different"
dis as text "covariance matrices for the same event study."

dis _n as result "DEMO_SINGLETON_DONE"
exit, clear
