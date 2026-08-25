* bonus2.do — case-6 channel: what does reghdfe itself return for the coarse
* regression under cluster(i)? (Mimic the traced internal form on identical data.)
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
clear
set seed 2010
qui set obs 30
gen i = _n
gen ttre = cond(mod(_n,2)==0, 6, .)
qui expand 6
bysort i: gen unit = _n
gen dose = exp(0.8*runiform() + 0.2*mod(i,5))
qui expand 12
bysort i unit: gen t = _n
gen y = 0.1*i + 0.05*t + rnormal(0,1000)
sort i unit t
* internal cell vars (mimic wooldid: _i = cohort for treated / 0 pooled controls;
* _t = calendar period; targetflag = treated-post): esf espre(2) keeps the full
* sample; approximate the cell structure for the mechanism check
gen _i = cond(ttre < ., i, 0)
gen tf = (ttre < .)
qui reghdfe y i._i#i.t#i.tf i._i#i.t#(c.dose), absorb(i.i i.t) vce(cluster i)
dis "B2 reghdfe rc=" _rc " df_r=" e(df_r)
tempname V
matrix `V' = e(V)
mata: V = st_matrix(st_local("V")); printf("B2 e(V): dim %g, max|V| = %12.6g, exactly zero = %1.0f, n nonzero diag = %g of %g\n", rows(V), max(abs(V)), max(abs(V))==0, sum(diagonal(V) :> 0), rows(V))
dis "B2_DONE"
exit, clear
