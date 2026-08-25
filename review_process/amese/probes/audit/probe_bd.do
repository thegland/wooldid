* probe_bd — CODE-PATH AUDIT: regimes A/B/D on ONE saved d1 realization.
* B: byCohort+cluster(i) total blackout — model-level vs margins-level vs wooldid-level.
* A: basic+cluster(i) same-file contrast (att fine, ame rank-1) + stale-ado sanity.
* D: byTime+cluster(i) — faithfulness of e(V_es_ame) = J e(V) J'.
clear all
set more off
set maxvar 20000
set graphics off
set linesize 255
local AUD /Users/hegland/.claude/jobs/640b9294/tmp/amese/audit
do "/Users/hegland/Documents/Research/wooldid_dev/test/suite/mc/amese_cases.doh"
do "`AUD'/mkvars.doh"

set seed 60000001
amese_mkdata d1
save "`AUD'/d1_r1.dta", replace

do "/Users/hegland/Documents/Research/wooldid_dev/wooldidinstall/wooldid.ado"

*========== PART B: wooldid on saved file, byCohort ci — e() surface (task 7) ==========
use "`AUD'/d1_r1.dta", clear
cap noi wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(byCohort) cluster(i)
dis "B_WOOLDID|rc=" _rc
foreach m in b V b_es_ame V_es_ame b_es_att V_es_att {
  cap confirm matrix e(`m')
  if _rc == 0 dis "B_E|`m'|EXISTS|" rowsof(e(`m')) "x" colsof(e(`m'))
  else dis "B_E|`m'|ABSENT"
}
matrix BWD = e(wdidmainresults)

*========== PART A: wooldid on the SAME file, basic ci ==========
use "`AUD'/d1_r1.dta", clear
cap noi wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(basic) cluster(i)
dis "A_WOOLDID|rc=" _rc
foreach m in b V b_es_ame V_es_ame b_es_att V_es_att {
  cap confirm matrix e(`m')
  if _rc == 0 dis "A_E|`m'|EXISTS|" rowsof(e(`m')) "x" colsof(e(`m'))
  else dis "A_E|`m'|ABSENT"
}
matrix AWD = e(wdidmainresults)
cap matrix AVES = e(V_es_ame)
if _rc == 0 {
  mata: ev = symeigenvalues(st_matrix("AVES")); printf("A_VESAME|eig_max=%12.0g|eig2=%12.0g|diag1=%15.0g\n", max(ev), sort(ev',1)[cols(ev)-1], st_matrix("AVES")[1,1])
}
dis "A same-file table (compare att SEs vs B blackout):"
matrix list AWD, format(%12.0g)

*========== PART D: wooldid on the SAME file, byTime ci — capture V_es_ame ==========
use "`AUD'/d1_r1.dta", clear
cap noi wooldid y i t ttre, att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod contreatcontrols(byTime) cluster(i)
dis "D_WOOLDID|rc=" _rc
cap confirm matrix e(V_es_ame)
dis "D_E|V_es_ame|" cond(_rc==0, "EXISTS", "ABSENT")
matrix DVES = e(V_es_ame)
matrix DBES = e(b_es_ame)
matrix DWD = e(wdidmainresults)

*========== PART B manual: reghdfe byCohort absorb, NOISILY ==========
use "`AUD'/d1_r1.dta", clear
mkvars
set emptycells drop
noi reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.i#c.(ringer)) cluster(i) tol(1e-9)
dis "B_MODEL|df_r=" e(df_r) "|df_m=" e(df_m) "|df_a=" e(df_a) "|rank=" e(rank) "|Nclust=" e(N_clust) "|N=" e(N)
local dfr = e(df_r)
mata:
V = st_matrix("e(V)")
nm = st_matrixcolstripe("e(V)")
K  = rows(V)
d  = diagonal(V)
isslope = strpos(nm[.,2], "dose") :> 0
isom    = strpos(nm[.,2], "o.") :> 0
npos=0; nzero=0; nneg=0; nmiss=0
for (j=1; j<=K; j++) {
  if (missing(d[j])) nmiss++
  else if (d[j]==0)  nzero++
  else if (d[j]<0)   nneg++
  else               npos++
}
printf("B_EVDIAG|K=%g|pos=%g|zero=%g|neg=%g|miss=%g\n", K, npos, nzero, nneg, nmiss)
printf("B_EVMISS|total missing entries=%g of %g\n", sum(rowmissing(V)), K*K)
// per-class diag summary (exclude omitted from the class mins)
live = 1 :- isom
sl = select(d, isslope :& live); du = select(d, (1:-isslope) :& live)
printf("B_EVDIAG|live slope coefs n=%g min=%21.15g max=%21.15g nmiss=%g\n", rows(sl), min(sl), max(sl), colmissing(sl))
printf("B_EVDIAG|live dummy coefs n=%g min=%21.15g max=%21.15g nmiss=%g\n", rows(du), min(du), max(du), colmissing(du))
printf("B_EVDIAG|omitted coefs n=%g\n", sum(isom))
// names of any missing-diag coefficients (up to 30)
cnt=0
for (j=1; j<=K; j++) {
  if (missing(d[j]) & cnt<30) { printf("B_EVMISSNAME|%s\n", nm[j,2]); cnt++ }
}
end

* binary att margins over ES window (exact wooldid L1771 form)
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "B_MARG_ATT|rc=" _rc
matrix RTb = r(table)
matrix RVb = r(V)
matrix RJb = r(Jacobian)
dis "B_MARG_ATT r(table):"
matrix list RTb, format(%14.0g)
mata:
d = diagonal(st_matrix("RVb"))
printf("B_MARG_ATT|rVdiag nmiss=%g nneg=%g\n", colmissing(d), sum((d:<0):*(1:-(d:>=.))))
printf("B_MARG_ATT rV diag raw:\n"); d'
JV = st_matrix("RJb")*st_matrix("e(V)")*st_matrix("RJb")'
printf("B_MARG_ATT J eV Jt diag:\n"); diagonal(JV)'
printf("B_MARG_ATT|max|rV-JVJt|=%12.0g\n", max(abs(st_matrix("RVb")-JV)))
end

* continuous ame margins over ES window (exact wooldid L2186 form)
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr')
dis "B_MARG_AME|rc=" _rc
matrix RTa = r(table)
matrix RVa = r(V)
matrix RJa = r(Jacobian)
matrix RBa = r(b)
dis "B_MARG_AME r(table):"
matrix list RTa, format(%14.0g)
mata:
d = diagonal(st_matrix("RVa"))
printf("B_MARG_AME|rVdiag nmiss=%g nneg=%g\n", colmissing(d), sum((d:<0):*(1:-(d:>=.))))
printf("B_MARG_AME rV diag raw:\n"); d'
JV = st_matrix("RJa")*st_matrix("e(V)")*st_matrix("RJa")'
printf("B_MARG_AME J eV Jt diag:\n"); diagonal(JV)'
end
* match manual ame b to wooldid's displayed ES_ame b (rows 13..20 of BWD)
mata:
BWD = st_matrix("BWD")
RB  = st_matrix("RBa")
if (rows(BWD)>=20 & cols(RB)>=8) printf("B_MATCH|max|manual_b - wooldid_b|=%12.0g\n", max(abs(BWD[13::20,1] - RB[1,1::8]')))
end

* B determinism probe: a second d1 realization, flags only
set seed 60000002
amese_mkdata d1
mkvars
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.i#c.(ringer)) cluster(i) tol(1e-9)
local dfr2 = e(df_r)
qui margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr2')
mata: st_local("m1", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfr2')
mata: st_local("m2", strofreal(colmissing(st_matrix("r(table)")[2,.]')))
dis "B_SEED2|attSEmiss=`m1'|ameSEmiss=`m2'"

*========== PART B-regress: same spec via plain regress (CRVE-generic vs reghdfe?) ==========
use "`AUD'/d1_r1.dta", clear
mkvars
set emptycells drop
cap qui regress y i.i i.t i.i#c.ringer i._i#i._t#i.targetflag##c.dose, vce(cluster i)
dis "BREG|rc=" _rc "|df_r=" e(df_r) "|rank=" e(rank)
local dfrr = e(df_r)
mata:
d = diagonal(st_matrix("e(V)"))
printf("BREG_EVDIAG|K=%g|miss=%g|zero=%g|neg=%g\n", rows(d), colmissing(d), sum(d:==0), sum((d:<0):*(1:-(d:>=.))))
end
cap noi margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrr')
dis "BREG_MARG_AME|rc=" _rc
matrix RTr = r(table)
dis "BREG_MARG_AME r(table) se row:"
mata: st_matrix("RTr")[1::2,.]
cap noi margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrr')
dis "BREG_MARG_ATT|rc=" _rc
mata: st_matrix("r(table)")[1::2,.]

*========== PART A manual on same file: basic absorb — contrast ==========
use "`AUD'/d1_r1.dta", clear
mkvars
set emptycells drop
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.null##c.(ringer)) cluster(i) tol(1e-9)
dis "A_MODEL|df_r=" e(df_r) "|rank=" e(rank)
local dfra = e(df_r)
mata:
V = st_matrix("e(V)")
nm = st_matrixcolstripe("e(V)")
d = diagonal(V)
isslope = strpos(nm[.,2], "dose") :> 0
isom    = strpos(nm[.,2], "o.") :> 0
live = 1 :- isom
sl = select(d, isslope :& live); du = select(d, (1:-isslope) :& live)
printf("A_EVDIAG|K=%g|miss=%g|zero=%g|neg=%g\n", rows(d), colmissing(d), sum(d:==0), sum((d:<0):*(1:-(d:>=.))))
printf("A_EVDIAG|live slope n=%g min=%21.15g max=%21.15g\n", rows(sl), min(sl), max(sl))
printf("A_EVDIAG|live dummy n=%g min=%21.15g max=%21.15g\n", rows(du), min(du), max(du))
end
qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfra')
mata:
d = diagonal(st_matrix("r(V)"))
printf("A_MARG_AME rV diag raw:\n"); d'
printf("A_MARG_AME|nmiss=%g|min=%21.15g|max=%21.15g\n", colmissing(d), min(d), max(d))
end
qui margins [iw=1], dydx(targetflag) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfra')
mata:
d = diagonal(st_matrix("r(V)"))
printf("A_MARG_ATT rV diag raw:\n"); d'
printf("A_MARG_ATT|nmiss=%g|min=%21.15g|max=%21.15g\n", colmissing(d), min(d), max(d))
end

*========== PART D manual: byTime absorb — faithfulness ==========
use "`AUD'/d1_r1.dta", clear
mkvars
set emptycells drop
qui reghdfe y i._i#i._t#i.targetflag##c.dose, absorb(i.i i.t i.t#c.(ringer)) cluster(i) tol(1e-9)
dis "D_MODEL|df_r=" e(df_r) "|rank=" e(rank)
local dfrd = e(df_r)
qui margins [iw=1], dydx(dose) subpop(if targetflag==1 & _t >= $HS & _t <= $HE) noestimcheck over(_t) df(`dfrd')
matrix DRV = r(V)
matrix DRB = r(b)
matrix DRJ = r(Jacobian)
mata:
RV = st_matrix("DRV")
JV = st_matrix("DRJ")*st_matrix("e(V)")*st_matrix("DRJ")'
printf("D_FAITH|max|rV - J eV Jt| = %12.0g\n", max(abs(RV-JV)))
DV = st_matrix("DVES")
printf("D_FAITH|dim e(V_es_ame)=%g|dim manual rV=%g\n", rows(DV), rows(RV))
if (rows(DV)==rows(RV)) printf("D_FAITH|max|e(V_es_ame) - manual rV| = %12.0g\n", max(abs(DV-RV)))
DB = st_matrix("DBES")
RB = st_matrix("DRB")
if (cols(DB)==cols(RB)) printf("D_FAITH|max|e(b_es_ame) - manual rb| = %12.0g\n", max(abs(DB-RB)))
ev = symeigenvalues(RV)
printf("D_RANK|eigmax=%12.0g|eigmin=%12.0g\n", max(ev), min(ev))
end
dis "PROBE_BD_DONE"
