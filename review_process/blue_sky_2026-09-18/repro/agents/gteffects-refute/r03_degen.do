* r03: (Q) does fastnose really withhold e(V_es_att)?  its own Caution says it does.
*      (R) the gteffects cell-SE degeneracy guard is EXACT-ZERO only (ado:13496
*          `se = (v < . & v > 0) ? sqrt(v) : .`), so a rank-deficient cluster-robust
*          VCE whose diagonal residue sits at ~1e-29 yields se ~ 1e-15, p ~ 1e-28 and a
*          zero-width confidence interval, printed as an ordinary row.
*          The gteffects help promises "degenerate zero variances are always reported as
*          missing, never zero".
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

cap prog drop sedegen
prog def sedegen
  syntax , tag(string)
  tempname G
  cap matrix `G' = e(gteffects)
  if _rc != 0 {
    dis "SED|`tag'|NOMAT"
    exit
  }
  local cn : colnames `G'
  local js : list posof "se" in cn
  local jp : list posof "p" in cn
  local jb : list posof "estimate" in cn
  local jl : list posof "lb_95ci" in cn
  local ju : list posof "ub_95ci" in cn
  mata: GG = st_matrix("`G'"); ///
        se = GG[.,`js']; pp = GG[.,`jp']; bb = GG[.,`jb']; ll = GG[.,`jl']; uu = GG[.,`ju']; ///
        nse = select(se, se :< .); ///
        st_local("nr3", strofreal(rows(se))); ///
        st_local("nmiss", strofreal(rows(se) - rows(nse))); ///
        st_local("nzero", strofreal(sum(se :== 0))); ///
        st_local("ntiny", strofreal(sum((se :< 1e-8) :& (se :> 0)))); ///
        st_local("minse", strofreal(rows(nse)>0 ? min(nse) : ., "%12.4e")); ///
        st_local("minp",  strofreal(min(select(pp, pp :< .)), "%12.4e")); ///
        st_local("minw",  strofreal(min(select(uu-ll, (uu :< .) :& (ll :< .))), "%12.4e"))
  dis "SED|`tag'|rows=`nr3'|se_missing=`nmiss'|se_exact0=`nzero'|se_lt_1e8=`ntiny'|min_se=`minse'|min_p=`minp'|min_ciwidth=`minw'"
end

* ---------------- Q: fastnose and e(V_es_att) ----------------
clear
set obs 90
gen long unit = _n
gen byte coh = 1 + mod(_n-1,3)
gen double ttre = .
replace ttre = 2004 if coh==1
replace ttre = 2007 if coh==2
gen double ufe = rnormal()
expand 10
bysort unit: gen int t = 1999 + _n
gen double tfe = 0.10*(t-1999)
gen byte trt = (ttre<.) & (t>=ttre)
gen double rel = t - ttre
gen double y = 2 + ufe + tfe + cond(trt==1, 1 + 0.30*rel, 0) + rnormal()*0.5
compress
save "r03_dat.dta", replace

dis "===== Q1: fastnose, esfixedbaseperiod ====="
use "r03_dat.dta", clear
cap noi wooldid y coh t ttre, att esprelength(3) espostlength(3) esfixedbaseperiod fastnose
dis "Q1_RC=" _rc
cap matrix drop zzb
cap matrix zzb = e(b_es_att)
dis "Q1|b_es_att_rc=" _rc
cap matrix drop zzv
cap matrix zzv = e(V_es_att)
dis "Q1|V_es_att_rc=" _rc
cap noi dis "Q1|V_es_att_rows=" rowsof(zzv) "|cols=" colsof(zzv)
cap noi matrix list zzv
dis "Q1|numpre_att=" e(numpre_att)
dis "Q1|wdid_vcvstate=" "`e(wdid_vcvstate)'"
dis "Q1|es_ndegen_att=" e(es_ndegen_att)
dis "--- ereturn matrix inventory ---"
ereturn list

dis "===== Q2: analytic control, same spec ====="
use "r03_dat.dta", clear
cap noi wooldid y coh t ttre, att esprelength(3) espostlength(3) esfixedbaseperiod
dis "Q2_RC=" _rc
cap matrix drop zzv2
cap matrix zzv2 = e(V_es_att)
dis "Q2|V_es_att_rc=" _rc
cap noi dis "Q2|V_es_att_rows=" rowsof(zzv2)
dis "Q2|wdid_vcvstate=" "`e(wdid_vcvstate)'"

* ---------------- R: degenerate cell SE map ----------------
dis "===== R: gteffects cell SE degeneracy map ====="
foreach arm in "timetrends cluster(coh)" "cluster(coh)" "timetrends robust" "timetrends" "timetrends cluster(unit)" "robust" "cluster(unit)" {
  use "r03_dat.dta", clear
  local tg = subinstr("`arm'"," ","_",.)
  local tg = subinstr("`tg'","(","",.)
  local tg = subinstr("`tg'",")","",.)
  if "`tg'" == "" local tg iid
  cap noi wooldid y coh t ttre, att esprelength(2) espostlength(2) gteffects `arm'
  dis "R_`tg'_RC=" _rc
  sedegen, tag(R3coh_`tg')
}

* wider design: 6 onset cohorts + never-treated  -> 7 clusters at i
clear
set obs 210
gen long unit = _n
gen byte coh = 1 + mod(_n-1,7)
gen double ttre = .
forval c = 1/6 {
  replace ttre = 2001 + `c' if coh == `c'
}
gen double ufe = rnormal()
expand 12
bysort unit: gen int t = 1999 + _n
gen double tfe = 0.10*(t-1999)
gen byte trt = (ttre<.) & (t>=ttre)
gen double rel = t - ttre
gen double y = 2 + ufe + tfe + cond(trt==1, 1 + 0.30*rel, 0) + rnormal()*0.5
compress
save "r03_dat7.dta", replace

foreach arm in "timetrends cluster(coh)" "cluster(coh)" "timetrends cluster(unit)" "timetrends robust" {
  use "r03_dat7.dta", clear
  local tg = subinstr("`arm'"," ","_",.)
  local tg = subinstr("`tg'","(","",.)
  local tg = subinstr("`tg'",")","",.)
  cap noi wooldid y coh t ttre, att esprelength(2) espostlength(2) gteffects `arm'
  dis "R7_`tg'_RC=" _rc
  sedegen, tag(R7coh_`tg')
}

dis "===== R03 DONE ====="
