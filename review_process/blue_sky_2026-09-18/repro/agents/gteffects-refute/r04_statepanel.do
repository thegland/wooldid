* r04: does the gteffects near-zero-SE degeneracy fire on the canonical applied spec --
*      a state panel with staggered adoption, i = state, cluster(state)?
*      Guard at ado:13496 is `se = (v < . & v > 0) ? sqrt(v) : .` -- EXACT zero only.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

cap prog drop sedegen
prog def sedegen
  syntax , tag(string) [show(integer 0)]
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
  local ji : list posof "i" in cn
  local jr : list posof "relyr" in cn
  mata: GG = st_matrix("`G'"); ///
        se = GG[.,`js']; pp = GG[.,`jp']; ll = GG[.,`jl']; uu = GG[.,`ju']; ///
        nse = select(se, se :< .); ///
        st_local("nr3", strofreal(rows(se))); ///
        st_local("nmiss", strofreal(rows(se) - rows(nse))); ///
        st_local("nzero", strofreal(sum(se :== 0))); ///
        st_local("ntiny", strofreal(sum((se :< 1e-8) :& (se :> 0)))); ///
        st_local("minse", strofreal(rows(nse)>0 ? min(nse) : ., "%12.4e")); ///
        st_local("minp",  strofreal(min(select(pp, pp :< .)), "%12.4e")); ///
        st_local("minw",  strofreal(min(select(uu-ll, (uu :< .) :& (ll :< .))), "%12.4e"))
  dis "SED|`tag'|rows=`nr3'|se_missing=`nmiss'|se_exact0=`nzero'|se_lt_1e8=`ntiny'|min_se=`minse'|min_p=`minp'|min_ciwidth=`minw'"
  if `show' > 0 {
    mata: GG = st_matrix("`G'"); bad = select(GG, (GG[.,`js'] :< 1e-8) :& (GG[.,`js'] :> 0)); k = min((rows(bad), `show')); for (r=1; r<=k; r++) printf("BADROW|`tag'|i=%g|relyr=%g|b=%12.8f|se=%12.4e|p=%12.4e|lb=%12.8f|ub=%12.8f\n", bad[r,`ji'], bad[r,`jr'], bad[r,`jb'], bad[r,`js'], bad[r,`jp'], bad[r,`jl'], bad[r,`ju'])
  }
end

* ---- 40-state panel, 20 years, three adoption waves, 10 never-treated states ----
clear
set obs 40
gen long state = _n
gen byte wave = mod(_n-1, 4)
gen double ttre = .
replace ttre = 2006 if wave==1
replace ttre = 2010 if wave==2
replace ttre = 2014 if wave==3
gen long cohid = wave + 1
gen double sfe = rnormal()*2
expand 20
bysort state: gen int year = 1999 + _n
gen double tfe = 0.05*(year-1999)
gen byte trt = (ttre<.) & (year >= ttre)
gen double rel = year - ttre
gen double y = 10 + sfe + tfe + cond(trt==1, 1.5 + 0.10*rel, 0) + rnormal()
compress
save "r04_states.dta", replace

dis "===== S1: i = state (unit level), cluster(state)  [canonical applied spec] ====="
use "r04_states.dta", clear
cap noi wooldid y state year ttre, att esprelength(3) espostlength(3) gteffects cluster(state)
dis "S1_RC=" _rc
sedegen, tag(S1_i_state_cl_state) show(6)
dis "S1|es_ndegen_att=" e(es_ndegen_att) "|vcvstate=`e(wdid_vcvstate)'"

dis "===== S2: i = state, robust ====="
use "r04_states.dta", clear
cap noi wooldid y state year ttre, att esprelength(3) espostlength(3) gteffects robust
dis "S2_RC=" _rc
sedegen, tag(S2_i_state_robust) show(6)

dis "===== S3: i = onset cohort, cluster(state)  [help's recommended shape] ====="
use "r04_states.dta", clear
cap noi wooldid y cohid year ttre, att esprelength(3) espostlength(3) gteffects cluster(state)
dis "S3_RC=" _rc
sedegen, tag(S3_i_coh_cl_state) show(6)

dis "===== S4: i = onset cohort, cluster(cohid)  [cluster AT i] ====="
use "r04_states.dta", clear
cap noi wooldid y cohid year ttre, att esprelength(3) espostlength(3) gteffects cluster(cohid)
dis "S4_RC=" _rc
sedegen, tag(S4_i_coh_cl_coh) show(8)
dis "S4|es_ndegen_att=" e(es_ndegen_att) "|vcvstate=`e(wdid_vcvstate)'"

dis "===== S5: i = state, cluster(region) with 8 regions nesting states ====="
use "r04_states.dta", clear
gen long region = 1 + mod(state-1, 8)
cap noi wooldid y state year ttre, att esprelength(3) espostlength(3) gteffects cluster(region)
dis "S5_RC=" _rc
sedegen, tag(S5_i_state_cl_region) show(6)

dis "===== S6: i = state, iid default ====="
use "r04_states.dta", clear
cap noi wooldid y state year ttre, att esprelength(3) espostlength(3) gteffects
dis "S6_RC=" _rc
sedegen, tag(S6_i_state_iid) show(6)

dis "===== R04 DONE ====="
