* ==============================================================================
* SINGULAR e(V_es_*) WITH A STRICTLY POSITIVE DIAGONAL
* sthlp: "...flagged by a Caution and by e(es_ndegen_<type>) (the count of
* degenerate rows, posted only when such rows exist - ITS ABSENCE MARKS A CLEAN
* COVARIANCE)".  e(es_ndegen_*) counts only ZERO/MISSING DIAGONAL entries
* (ado:4130-4136), so a block that is singular with every diagonal entry strictly
* positive is not flagged at all - and honestdid / pretrends / test invert it.
* ==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid

cap prog drop BLOCKRANK
prog def BLOCKRANK
  args tag
  local sgl `e(wdid_sglist)'
  local keys att att_it ame ame_it
  foreach v of local sgl {
    local keys `keys' att_s`v' att_it_s`v' ame_s`v' ame_it_s`v'
  }
  local rkv = .
  cap mata: st_local("rkv", strofreal(rank(st_matrix("e(V)"))))
  dis "VTOT|`tag'|K=" colsof(e(b)) "|rank(e(V))=`rkv'|vcvstate=`e(wdid_vcvstate)'|vblocks=`e(wdid_vblocks)'"
  foreach k of local keys {
    cap confirm matrix e(V_es_`k')
    if _rc continue
    tempname VB
    matrix `VB' = e(V_es_`k')
    local H = colsof(`VB')
    local rk = .
    local mn = .
    local mind = .
    cap mata: st_local("rk", strofreal(rank(st_matrix("e(V_es_`k')"))))
    cap mata: st_local("mn", strofreal(min(Re(eigenvalues(st_matrix("e(V_es_`k')")))), "%12.4e"))
    cap mata: st_local("mind", strofreal(min(diagonal(st_matrix("e(V_es_`k')"))), "%12.4e"))
    local nd = .
    cap local nd = e(es_ndegen_`k')
    local np = .
    cap local np = e(numpre_`k')
    local cn : colnames `VB'
    dis "BLK|`tag'|`k'|H=`H'|rank=`rk'|mineig=`mn'|mindiag=`mind'|es_ndegen=`nd'|numpre=`np'|cols=`cn'"
    matrix drop `VB'
  }
end

* ---------- geometry A: one subgroup adopts almost immediately ----------
clear
set obs 100
gen long id = _n
gen byte sg = mod(id,2)
gen double ttre = .
replace ttre = 2 if sg==0 & id<=50
replace ttre = 9 if sg==1 & id<=50
gen long clu = ceil(id/5)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .2*t + 1.1*post + rnormal()
save S1.dta, replace
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(6) espostlength(6) subgroups(sg) cluster(clu)
dis "A_RC=" _rc
BLOCKRANK A_shortpre6
* consequence 1: a joint test of every pre-period coefficient
cap noi test ES_att_SG1_pre7 ES_att_SG1_pre6 ES_att_SG1_pre5 ES_att_SG1_pre4 ES_att_SG1_pre3 ES_att_SG1_pre2
dis "A_TESTPRE requested=6 df=" r(df) " p=" r(p) " F=" r(F) " rc=" _rc
cap noi test ES_att_pre7 ES_att_pre6 ES_att_pre5 ES_att_pre4 ES_att_pre3 ES_att_pre2
dis "A_TESTPRE_POOLED requested=6 df=" r(df) " p=" r(p) " rc=" _rc
* consequence 2: the honestdid handoff the help advertises
cap noi honestdid, pre(1/6) post(7/12) b(e(b_es_att)) vcv(e(V_es_att)) parallel(0)
dis "A_HONESTDID_POOLED_RC=" _rc
cap noi honestdid, pre(1/6) post(7/10) b(e(b_es_att_s1)) vcv(e(V_es_att_s1)) parallel(0)
dis "A_HONESTDID_SG1_RC=" _rc
* consequence 3: invert it
cap mata: st_local("nm", strofreal(sum(invsym(st_matrix("e(V_es_att)")) :>= .)))
dis "A_INVSYM_MISS=`nm'"
cap mata: st_local("ch", strofreal(sum(cholesky(makesymmetric(st_matrix("e(V_es_att)"))) :>= .)))
dis "A_CHOL_MISS=`ch' rc=" _rc

* ---------- geometry B: a subgroup present in only one cohort ----------
clear
set obs 120
gen long id = _n
gen double ttre = .
replace ttre = 4 if inrange(id,1,30)
replace ttre = 7 if inrange(id,31,60)
replace ttre = 10 if inrange(id,61,80)
gen byte sg = 0
replace sg = 1 if inrange(id,1,30)
gen long clu = ceil(id/6)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .2*t + 1.0*post + rnormal()
save S2.dta, replace
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) subgroups(sg) cluster(clu)
dis "B_RC=" _rc
BLOCKRANK B_sgonecohort
cap noi test ES_att_SG1_pre3 ES_att_SG1_pre2
dis "B_TESTPRE requested=2 df=" r(df) " p=" r(p) " rc=" _rc

* ---------- geometry C: the same WITHOUT subgroups (is it the subgroup or the design?) ----------
use S2.dta, clear
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) cluster(clu)
dis "C_RC=" _rc
BLOCKRANK C_nosg
use S1.dta, clear
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(6) espostlength(6) cluster(clu)
dis "D_RC=" _rc
BLOCKRANK D_nosg_shortpre
* and a plain well-behaved panel, for the baseline
use D1.dta, clear
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) cluster(clu)
dis "E_RC=" _rc
BLOCKRANK E_baseline
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) subgroups(sg) cluster(clu)
dis "F_RC=" _rc
BLOCKRANK F_baseline_sg
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3) robust
dis "G_RC=" _rc
BLOCKRANK G_baseline_robust

dis "SINGULAR_DONE"
