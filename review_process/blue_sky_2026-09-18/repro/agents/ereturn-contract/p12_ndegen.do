* ==============================================================================
* SELF-CONTAINED REPRO: e(V_es_att) singular with a strictly positive diagonal,
* e(es_ndegen_att) ABSENT, on a canonical 3-cohort staggered panel.
* sthlp (the e(b_es)/e(V_es) entry): "...flagged by a Caution and by
* e(es_ndegen_<type>) (the count of degenerate rows, posted only when such rows
* exist - ITS ABSENCE MARKS A CLEAN COVARIANCE)".
* ado:4130-4136 computes that count by scanning ONLY the diagonal, so a block that
* is singular off the diagonal is never flagged.
* ==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/ereturn-contract"
which wooldid

* --- a completely ordinary staggered panel: 120 units, 12 periods, 3 cohorts ---
clear
set obs 120
gen long id = _n
gen double ttre = .
replace ttre = 4  if inrange(id,1,30)
replace ttre = 7  if inrange(id,31,60)
replace ttre = 10 if inrange(id,61,80)
gen long clu = ceil(id/6)
gen double ufe = rnormal()
expand 12
bys id: gen int t = _n
gen byte post = (t>=ttre) & ttre<.
gen double y = ufe + .2*t + 1.0*post + rnormal()

cap prog drop SHOW
prog def SHOW
  args tag fam
  cap confirm matrix e(V_es_`fam')
  if _rc {
    dis "SHOW|`tag'|`fam'|ABSENT"
    exit
  }
  tempname V B
  matrix `V' = e(V_es_`fam')
  matrix `B' = e(b_es_`fam')
  local H = colsof(`V')
  local rk = .
  local mn = .
  local md = .
  local ev = ""
  cap mata: st_local("rk", strofreal(rank(st_matrix("e(V_es_`fam')"))))
  cap mata: st_local("mn", strofreal(min(Re(eigenvalues(st_matrix("e(V_es_`fam')")))), "%12.4e"))
  cap mata: st_local("md", strofreal(min(diagonal(st_matrix("e(V_es_`fam')"))), "%12.4e"))
  cap mata: st_local("ev", invtokens(strofreal(sort(Re(eigenvalues(st_matrix("e(V_es_`fam')")))', 1)', "%9.3e")))
  local nd = .
  cap local nd = e(es_ndegen_`fam')
  dis "NDG|`tag'|`fam'|H=`H'|rank=`rk'|mineig=`mn'|mindiag=`md'|es_ndegen=`nd'|numpre=" e(numpre_`fam')
  dis "EIG|`tag'|`fam'|`ev'"
  * the matrix consumers the help points at
  cap mata: st_local("ch", strofreal(sum(cholesky(makesymmetric(st_matrix("e(V_es_`fam')"))) :>= .)))
  dis "CHOL|`tag'|`fam'|missingcells=`ch'|rc=" _rc
  matrix WDB = e(b_es_`fam')
  matrix WDV = e(V_es_`fam')
  cap noi honestdid, pre(1/`=e(numpre_`fam')') post(`=e(numpre_`fam')+1'/`H') b(WDB) vcv(WDV) parallel(0)
  dis "HONESTDID|`tag'|`fam'|rc=" _rc
  cap noi pretrends power 0.5, pre(1/`=e(numpre_`fam')') post(`=e(numpre_`fam')+1'/`H') b(WDB) v(WDV)
  dis "PRETRENDS|`tag'|`fam'|rc=" _rc
end

dis "########## 1. clustered ##########"
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) cluster(clu)
dis "C1_RC=" _rc
SHOW clustered att

dis "########## 2. plain (iid) ##########"
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5)
dis "C2_RC=" _rc
SHOW iid att

dis "########## 3. robust ##########"
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) robust
dis "C3_RC=" _rc
SHOW robust att

dis "########## 4. att only ##########"
cap noi wooldid y id t ttre, att esfixedbaseperiod esprelength(3) espostlength(5) cluster(clu)
dis "C4_RC=" _rc
SHOW attonly att

dis "########## 5. a window the design supports (control) ##########"
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(2) espostlength(2) cluster(clu)
dis "C5_RC=" _rc
SHOW narrow att

dis "########## 6. pretrendpower(), which consumes the same matrix ##########"
cap noi wooldid y id t ttre, att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(5) cluster(clu) pretrendpower(0.5)
dis "C6_RC=" _rc " roth_ok=" e(roth_ok) " mdt=" e(roth_mdt_att) " slope=" e(roth_slope_att)
SHOW ptp att

dis "NDEGEN_DONE"
