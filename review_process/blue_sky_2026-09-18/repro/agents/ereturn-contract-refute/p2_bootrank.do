* SKEPTIC probe of EC-1: is rank(e(V)) under bootstrap() a DEFECT or true structure?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 424242
set sortseed 424242
which wooldid

cap prog drop rkshow
prog def rkshow
  args mname tag
  cap confirm matrix `mname'
  if _rc != 0 {
    dis "`tag': MATRIX ABSENT"
    exit
  }
  mata: st_local("rk", strofreal(rank(st_matrix("`mname'"))))
  mata: st_local("kk", strofreal(cols(st_matrix("`mname'"))))
  mata: st_local("nn", strofreal(rows(st_matrix("`mname'"))))
  dis "`tag': rows=`nn' cols=`kk' rank=`rk'"
end

* ============ DESIGN A: deliberately IRREGULAR (unequal cohort sizes, ragged panel) ====
clear
set obs 130
gen long id = _n
gen byte clu = ceil(id/5)
gen double ui = rnormal()
gen int ttre = .
* unequal cohort sizes: 17 / 29 / 41 treated, rest never-treated
replace ttre = 2004 if id<=17
replace ttre = 2007 if id>17 & id<=46
replace ttre = 2010 if id>46 & id<=87
expand 13
bys id: gen int t = 2000 + _n - 1
gen double eps = rnormal()
* heterogeneous dynamic effects by cohort so the 4 estimands genuinely differ
gen double trteff = 0
replace trteff = (0.4 + 0.25*(ttre-2004)/3)*(1 + 0.35*(t-ttre)) if t>=ttre & ttre<.
gen double y = 1 + ui + 0.06*(t-2000) + trteff + eps
* ragged: drop a scattering of rows (administrative gaps)
set seed 9911
gen double u = runiform()
drop if u < 0.04
drop u
tab ttre, m

local SPEC "att att_i att_it att_itime esfixedbaseperiod esprelength(3) espostlength(3)"

foreach R in 5 8 12 20 40 60 {
  cap noisily wooldid y id t ttre, `SPEC' cluster(clu) bootstrap(standard) reps(`R') rseed(123)
  dis "A_RC_`R'=" _rc
  if _rc==0 {
    dis "A_R`R'  K=" colsof(e(b)) " completed=" e(boot_completed) " failed=" e(boot_failed)
    mat V = e(V)
    rkshow V "A_R`R'_eV"
    cap mat VB = e(V_boot)
    rkshow VB "A_R`R'_eVboot"
    cap mat BD = e(bootdist)
    cap rkshow BD "A_R`R'_bootdist_raw"
    * centered replication matrix rank
    cap mata: X = st_matrix("BD"); X = X :- mean(X); st_local("rkc", strofreal(rank(X)))
    dis "A_R`R'_bootdist_centered_rank=`rkc'"
  }
}

* which linear combinations are in the null space?  Look at the last fit.
dis "=== NAMES ==="
dis "`: colnames e(b)'"
mat b = e(b)
mat list b, format(%12.6f)
mat V = e(V)
mata: V = st_matrix("V"); symeigensystem(V, X=., L=.); L
dis "=== correlation of bootdist columns (|r|>0.99999 pairs) ==="
cap mat BD = e(bootdist)
mata:
  Z = st_matrix("BD")
  nm = st_matrixcolstripe("BD")
  C = correlation(Z)
  for (a=1; a<=cols(C); a++) {
    for (bq=a+1; bq<=cols(C); bq++) {
      if (abs(C[a,bq]) > 0.99999) printf("NEARDUP %s ~ %s r=%12.10f\n", nm[a,2], nm[bq,2], C[a,bq])
    }
  }
end

* ============ DESIGN B: the REGULAR shape the prober used (equal cohorts, balanced) ====
clear
set obs 120
gen long id = _n
gen byte clu = ceil(id/6)
gen double ui = rnormal()
gen int ttre = .
replace ttre = 2004 if mod(id,5)==0
replace ttre = 2006 if mod(id,5)==1
replace ttre = 2008 if mod(id,5)==2
replace ttre = 2010 if mod(id,5)==3
expand 12
bys id: gen int t = 2000 + _n - 1
gen double y = 1 + ui + 0.05*(t-2000) + rnormal()
replace y = y + 2 + 0.3*(t-ttre) if t>=ttre & ttre<.

foreach R in 5 12 40 {
  cap noisily wooldid y id t ttre, `SPEC' cluster(clu) bootstrap(standard) reps(`R') rseed(77)
  dis "B_RC_`R'=" _rc
  if _rc==0 {
    dis "B_R`R'  K=" colsof(e(b))
    mat V = e(V)
    rkshow V "B_R`R'_eV"
    cap mat BD = e(bootdist)
    cap mata: X = st_matrix("BD"); X = X :- mean(X); st_local("rkc", strofreal(rank(X)))
    dis "B_R`R'_bootdist_centered_rank=`rkc'"
  }
}
dis "=== B NAMES/points ==="
mat b = e(b)
mat list b, format(%14.9f)

* ============ C: is the low-rep F inflation a BIAS or just noise?  5 seeds. ====
foreach S in 11 22 33 44 55 {
  cap noisily wooldid y id t ttre, `SPEC' cluster(clu) bootstrap(standard) reps(5) rseed(`S')
  if _rc==0 {
    cap test ES_att_post0 = 0
    local f5 = r(F)
    local d5 = r(df)
    cap test ES_att_post0 ES_att_post1 ES_att_post2
    dis "SEED`S' reps5  F(post0)=" `f5' " ; joint F=" r(F) " df=" r(df) " p=" r(p)
  }
}
cap noisily wooldid y id t ttre, `SPEC' cluster(clu) bootstrap(standard) reps(200) rseed(99)
if _rc==0 {
  cap test ES_att_post0 ES_att_post1 ES_att_post2
  dis "REPS200 joint F=" r(F) " df=" r(df) " p=" r(p)
  mat V = e(V)
  rkshow V "R200_eV"
}
dis "ALLDONE_P2"
