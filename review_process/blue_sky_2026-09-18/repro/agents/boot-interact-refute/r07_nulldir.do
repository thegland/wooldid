* r07_nulldir.do -- SELF-CONTAINED.  DECIDES the BOOTX-2 "the bootstrap hides a real
* rank deficiency" question.  Take v = the eigenvector of the ANALYTIC pre-period block
* whose eigenvalue is ~0.  If the zero direction is a REAL property of the estimator,
* v'b_r must be (near) constant across bootstrap replications, and the bootstrap
* covariance's full numerical rank is spurious.  If v'S_boot v is an ordinary-sized
* variance, the analytic zero is an artifact of the sandwich VCV and the bootstrap's
* full rank is CORRECT.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
clear
set obs 30
gen state = _n
gen onset = .
replace onset = 7 + mod(_n, 6) if mod(_n,3)==0
expand 2
bysort state: gen cty = _n
gen county = state*10 + cty
expand 16
bysort county: gen t = _n
bysort county (t): gen double ushock = rnormal(0,.3) if _n==1
bysort county (t): replace ushock = ushock[_n-1] if _n>1
gen double e = rnormal(0,1)
gen double rel = t - onset
gen double te = cond(onset<. & rel>=0, .8 + .12*rel, 0)
gen double y = 2 + .05*t + ushock + te + e
drop rel te

local OPT "cluster(state) esfixedbaseperiod esprelength(6) espostlength(4) pretrendpower(80)"

qui wooldid y county t onset, `OPT'
matrix Va = e(V_es_att)
matrix ba = e(b_es_att)
local ln : colnames ba
local k 0
local idx
local j 0
foreach l of local ln {
  local ++j
  if substr("`l'",1,3)=="pre" {
    local ++k
    local idx `idx' `j'
  }
}
matrix Sa = J(`k',`k',0)
local r 0
foreach i1 of local idx {
  local ++r
  local c 0
  foreach i2 of local idx {
    local ++c
    matrix Sa[`r',`c'] = Va[`i1',`i2']
  }
}
matrix Sa = (Sa+Sa')/2
mata:
  S = st_matrix("Sa")
  X = .; L = .
  symeigensystem(S, X, L)
  // last column of X is the smallest eigenvalue's eigenvector
  ord = order(L', 1)
  v = X[., ord[1]]
  st_matrix("vnull", v)
  printf("ANALYTIC smallest eigenvalue = %12.4e\n", L[ord[1]])
  printf("ANALYTIC v'Sa v            = %12.4e\n", v'*S*v)
  printf("ANALYTIC trace(Sa)/k       = %12.4e\n", trace(S)/rows(S))
end
matrix list vnull

qui wooldid y county t onset, `OPT' bootstrap(standard) reps(200) rseed(777)
matrix Vb = e(V_es_att)
matrix Sb = J(`k',`k',0)
local r 0
foreach i1 of local idx {
  local ++r
  local c 0
  foreach i2 of local idx {
    local ++c
    matrix Sb[`r',`c'] = Vb[`i1',`i2']
  }
}
matrix Sb = (Sb+Sb')/2
mata:
  v = st_matrix("vnull")
  Sb = st_matrix("Sb")
  printf("BOOT   v'Sb v              = %12.4e\n", v'*Sb*v)
  printf("BOOT   trace(Sb)/k         = %12.4e\n", trace(Sb)/rows(Sb))
  printf("BOOT   ratio v'Sbv / (tr/k)= %12.4e\n", (v'*Sb*v)/(trace(Sb)/rows(Sb)))
end

* and directly: the variance of v'b_r across the replication draws
cap confirm matrix e(bootdist)
dis "BOOTDIST_RC=" _rc
if _rc==0 {
  matrix BD = e(bootdist)
  local bdn : colnames BD
  dis "BOOTDIST rows=" rowsof(BD) " cols=" colsof(BD)
  mata:
    BD  = st_matrix("BD")
    nm  = st_matrixcolstripe("BD")[.,2]
    lab = ("pre7","pre6","pre5","pre4","pre3","pre2")
    id  = J(1, cols(lab), .)
    for (a=1; a<=cols(lab); a++) {
      for (b=1; b<=rows(nm); b++) {
        if (nm[b] == "ES_att_" + lab[a]) id[a] = b
      }
    }
    id
    C = BD[., id']
    ok = rowmissing(C) :== 0
    C = select(C, ok)
    printf("listwise complete replications for the 6 pre columns = %f\n", rows(C))
    v = st_matrix("vnull")
    z = C * v
    printf("SD of v'b_r across replications      = %12.4e\n", sqrt(variance(z)))
    // compare with the SD of the largest-eigenvalue direction
    S = st_matrix("Sa"); X=.; L=.
    symeigensystem(S, X, L)
    ord = order(L', 1)
    w = X[., ord[cols(ord)]]
    printf("SD of w'b_r (largest analytic dir)   = %12.4e\n", sqrt(variance(C*w)))
    printf("analytic sqrt(v'Sa v)                = %12.4e\n", sqrt(v'*S*v))
    printf("analytic sqrt(w'Sa w)                = %12.4e\n", sqrt(w'*S*w))
  end
}
dis "DONE"
* NOTE (author of this probe): the two lines printing "w'b_r (largest analytic dir)" and
* "analytic sqrt(w'Sa w)" are WRONG and should be ignored - `ord' is a colvector so
* cols(ord) is 1 and w ends up equal to v; and sqrt() of a 1x1 matrix expression returns
* missing.  The load-bearing numbers are v'Sa v, v'Sb v, trace(Sb)/k and the SD of v'b_r,
* all of which are correct.
