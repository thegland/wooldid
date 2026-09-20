*=====================================================================
* LENS 2 (realism & impact) adjudication of FCM-4.
* Q1 trigger breadth : does a missing merged CONTROL (not missing y)
*                      pad the resampling-cluster count?  And a cluster
*                      with NO missing values at all (singleton units)?
* Q2 opt-in          : is e(df_r) padded on the DEFAULT invocation
*                      (cluster() + bootstrap() + reps(), nothing else)?
* Q3 reach           : does the padded df reach lincom / test?
* Q4 notice          : are the few-cluster cautions silenced exactly
*                      where the df error is largest?
* Q5 guard           : sthlp says bootstrap() "requires at least two
*                      clusters IN THE ESTIMATION SAMPLE".  Does a run
*                      with ONE live cluster and 23 dead ones refuse?
*=====================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 20260917
set sortseed 20260917

prog def mkpanel
  syntax , [seed(int 7001) nclus(int 24) nfirm(int 5) nyr(int 10)]
  clear
  set seed `seed'
  set sortseed `seed'
  set obs `=`nclus'*`nfirm''
  gen long firm = _n
  gen long ind  = ceil(_n/`nfirm')
  gen int onset = .
  replace onset = 2004 if mod(firm,4)==1
  replace onset = 2006 if mod(firm,4)==2
  replace onset = 2008 if mod(firm,4)==3
  bysort ind (firm): gen double ishk = rnormal() if _n==1
  by ind: replace ishk = ishk[1]
  gen double ai = rnormal()
  expand `nyr'
  bysort firm (ind): gen int year = 1999 + _n
  gen double x = 0.5*ai + rnormal()
  gen double y = ai + ishk + 0.04*(year-2000) + 1.2*((year>=onset) & !missing(onset)) + 0.3*x + rnormal()*0.4
  sort firm year
end

prog def rep
  syntax , tag(string)
  dis "RES|`tag'|N=" e(N) "|N_clust=" e(N_clust) "|df_r=" e(df_r) "|boot_df=" e(boot_df) "|df_r_analytic=" e(df_r_analytic)
  tempname b V
  mat `b' = e(b)
  mat `V' = e(V)
  local nm : colnames `b'
  forvalues q=1/`=colsof(`b')' {
     local this : word `q' of `nm'
     if "`this'"=="main_att_i" {
        local se = sqrt(`V'[`q',`q'])
        dis "RES|`tag'|att_i b=" %14.10f `b'[1,`q'] " se=" %14.10f `se'
     }
  }
end

*---------------------------------------------------------------------
* Q1a  TRIGGER = merged control with partial coverage.  12 of 24
*      industries have no x at all.  Nothing is wrong with y.
*---------------------------------------------------------------------
dis "===== A1 padded: controls(x), x missing in 12 of 24 industries ====="
mkpanel
qui replace x = . if ind <= 12
cap noi wooldid y firm year onset , cluster(ind) att att_i controls(x) bootstrap(standard) reps(40) rseed(11) bootdf(clusters) bootdisplay(wald) dots(0)
dis "A1_RC=" _rc
rep , tag(A1_PAD_CTRL)

dis "===== A2 control: same estimation sample, the 12 dead industries dropped ====="
mkpanel
qui drop if ind <= 12
cap noi wooldid y firm year onset , cluster(ind) att att_i controls(x) bootstrap(standard) reps(40) rseed(11) bootdf(clusters) bootdisplay(wald) dots(0)
dis "A2_RC=" _rc
rep , tag(A2_DROP_CTRL)

*---------------------------------------------------------------------
* Q2/Q3  DEFAULT invocation - no bootdf(), no bootdisplay().
*---------------------------------------------------------------------
dis "===== B padded, DEFAULT invocation (no bootdf, no bootdisplay) ====="
mkpanel
qui replace x = . if ind <= 12
cap noi wooldid y firm year onset , cluster(ind) att att_i controls(x) bootstrap(standard) reps(40) rseed(11) dots(0)
dis "B_RC=" _rc
rep , tag(B_PAD_DEFAULT)
dis "----- Q3: what df do lincom and test use? -----"
cap noi lincom _b[main_att_i]
dis "LINCOM_RC=" _rc
dis "RES|B|lincom df from e(df_r)=" e(df_r) " ; visible cluster count e(N_clust)=" e(N_clust)
cap noi test _b[main_att_i] = 0
dis "TEST_RC=" _rc
dis "RES|B|r(df_r)=" r(df_r) " r(F)=" r(F) " r(p)=" r(p)
dis "RES|B|p at padded df=" %12.10f 2*ttail(e(df_r), abs(_b[main_att_i]/sqrt(e(V)[1,1]))) ///
    " vs at live df=" %12.10f 2*ttail(e(N_clust)-1, abs(_b[main_att_i]/sqrt(e(V)[1,1])))

dis "===== B2 control for B: dead dropped, default ====="
mkpanel
qui drop if ind <= 12
cap noi wooldid y firm year onset , cluster(ind) att att_i controls(x) bootstrap(standard) reps(40) rseed(11) dots(0)
dis "B2_RC=" _rc
rep , tag(B2_DROP_DEFAULT)

*---------------------------------------------------------------------
* Q4  small G_live: 4 live industries of 24.  Are the cautions silenced?
*---------------------------------------------------------------------
dis "===== D1 padded, only 4 live industries of 24 ====="
mkpanel
qui replace x = . if ind > 4
cap noi wooldid y firm year onset , cluster(ind) att att_i controls(x) bootstrap(standard) reps(40) rseed(11) dots(0)
dis "D1_RC=" _rc
rep , tag(D1_PAD_G4)

dis "===== D2 control, the 20 dead industries dropped ====="
mkpanel
qui keep if ind <= 4
cap noi wooldid y firm year onset , cluster(ind) att att_i controls(x) bootstrap(standard) reps(40) rseed(11) dots(0)
dis "D2_RC=" _rc
rep , tag(D2_DROP_G4)

*---------------------------------------------------------------------
* Q5  ONE live industry, 23 dead.  sthlp: bootstrap() "requires at
*     least two clusters in the estimation sample".
*---------------------------------------------------------------------
dis "===== E one live industry of 24 - does the >=2 refusal fire? ====="
mkpanel
qui replace x = . if ind > 1
cap noi wooldid y firm year onset , cluster(ind) att att_i controls(x) bootstrap(standard) reps(40) rseed(11) dots(0)
dis "E_RC=" _rc
if _rc == 0 rep , tag(E_PAD_G1)

dis "===== E2 control: the single live industry alone ====="
mkpanel
qui keep if ind == 1
cap noi wooldid y firm year onset , cluster(ind) att att_i controls(x) bootstrap(standard) reps(40) rseed(11) dots(0)
dis "E2_RC=" _rc

*---------------------------------------------------------------------
* Q1b  a dead cluster with NO missing values anywhere: every firm of
*      industry 1 observed in a single year (short-lived entrants), so
*      reghdfe drops them as singletons.
*---------------------------------------------------------------------
dis "===== F dead cluster by singleton drop, no missing values at all ====="
mkpanel
qui drop if ind == 1 & year != 2005
qui count if missing(y) | missing(x)
dis "F|rows with any missing = " r(N)
cap noi wooldid y firm year onset , cluster(ind) att att_i controls(x) bootstrap(standard) reps(40) rseed(11) dots(0)
dis "F_RC=" _rc
rep , tag(F_SINGLETON)

dis "L2_FCM4_DONE"
