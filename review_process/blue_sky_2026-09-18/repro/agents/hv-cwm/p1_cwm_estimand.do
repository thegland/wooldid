* ============================================================================
* p1_cwm_estimand.do
* Q1: does customweightmultiplier() move att_it / att_i / att_itime off the
*     documented "simple, unweighted mean of the cell effects" estimand?
* Q2: does the SAME multiplier applied as a sample weight (which is what
*     bootstrap(bayes) does, ado:9449-9451: repw = basew*drawv) cancel out,
*     as the bootctl cell-weighted gate comment (ado:6430-6439) asserts?
* Self-contained: generates its own panel; no project test data used.
* ============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916

which wooldid

mata:
void wd_chk(string scalar arm)
{
    real matrix G, M
    string colvector gcn, mcn
    string matrix mrs
    real colvector relyr, est, post, wcol, ii
    real scalar ci, k, wm, sm, bb
    string rowvector ests
    G   = st_matrix("G")
    gcn = st_matrixcolstripe("G")[.,2]
    M   = st_matrix("M")
    mrs = st_matrixrowstripe("M")
    mcn = st_matrixcolstripe("M")[.,2]
    ci  = selectindex(mcn:=="estimate")[1]
    relyr = G[., selectindex(gcn:=="relyr")[1]]
    est   = G[., selectindex(gcn:=="estimate")[1]]
    post  = selectindex(relyr:>=0)
    ests  = ("att","att_it","att_i","att_itime")
    for (k=1; k<=4; k++) {
        wcol = G[., selectindex(gcn:=="w_"+ests[k])[1]]
        wm   = sum(est[post]:*wcol[post])/sum(wcol[post])
        sm   = mean(est[post])
        ii   = selectindex((mrs[.,1]:=="main"):*(mrs[.,2]:==ests[k]))
        if (length(ii)==0) bb = .
        else bb = M[ii[1], ci]
        printf("CHK %s %s rep=%12.8f wtd=%12.8f simple=%12.8f rdW=%g rdS=%g\n", arm, ests[k], bb, wm, sm, abs(wm-bb)/abs(bb), abs(sm-bb)/abs(bb))
        printf("HEX %s %s rep=%s wtd=%s simple=%s\n", arm, ests[k], strofreal(bb,"%21x"), strofreal(wm,"%21x"), strofreal(sm,"%21x"))
    }
}
end

* ---------------------------------------------------------------- data ----
* 8 units (4 treated, staggered; 4 never treated), 10 periods,
* 5 individuals per unit-period cell.  Sample weights vary WITHIN cell.
clear
set obs 8
gen int i = _n
gen int ttre = .
replace ttre = 4 if i==1
replace ttre = 5 if i==2
replace ttre = 6 if i==3
replace ttre = 7 if i==4
gen double pop = cond(inlist(i,3,4), 50, 1)   // cell-constant "population"
expand 10
bys i: gen int t = _n
expand 5
bys i t: gen byte k = _n
gen long obsid = _n
gen double w      = 0.5 + runiform()
gen double cwmind = 0.2 + 2*runiform()        // cwm varying WITHIN cell
gen double mfine  = 0.5 + 1.5*runiform()      // weight mult varying WITHIN cell
gen byte D  = (ttre < .) & (t >= ttre)
gen double tau = cond(i==1,1.0,cond(i==2,2.0,cond(i==3,3.0,cond(i==4,4.0,0))))
replace tau = tau + 0.2*(t-ttre) if D==1
replace tau = 0 if D==0
gen double y = 0.3*i + 0.1*t + tau*D + 0.05*rnormal()
gen double wm  = w*pop
gen double wm2 = w*mfine
compress
save "panel.dta", replace

local A_w  w
local A_o
local A_c  1
local B_w  w
local B_o  customweightmultiplier(pop)
local B_c  pop
local C_w  wm
local C_o
local C_c  1
local D_w  w
local D_o  customweightmultiplier(cwmind)
local D_c  cwmind
local E_w  wm2
local E_o
local E_c  1

foreach arm in A B C D E {
  dis _n(2) "======================== ARM `arm' ========================"
  use "panel.dta", clear
  local WV  "``arm'_w'"
  local OPT "``arm'_o'"
  local CWX "``arm'_c'"
  dis as text "weight = `WV' ; option = `OPT' ; oracle cwm = `CWX'"
  cap noi wooldid y i t ttre [aw=`WV'], att att_it att_i att_itime gteffects `OPT'
  local rc = _rc
  dis "ARM`arm'_RC=`rc'"
  if `rc' == 0 {
    dis "ARM`arm'_N=" e(N)
    mat G = e(gteffects)
    mat M = e(wdidmainresults)
    if "`arm'" == "A" mat list M
    preserve
      gen double _w  = `WV'
      gen double _wc = `WV' * (`CWX')
      collapse (sum) _w _wc, by(i t)
      gen double cellwt = _wc/_w
      tempfile ORC
      save `ORC'
      clear
      svmat double G, names(col)
      rename caltime t
      merge m:1 i t using `ORC', keep(match master) nogen
      gen double reld_it = abs(w_att_it - cellwt)/max(abs(cellwt),1e-30)
      dis as result "---- ARM `arm': w_att_it vs sum(w*cwm)/sum(w) oracle ----"
      qui su reld_it
      dis "ARM`arm'_MAXRELD=" r(max)
      qui su w_att_it
      dis "ARM`arm'_WIT min=" r(min) " max=" r(max) " sd=" r(sd)
      sort i t
      list i t relyr estimate w_att w_att_it w_att_i w_att_itime cellwt if relyr>=0, ///
           sepby(i) noobs
      dis as result "---- ARM `arm': pre rows, i=1 and i=3 ----"
      list i t relyr w_att w_att_it w_att_i w_att_itime if relyr<0 & inlist(i,1,3), ///
           sepby(i) noobs
    restore
    mata: wd_chk("`arm'")
  }
}
dis "DONE_P1"
