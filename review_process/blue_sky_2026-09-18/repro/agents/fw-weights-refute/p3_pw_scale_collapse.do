* SELF-CONTAINED.  Three questions in one run.
*  (1) FW-02: under [pw=w] with no vce option, is e(V) the robust covariance, and what do
*      e(vce)/e(vcetype) say it is?
*  (2) NEW: weight-scale invariance -- ado:1655 asserts "weight scale here actually doesn't
*      matter". Test [aw=w] vs [aw=1000*w] (pre-generated variable, not an expression),
*      and the same under [pw].
*  (3) NEW: default cell collapse vs nofastcollapse when the weight VARIES WITHIN the
*      cohort-by-period cell (coarse i = state, many counties per state-year).
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
which wooldid

* --- 30 states x 12 years x 8 counties; i = state (coarse cohort) ---
set obs 30
gen long i = _n
gen ttre = .
replace ttre = 5 if mod(_n,4)==1
replace ttre = 8 if mod(_n,4)==2
replace ttre = 10 if mod(_n,4)==3
gen double ai = rnormal()
expand 12
bys i: gen int t = _n
gen double at = 0.2*t
expand 8
bys i t: gen byte c = _n
gen double aic = rnormal()*0.5
gen double y = ai + at + aic + rnormal() + 1.2*(t>=ttre & ttre<.)
* county population: varies a lot within state-year
gen double pop = exp(rnormal()*0.8 + 9)
gen double pop1000 = 1000*pop
compress

prog def grab
  args tag
  tempname B V
  matrix `B' = e(b)
  matrix `V' = e(V)
  foreach c in att att_it att_i att_itime {
    local cn = colnumb(`B', "main_`c'")
    if `cn' < . {
      dis "VAL_`tag'_`c' = " %21x (`B'[1,`cn']) "  SE " %21x (sqrt(`V'[`cn',`cn'])) "  dec " %20.15f (`B'[1,`cn']) " " %20.15f (sqrt(`V'[`cn',`cn']))
    }
  }
  dis "EVCE_`tag'=[" e(vce) "] EVCETYPE_`tag'=[" e(vcetype) "] N_`tag'=" e(N)
end

*******************  (1) FW-02  *******************
dis _n "=== PW1: [pw=pop], no vce option ==="
capture noisily wooldid y i t ttre [pw=pop], att
dis "PW1_RC=" _rc
grab PW1
tempname Vpw
matrix `Vpw' = e(V)

dis _n "=== PW2: [aw=pop], robust ==="
capture noisily wooldid y i t ttre [aw=pop], att robust
dis "PW2_RC=" _rc
grab PW2
tempname Vawr
matrix `Vawr' = e(V)

dis _n "=== PW3: [aw=pop], iid default ==="
capture noisily wooldid y i t ttre [aw=pop], att
dis "PW3_RC=" _rc
grab PW3
tempname Vawi
matrix `Vawi' = e(V)

dis "VDIFF_pw_vs_awrobust = " %21x (mreldif(`Vpw', `Vawr'))
dis "VDIFF_pw_vs_awiid    = " %21x (mreldif(`Vpw', `Vawi'))

dis _n "=== reghdfe direct with [pw=pop] ==="
gen byte post = (t>=ttre & ttre<.)
capture noisily reghdfe y post [pw=pop], absorb(i t)
dis "RH_RC=" _rc
dis "RH_EVCE=[" e(vce) "] RH_EVCETYPE=[" e(vcetype) "]"
capture noisily regress y post i.i i.t [pw=pop]
dis "REG_EVCE=[" e(vce) "] REG_EVCETYPE=[" e(vcetype) "]"

dis _n "=== does esttab/estout read e(vcetype)? ==="
qui wooldid y i t ttre [pw=pop], att
capture noisily est store wpw
capture noisily esttab wpw, se nostar
dis "ESTTAB_RC=" _rc

*******************  (2) scale invariance  *******************
dis _n "=== SC1: [aw=pop]  4 estimands, robust ==="
capture noisily wooldid y i t ttre [aw=pop], att att_it att_i att_itime robust
dis "SC1_RC=" _rc
grab SC1
dis _n "=== SC2: [aw=pop1000] (=1000*pop) ==="
capture noisily wooldid y i t ttre [aw=pop1000], att att_it att_i att_itime robust
dis "SC2_RC=" _rc
grab SC2
dis _n "=== SC3: [pw=pop] ==="
capture noisily wooldid y i t ttre [pw=pop], att att_it att_i att_itime
dis "SC3_RC=" _rc
grab SC3
dis _n "=== SC4: [pw=pop1000] ==="
capture noisily wooldid y i t ttre [pw=pop1000], att att_it att_i att_itime
dis "SC4_RC=" _rc
grab SC4

*******************  (3) collapse vs nofastcollapse  *******************
foreach v in "" "robust" "cluster(i)" {
  local tg = cond("`v'"=="", "iid", cond("`v'"=="robust","rob","clu"))
  dis _n "=== FC_`tg': default collapse, `v' ==="
  capture noisily wooldid y i t ttre [aw=pop], att att_it att_i att_itime `v'
  dis "FCON_`tg'_RC=" _rc
  grab FCON_`tg'
  dis _n "=== FC_`tg': nofastcollapse, `v' ==="
  capture noisily wooldid y i t ttre [aw=pop], att att_it att_i att_itime `v' nofastcollapse
  dis "FCOFF_`tg'_RC=" _rc
  grab FCOFF_`tg'
}
dis "P3_DONE"
