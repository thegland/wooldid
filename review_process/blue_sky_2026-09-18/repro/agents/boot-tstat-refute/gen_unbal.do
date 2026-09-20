* gen_unbal.do -- MY OWN generator: deliberately UNBALANCED staggered panel.
*   varying #units per state, ragged panel (random row deletions), irregular cohort
*   sizes, and a subgroup indicator that is UNEVENLY distributed across cohorts.
* knobs: uG (states, 30), uT (periods, 8), uSkew (0/1), uRagged (0/1), uSGmode
*   uSGmode = "bal"   -> sg = mod(state,2)            (balanced across cohorts)
*   uSGmode = "unbal" -> sg = 1 for a lopsided, cohort-uneven set of states
capture drop _all
if "`uG'"      == "" local uG 30
if "`uT'"      == "" local uT 8
if "`uSkew'"   == "" local uSkew 0
if "`uRagged'" == "" local uRagged 1
if "`uSGmode'" == "" local uSGmode unbal
clear
set obs `uG'
gen int state = _n
gen int onset = .
* irregular cohort sizes: about 70% treated in three uneven cohorts, rest never
local uTR = round(0.70*`uG')
local uC1 = round(0.24*`uG')
local uC2 = round(0.48*`uG')
qui replace onset = 2002 if state <= `uC1'
qui replace onset = 2003 if state >  `uC1' & state <= `uC2'
qui replace onset = 2004 if state >  `uC2' & state <= `uTR'
gen double a_i = rnormal()*1.0
* VARYING number of units per state: 2,3,4,5 cycling
gen int nun = 2 + mod(state, 4)
expand nun
bys state: gen int unit = _n + 100*state
gen double u_i = rnormal()*0.6
expand `uT'
bys unit: gen int t = 1999 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = 1.0 + 0.2*(t - onset)
qui replace tau = 0 if !post
if `uSkew' gen double eps = (exp(rnormal()*1.2) - exp(1.2^2/2))*0.7
else       gen double eps = rnormal()*0.7
gen double y = a_i + u_i + tshock + tau*post + eps
gen int ttre = onset
* subgroup
gen byte sg = 0
if "`uSGmode'" == "bal"   qui replace sg = mod(state, 2)
if "`uSGmode'" == "unbal" qui replace sg = inlist(state,1,2,3,6,7,14,22,23,24) | (state > `uTR' & mod(state,3)==0)
* ragged panel: delete ~8% of rows at random (administrative-style gaps)
if `uRagged' {
  gen double _rr = runiform()
  qui drop if _rr < 0.08
  drop _rr
}
