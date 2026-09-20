* gen_small.do -- small staggered panel, i = state (cluster == i), G states.
* knobs (set before include): gG (states, default 20), gUnits (per state, 3),
*   gT0 2000, gT 6, gonsets "2002 2003 2004", gnever (# never-treated states, 5),
*   gysd .7, gtauslope .2, gtau0 1.0, gskew 0 (1 = lognormal error -> skewed t)
capture drop _all
if "`gG'"      == "" local gG 20
if "`gUnits'"  == "" local gUnits 3
if "`gT0'"     == "" local gT0 2000
if "`gT'"      == "" local gT 6
if "`gonsets'" == "" local gonsets "2002 2003 2004"
if "`gnever'"  == "" local gnever 5
if "`gysd'"    == "" local gysd 0.7
if "`gtau0'"   == "" local gtau0 1.0
if "`gtauslope'" == "" local gtauslope 0.2
if "`gskew'"   == "" local gskew 0
local nons : word count `gonsets'
clear
set obs `gG'
gen int state = _n
gen int onset = .
local tr = `gG' - `gnever'
forval s = 1/`tr' {
  local w = mod(`s'-1, `nons') + 1
  local o : word `w' of `gonsets'
  qui replace onset = `o' in `s'
}
gen double a_i = rnormal()*1.0
expand `gUnits'
bys state: gen int unit = _n + 1000*state
gen double u_i = rnormal()*0.5
expand `gT'
bys unit: gen int t = `gT0' - 1 + _n
gen double tshock = 0
qui bys t (unit): replace tshock = rnormal() if _n==1
qui bys t (unit): replace tshock = tshock[1]
gen byte post = (t >= onset) & !missing(onset)
gen double tau = `gtau0' + `gtauslope'*(t - onset)
qui replace tau = 0 if !post
if `gskew' gen double eps = (exp(rnormal()*1.2) - exp(1.2^2/2))*`gysd'
else       gen double eps = rnormal()*`gysd'
gen double y = a_i + u_i + tshock + tau*post + eps
gen int ttre = onset
