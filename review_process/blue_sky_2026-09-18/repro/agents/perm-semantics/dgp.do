*! dgp.do -- self-contained staggered DiD generator for the perm-semantics probes.
*! Knobs (set as locals BEFORE `include'): dN (units), dT (periods), dons (onset list),
*! dnever (never-treated count), dtau (true effect), dsd (noise sd), dseed.
if "`dN'"    == "" local dN    24
if "`dT'"    == "" local dT    8
if `"`dons'"' == "" local dons 2004 2006
if "`dnever'" == "" local dnever 8
if "`dtau'"  == "" local dtau  1.0
if "`dsd'"   == "" local dsd   1.0
if "`dseed'" == "" local dseed 424242
clear
set seed `dseed'
set sortseed `dseed'
local K : word count `dons'
local ntr = `dN' - `dnever'
set obs `dN'
gen long unit = _n
gen double u = runiform()
sort u
gen long grp = .
* assign the first `ntr' units round-robin to the K onsets, the rest never-treated
forvalues j = 1/`ntr' {
  local kk = mod(`j'-1, `K') + 1
  local on : word `kk' of `dons'
  replace grp = `on' in `j'
}
gen double ttre = grp
gen double ai = rnormal()
keep unit ttre ai
expand `dT'
bysort unit: gen int t = 2000 + _n - 1
gen double tfe = 0.2*(t - 2000)
gen byte post = (t >= ttre) & !missing(ttre)
gen int h = t - ttre
gen double truetau = cond(post, `dtau'*(1 + 0.3*(h)), 0)
gen double region = mod(unit, 4) + 1
gen byte sg = mod(unit, 2)
gen double y = ai + tfe + truetau + `dsd'*rnormal()
label var y "outcome"
