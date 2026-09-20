* gen_mypanel.do -- my own staggered-DiD generator (no repo data, no repo generator).
* knobs read from caller scope: mN mT mt0 monsets mnevshare mysd maisd mtsd mtau mdose mseed
if "`mseed'"      == "" local mseed 71717
if "`mN'"         == "" local mN 120
if "`mT'"         == "" local mT 10
if "`mt0'"        == "" local mt0 2000
if `"`monsets'"'  == "" local monsets "2003 2005 2007"
if "`mnevshare'"  == "" local mnevshare 0.35
if "`mysd'"       == "" local mysd 0.30
if "`maisd'"      == "" local maisd 1
if "`mtsd'"       == "" local mtsd 0.50
if "`mtau'"       == "" local mtau 2
clear
set seed `mseed'
set sortseed `mseed'
local mK : word count `monsets'
qui set obs `mN'
gen long unit = _n
gen double u = runiform()
* assign: never-treated share, then even split across cohorts
gen byte grp = 0
qui replace grp = ceil((u - `mnevshare')/((1-`mnevshare')/`mK') ) if u > `mnevshare'
qui replace grp = `mK' if grp > `mK'
qui replace grp = 0 if u <= `mnevshare'
gen double ons = .
forvalues k = 1/`mK' {
  local o : word `k' of `monsets'
  qui replace ons = `o' if grp == `k'
}
gen double a_i = rnormal(0,`maisd')
gen double dose_u = 0.5 + runiform()*2
qui expand `mT'
bysort unit: gen int tt = `mt0' + _n - 1
* common period shock
tempvar shk
qui gen double `shk' = .
forvalues j = 0/`=`mT'-1' {
  local s = rnormal(0,`mtsd')
  qui replace `shk' = `s' if tt == `mt0' + `j'
}
gen byte D = (ons < . & tt >= ons)
gen double x1 = rnormal(0,1)
gen double y = a_i + `shk' + `mtau'*D + 0.5*x1 + rnormal(0,`mysd') + 10
gen double dose = dose_u*D
gen double wt = 1 + runiform()*2
gen byte region = mod(unit,8) + 1
rename tt t
rename unit i
rename ons ttre
drop u grp a_i dose_u
order y i t ttre
