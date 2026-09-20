clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 909090
set sortseed 909090

cap prog drop mknoisy
prog def mknoisy
  args nu nt
  clear
  set obs `=4*`nu'*`nt''
  gen long unit = ceil(_n/`nt')
  gen int  t    = _n - (unit-1)*`nt'
  gen byte i    = ceil(unit/`nu')
  gen double ttre = .
  qui replace ttre = 4 if i == 1
  qui replace ttre = 6 if i == 2
  gen double dose = 1 + 2*runiform()
  gen double ufe = 0
  bysort unit (t): replace ufe = rnormal(0,2) if _n==1
  bysort unit (t): replace ufe = ufe[1]
  gen double eps = rnormal(0,1)
  gen double eff = 0
  qui replace eff = (1.0 + 0.6*dose) if ttre<. & t>=ttre
  gen double y = 5 + ufe + 0.4*t + 0.9*dose + eff + eps
end

foreach nu in 5 12 30 {
 foreach ln in 2 3 5 {
  foreach vc in "robust" "cluster(unit)" {
    set seed 909090
    qui mknoisy `nu' 8
    cap noisily qui wooldid y i t ttre , att contreat(dose) contreatcontrols(basic) lattice(all `ln') `vc' regtol(14)
    local rc = _rc
    if `rc'==0 {
      mat M = e(wdidmainresults)
      local se_att = M[1,2]
      local se_ame = M[2,2]
      dis "LAT|nu=`nu'|lat=`ln'|`vc'|N=" e(N) "|df_r=" e(df_r) "|rss=" %10.4e e(rss) "|rmse=" %10.4e e(rmse) "|r2=" %8.6f e(r2) "|se_att=" %12.6g `se_att' "|se_ame=" %12.6g `se_ame'
    }
    else dis "LAT|nu=`nu'|lat=`ln'|`vc'|RC=`rc'"
  }
 }
}
dis "DONE"
