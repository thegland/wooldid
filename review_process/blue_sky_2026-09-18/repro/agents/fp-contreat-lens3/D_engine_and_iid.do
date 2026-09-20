* LENS3 probe D:
*  (1) is the lattice ame POINT ESTIMATE engine-dependent (reghdfe vs regress)?
*  (2) is the all-zero VCV specific to the sandwich VCE, and does the iid
*      default path deliver honest coverage where robust delivers nothing?
clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
run "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat-lens3/mkfix.do"

tempname P
postfile `P' int rep double aIID double sIID double aROB double sROB double aHC1 double sHC1 using Dres , replace
forval r = 1/25 {
  local s = 70000 + `r'
  set seed `s'
  set sortseed `s'
  mkfix . 0
  tempfile G
  qui save `G'
  local a1 = .
  local s1 = .
  local a2 = .
  local s2 = .
  local a3 = .
  local s3 = .
  use `G', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3)
  if _rc==0 {
    mat M = e(wdidmainresults)
    local a1 = M[2,1]
    local s1 = M[2,2]
  }
  use `G', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) robust
  if _rc==0 {
    mat M = e(wdidmainresults)
    local a2 = M[2,1]
    local s2 = M[2,2]
  }
  use `G', clear
  cap qui wooldid y_dose st year ttre , att contreat(dose_it) contreatcontrols(basic) lattice(all 3) customvce(robust)
  if _rc==0 {
    mat M = e(wdidmainresults)
    local a3 = M[2,1]
    local s3 = M[2,2]
  }
  local d = abs(`a2'-`a3')
  dis "D rep`r' IIDa=" %12.5f `a1' " IIDse=" %10.5f `s1' " | ROBa=" %12.5f `a2' " ROBse=" %10.5f `s2' " | HC1a=" %12.5f `a3' " HC1se=" %10.5f `s3' "  ENGINEGAP=" %12.6g `d'
  post `P' (`r') (`a1') (`s1') (`a2') (`s2') (`a3') (`s3')
}
postclose `P'
use Dres, clear
qui gen double gap = abs(aROB-aHC1)
qui gen double relgap = gap/max(abs(aROB),1)
qui count if relgap > 1e-8
dis "DSUM engine mismatch (rel>1e-8) in " r(N) " of 25 draws"
qui sum gap
dis "DSUM engine gap: max=" %14.6f r(max) " mean=" %14.6f r(mean)
qui count if sROB < .
dis "DSUM draws where DEFAULT+robust returned a finite ame se: " r(N) " of 25"
qui count if sIID < .
dis "DSUM draws where DEFAULT+iid returned a finite ame se: " r(N) " of 25"
qui gen byte covIID = (aIID-1.96*sIID<=0.4) & (aIID+1.96*sIID>=0.4) if sIID<.
qui gen byte covHC1 = (aHC1-1.96*sHC1<=0.4) & (aHC1+1.96*sHC1>=0.4) if sHC1<.
qui count if covIID==1
dis "DSUM IID  95% CI covers truth 0.400 in " r(N) " of 25"
qui count if covHC1==1
dis "DSUM HC1  95% CI covers truth 0.400 in " r(N) " of 25"
qui sum aIID
dis "DSUM IID ame  mean=" %12.5f r(mean) " sd=" %12.5f r(sd)
qui sum sIID
dis "DSUM IID reported se mean=" %12.5f r(mean) " p50 next"
qui sum sIID, detail
dis "DSUM IID reported se p50=" %12.5f r(p50)
qui sum sHC1
dis "DSUM HC1 reported se mean=" %12.5f r(mean)

* ---- zoom on the mismatching draw: ranks and dropped columns --------------
set seed 70007
set sortseed 70007
mkfix . 0
gen byte post = (ttre<. & year>=ttre)
qui egen cell = group(st year) if post==1
qui replace cell = 0 if cell>=.
qui xtile lat = dose_it , nq(3)
gen double ringer = dose_it
gen byte nullv = 1
qui reghdfe y_dose ib0.cell##i.lat##c.dose_it , absorb(i.st i.year i.nullv##c.ringer i.lat##c.ringer) vce(robust)
dis "ZOOM reghdfe rank=" e(rank) " df_r=" e(df_r) " rss=" %14.6f e(rss) " N=" e(N)
qui margins , dydx(dose_it) subpop(if post==1) noestimcheck nose
mat b = r(b)
dis "ZOOM reghdfe margins ame=" %16.8f b[1,1]
qui reg y_dose i.st i.year i.nullv##c.ringer i.lat##c.ringer ib0.cell##i.lat##c.dose_it , vce(robust)
dis "ZOOM regress rank=" e(rank) " df_r=" e(df_r) " rss=" %14.6f e(rss) " N=" e(N)
qui margins , dydx(dose_it) subpop(if post==1) noestimcheck
mat b = r(b)
mat Vm = r(V)
dis "ZOOM regress margins ame=" %16.8f b[1,1] "  se=" %14.6g sqrt(Vm[1,1])
dis "DONE_D"
