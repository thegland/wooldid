clear all
set more off
set processors 2
set maxvar 20000
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
run "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat-lens3/mkfix.do"
set seed 70007
set sortseed 70007
mkfix . 0
gen byte post = (ttre<. & year>=ttre)
qui egen cell = group(st year) if post==1
qui replace cell = 0 if cell>=.
qui xtile lat = dose_it , nq(3)
gen double ringer = dose_it
gen byte nullv = 1
qui reg y_dose i.st i.year i.nullv##c.ringer i.lat##c.ringer ib0.cell##i.lat##c.dose_it
dis "ZOOM2 regress(ols) rank=" e(rank) " df_r=" e(df_r) " params=" e(N)-e(df_r) " rss=" %16.8f e(rss)
mat bR = e(b)
qui reghdfe y_dose ib0.cell##i.lat##c.dose_it , absorb(i.st i.year i.nullv##c.ringer i.lat##c.ringer)
dis "ZOOM2 reghdfe(ols) rank=" e(rank) " df_r=" e(df_r) " rss=" %16.8f e(rss)
mat bH = e(b)
* names present (non-omitted) in each
local nR : colnames bR
local nH : colnames bH
local cR = colsof(bR)
local cH = colsof(bH)
dis "ZOOM2 regress colsof=`cR'   reghdfe colsof=`cH'"
* count omitted (o.) entries in regress
local omitR = 0
forval j = 1/`cR' {
  local nm : word `j' of `nR'
  if strpos("`nm'","o.")==1 | strpos("`nm'",".o.")>0 local ++omitR
}
dis "ZOOM2 regress omitted-coef count=`omitR'"
* find slope columns present in reghdfe but omitted/absent in regress
local shown = 0
forval j = 1/`cH' {
  local nm : word `j' of `nH'
  if strpos("`nm'","dose_it")>0 {
    cap local v = bR[1,"`nm'"]
    if _rc != 0 {
      local ++shown
      if `shown' <= 12 dis "ZOOM2 in reghdfe only: `nm'  b=" %14.6g bH[1,`j']
    }
  }
}
dis "ZOOM2 dose slope cols in reghdfe but NOT in regress e(b): `shown'"
dis "DONE_E"
