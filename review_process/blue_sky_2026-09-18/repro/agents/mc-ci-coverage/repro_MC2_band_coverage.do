*=========================================================================
* SELF-CONTAINED REPRO  MC2 -- the analytic simultaneous (sup-t) band
* under-covers because its critical value is simulated from NORMAL draws
* while the pointwise intervals printed beside it use invttail(G-1).
* 12 clusters, 6 horizons, 200 Monte Carlo replications (~4 min).
*
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_MC2_band_coverage.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
cap prog drop mkB
prog def mkB
  args sd
  clear
  set obs 12
  gen long unit = _n
  gen byte cg = mod(_n-1,4)+1
  gen double ai = 0.9*mod(_n*13,17)/17 - 0.45
  expand 10
  bys unit: gen int t = 1999 + _n
  gen double ttre = .
  replace ttre = 2002 if cg==1
  replace ttre = 2004 if cg==2
  replace ttre = 2005 if cg==3
  gen double truetau = 0
  replace truetau = 1.0 + 0.25*(t-ttre) if cg<4 & t>=ttre & ttre<.
  gen double base = ai + 0.35*(t-1999) + truetau
  set seed `sd'
  set sortseed `sd'
  by unit: gen double ui = rnormal(0,0.40) if _n==1
  by unit: replace ui = ui[1]
  gen double y = base + ui + rnormal(0,0.60)
end
* target: the noiseless event-study coefficients
mkB 1
replace y = base
qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3) level(95) simultaneousbands
mat TRUE = e(b_es_att)
local hn : colnames TRUE
di as txt "design truth of the 6 horizons (`hn'):"
mat list TRUE, format(%12.8f)

local R 200
foreach L in 90 95 {
  local nband 0
  local npt 0
  local n 0
  local scrit 0
  forvalues r = 1/`R' {
    mkB `=1700000+`r''
    qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3) level(`L') simultaneousbands
    cap confirm matrix e(esband_att)
    if _rc continue
    mat BB = e(esband_att)
    mat MR = e(wdidmainresults)
    local scrit = `scrit' + e(esband_crit_att)
    local ++n
    local bandok 1
    local ptok 1
    local bc : colnames BB
    local rn : rownames MR
    local re : roweq MR
    forvalues h = 1/`=colsof(BB)' {
      local hh : word `h' of `bc'
      local tt = TRUE[1,`h']
      if !(BB[1,`h'] <= `tt' & `tt' <= BB[2,`h']) local bandok 0
      forvalues i = 1/`=rowsof(MR)' {
        local rr : word `i' of `rn'
        local ee : word `i' of `re'
        if "`ee'" == "ES_att" & "`rr'" == "`hh'" {
          if !(MR[`i',5] <= `tt' & `tt' <= MR[`i',6]) local ptok 0
        }
      }
    }
    local nband = `nband' + `bandok'
    local npt   = `npt'   + `ptok'
  }
  di ""
  di as txt "level(`L'), G=12, 6 horizons, R=`n'"
  di as txt "  SIMULTANEOUS band covers all 6 horizons : " as res %6.3f (`nband'/`n') ///
     as txt "   (nominal " as res %5.2f (`L'/100) as txt ")"
  di as txt "  the 6 POINTWISE intervals all cover     : " as res %6.3f (`npt'/`n') ///
     as txt "   (no claim attaches to this; shown for contrast)"
  di as txt "  mean band critical value                : " as res %9.5f (`scrit'/`n')
  di as txt "  invttail(11,(100-`L')/200)              : " as res %9.5f invttail(11,(100-`L')/200)
  di as txt "  invnormal(1-(100-`L')/200)              : " as res %9.5f invnormal(1-(100-`L')/200)
}
