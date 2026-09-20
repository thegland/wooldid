clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* Same simultaneous-band coverage question at FEW clusters (G=12), where the
* t(G-1) vs normal critical-value gap is 8-15%.  The band's critical value is
* simulated from NORMAL draws (wooldid.ado:16725 rnormal) while the pointwise
* intervals printed beside it use invttail(G-1).
local R 400
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
mkB 1
replace y = base
qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3) level(95) simultaneousbands
mat TB = e(b_es_att)
local tn : colnames TB
tempname fh
file open `fh' using "$WD/mc_band10.txt", write replace text
file write `fh' "draw|lvl|horizon|lb|ub|truth|crit" _n
forvalues c = 1/`=colsof(TB)' {
  local nm : word `c' of `tn'
  file write `fh' "0|.|`nm'|.|.|" (TB[1,`c']) "|." _n
}
forvalues r = 1/`R' {
  mkB `=1700000+`r''
  foreach L in 90 95 {
    cap qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3) level(`L') simultaneousbands
    if _rc==0 {
      cap confirm matrix e(esband_att)
      if _rc==0 {
        mat BB = e(esband_att)
        local bc : colnames BB
        local cr = e(esband_crit_att)
        forvalues h = 1/`=colsof(BB)' {
          local hn : word `h' of `bc'
          file write `fh' "`r'|`L'|`hn'|" (BB[1,`h']) "|" (BB[2,`h']) "|.|" (`cr') _n
        }
        * pointwise analytic rows for the same horizons, tagged lvl+1000
        mat MR = e(wdidmainresults)
        local rn : rownames MR
        local re : roweq MR
        forvalues i = 1/`=rowsof(MR)' {
          local rr : word `i' of `rn'
          local ee : word `i' of `re'
          if "`ee'" == "ES_att" {
            file write `fh' "`r'|`=`L'+1000'|`rr'|" (MR[`i',5]) "|" (MR[`i',6]) "|.|." _n
          }
        }
      }
    }
  }
}
file close `fh'
dis "MC_BAND10_DONE"
