clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* ============================================================================
* probeE: does a simultaneous band ever come out NARROWER than the pointwise
* interval printed beside it?  The help (sthlp, simultaneousbands) states that
* under boottstat "the band can never come out narrower than the interval it
* accompanies."  The band is SYMMETRIC (b +/- c*se_obs) while the boottstat
* pointwise interval is an EQUAL-TAILED rank interval in the studentized draws,
* so a skewed t-draw distribution can make the pointwise interval wider.
* Setting: 12 units (few clusters), 4 cohorts, ONE variance-outlier unit.
* ============================================================================
cap prog drop bandcmp
prog def bandcmp
  args tag
  local eng = e(esband_engine)
  mat MR = e(wdidmainresults)
  local rn : rownames MR
  local re : roweq MR
  local cn : colnames MR
  local j = 0
  local kl = 0
  foreach c of local cn {
    local ++j
    if substr("`c'",1,3)=="lb_" {
      local kl = `j'
    }
    if substr("`c'",1,3)=="ub_" {
      local ku = `j'
    }
  }
  foreach f in att att_it {
    cap confirm matrix e(esband_`f')
    if _rc == 0 {
      mat BB = e(esband_`f')
      local bc : colnames BB
      local crit = e(esband_crit_`f')
      forvalues h = 1/`=colsof(BB)' {
        local hn : word `h' of `bc'
        local blo = BB[1,`h']
        local bhi = BB[2,`h']
        local wband = `bhi'-`blo'
        * find the matching pointwise row  ES_<f>_<hn>
        local wpt = .
        local plo = .
        local phi = .
        forvalues i = 1/`=rowsof(MR)' {
          local rr : word `i' of `rn'
          local ee : word `i' of `re'
          if "`ee'" == "ES_`f'" & "`rr'" == "`hn'" {
            local plo = MR[`i',`kl']
            local phi = MR[`i',`ku']
            local wpt = `phi'-`plo'
          }
        }
        if `wpt' < . & `wband' < . {
          local ratio = `wband'/`wpt'
          dis "BANDCMP|`tag'|`eng'|`f'|`hn'|crit=" %9.6f `crit' "|wband=" %10.6f `wband' "|wpt=" %10.6f `wpt' "|ratio=" %8.5f `ratio' "|" cond(`wband' < `wpt' - 1e-10, "NARROWER", "ok")
        }
      }
    }
  }
end

forvalues s = 1/12 {
  clear
  set seed `=90000+`s''
  set sortseed `=90000+`s''
  set obs 12
  gen long unit = _n
  gen byte cg = mod(_n-1,4)+1
  gen double ai = rnormal(0,1)
  gen double sdmult = cond(_n==1, 6, 1)
  expand 10
  bys unit: gen int t = 1999 + _n
  gen double ttre = .
  replace ttre = 2002 if cg==1
  replace ttre = 2004 if cg==2
  replace ttre = 2006 if cg==3
  gen double gt = 0.35*(t-1999)
  gen double truetau = 0
  replace truetau = 1.0 + 0.25*(t-ttre) if cg<4 & t>=ttre & ttre<.
  gen double y = ai + gt + truetau + rnormal(0,0.5)*sdmult
  foreach cfg in stud plain {
    if "`cfg'"=="stud"  local extra "boottstat"
    if "`cfg'"=="plain" local extra ""
    cap noi wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3) ///
        bootstrap(standard) reps(299) rseed(`=700000+`s'') level(95) ///
        simultaneousbands `extra' dots(0)
    dis "RC|E_`s'_`cfg'|" _rc
    if _rc == 0 {
      bandcmp E_`s'_`cfg'
    }
  }
}
dis "PROBEE_DONE"
