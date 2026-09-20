clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
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

* probeE2: same question at level(99) (where the equal-tailed pointwise endpoints are the
* 2nd order statistics of 199 studentized draws and the band's critical value is the 99th
* percentile of the max) and with an ASYMMETRIC outlier cluster.
forvalues s = 1/8 {
  clear
  set seed `=95000+`s''
  set sortseed `=95000+`s''
  set obs 12
  gen long unit = _n
  gen byte cg = mod(_n-1,4)+1
  gen double ai = rnormal(0,1)
  expand 10
  bys unit: gen int t = 1999 + _n
  gen double ttre = .
  replace ttre = 2002 if cg==1
  replace ttre = 2004 if cg==2
  replace ttre = 2006 if cg==3
  gen double gt = 0.35*(t-1999)
  gen double truetau = 0
  replace truetau = 1.0 + 0.25*(t-ttre) if cg<4 & t>=ttre & ttre<.
  * ASYMMETRIC contamination: unit 1 takes only large NEGATIVE shocks
  gen double e = rnormal(0,0.5)
  replace e = -abs(rnormal(0,4)) if unit == 1
  gen double y = ai + gt + truetau + e
  foreach L in 99 95 {
    cap noi wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3) ///
        bootstrap(standard) reps(199) rseed(`=750000+`s'') level(`L') ///
        simultaneousbands boottstat dots(0)
    dis "RC|E2_`s'_`L'|" _rc
    if _rc == 0 {
      bandcmp E2_`s'_`L'
    }
  }
  * bca on the same contaminated data: half-interval / guard check
  cap noi wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(3) ///
      bootstrap(standard) reps(199) rseed(`=750000+`s'') level(99) bootdisplay(bca) dots(0)
  dis "RC|E2_`s'_bca|" _rc
  if _rc == 0 {
    foreach fam in normal percentile bc bca {
      cap confirm matrix e(boot_ci_`fam')
      if _rc==0 {
        mat CC = e(boot_ci_`fam')
        mat BO = e(b)
        local cn : colnames CC
        forvalues k = 1/`=colsof(CC)' {
          local nm : word `k' of `cn'
          local lo = CC[1,`k']
          local hi = CC[2,`k']
          local bb = BO[1,`k']
          if (`lo'<. & `hi'>=.) | (`lo'>=. & `hi'<.) dis "VIOL_CIHALF|E2_`s'|`fam'|`nm'|lb=" `lo' "|ub=" `hi'
          if `lo'<. & `hi'<. {
            if !(`lo'<`hi') dis "VIOL_CILBUB|E2_`s'|`fam'|`nm'|" `lo' "|" `hi'
            if `bb'<. & (`bb'<`lo' | `bb'>`hi') dis "VIOL_BOUT|E2_`s'|`fam'|`nm'|b=" `bb' "|[" `lo' "," `hi' "]"
          }
        }
      }
    }
    cap confirm matrix e(boot_accel)
    if _rc==0 {
      mat list e(boot_accel), format(%14.9f)
    }
  }
}
dis "PROBEE2_DONE"
