clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
global WD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/mc-ci-coverage"
* ===========================================================================
* probeF: what does the Ibragimov-Muller cross-cluster standard error under
* `impvalues' actually measure?
*
* If the per-cluster estimates were separate per-cluster fits (as IM requires),
* their dispersion would SHRINK with the idiosyncratic error.  If instead they
* are deterministic re-weightings of ONE pooled coefficient vector, their
* dispersion is driven by cluster COMPOSITION and the true cell effects, and
* is (nearly) invariant to the noise level.
*
* Two geometries:
*   A. i = cg (cohort) -- SEVERAL clusters share each level of i
*   B. i = unit        -- each cluster nests whole levels of i
* Noise sd swept over 0.60, 0.06, 0.006.
* ===========================================================================
cap prog drop imrep
prog def imrep
  args tag
  mat Q = e(wdidmainresults)
  local rn : rownames Q
  local re : roweq Q
  forvalues i = 1/`=rowsof(Q)' {
    local rr : word `i' of `rn'
    local ee : word `i' of `re'
    dis "IMROW|`tag'|`ee':`rr'|b=" %14.9f Q[`i',1] "|se=" %14.9f Q[`i',2] "|lb=" %14.9f Q[`i',5] "|ub=" %14.9f Q[`i',6]
  }
  dis "IMMETA|`tag'|df_r=" e(df_r) "|G=" e(N_clust) "|N=" e(N)
end

foreach SIG in 0.60 0.06 0.006 {
  foreach COMP in balanced unbalanced {
    clear
    set seed 4242
    set sortseed 4242
    set obs 120
    gen long unit = _n
    gen byte cg = mod(_n-1,4)+1
    if "`COMP'" == "balanced" {
      gen int reg = mod(_n-1,20)+1
    }
    if "`COMP'" == "unbalanced" {
      gen int reg = 1 + int(runiform()*20)
      replace reg = 20 if reg > 20
    }
    gen double ai = rnormal(0,1)
    expand 8
    bys unit: gen int t = 1999 + _n
    gen double ttre = .
    replace ttre = 2002 if cg==1
    replace ttre = 2004 if cg==2
    replace ttre = 2005 if cg==3
    gen double gt = 0.4*(t-1999)
    gen double truetau = 0
    replace truetau = 1.0 + 0.30*(t-ttre) if cg==1 & t>=ttre
    replace truetau = 0.6 + 0.15*(t-ttre) if cg==2 & t>=ttre
    replace truetau = 1.4 - 0.10*(t-ttre) if cg==3 & t>=ttre
    gen double y = ai + gt + truetau + rnormal(0,`SIG')
    qui tab reg cg
    cap noi wooldid y cg t ttre, att cluster(reg) impvalues level(95)
    dis "RC|F_icg_`COMP'_`SIG'|" _rc
    if _rc==0 imrep F_icg_`COMP'_`SIG'
    cap noi wooldid y unit t ttre, att cluster(reg) impvalues level(95)
    dis "RC|F_iunit_`COMP'_`SIG'|" _rc
    if _rc==0 imrep F_iunit_`COMP'_`SIG'
    * reference: analytic clustered SE at the same geometry, same data
    cap noi wooldid y cg t ttre, att cluster(reg) level(95)
    if _rc==0 {
      mat Q = e(wdidmainresults)
      dis "ANALYTIC|F_icg_`COMP'_`SIG'|b=" %14.9f Q[1,1] "|se=" %14.9f Q[1,2]
    }
  }
}
dis "PROBEF_DONE"
