*=========================================================================
* repro_S2_band.do -- SELF-CONTAINED.
*
* CLAIM DEMONSTRATED: the analytic `simultaneousbands' critical value is the
* level() quantile of the sup of a STANDARD NORMAL vector drawn from the
* estimated correlation matrix, while the quantity it has to bound is a sup
* of t-type ratios on e(df_r) = G-1 degrees of freedom.  The band is about
* 22% too narrow at G = 12 and its measured simultaneous coverage is ~0.79
* at a nominal 0.90, while the POINTWISE intervals on the same runs are
* correctly calibrated.
*
* The DGP makes the true effect a function of the RELATIVE period only,
* identically across cohorts and units, so every weighting scheme gives the
* same horizon truth and the coverage measurement carries no estimand
* ambiguity: pre-horizon truth is 0, post-horizon truth is 0.60 / 0.90.
*
* ~8 minutes.
* Run: "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" -e do repro_S2_band.do
*=========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

local R 200

cap prog drop mkdat
prog def mkdat
  args nunit sd
  clear
  set seed `sd'
  set sortseed `sd'
  set obs `nunit'
  gen long unit = _n
  gen byte cg = mod(_n-1,4)+1
  gen double ttre = .
  replace ttre = 2004 if cg==1
  replace ttre = 2006 if cg==2
  replace ttre = 2008 if cg==3
  gen double ui = rnormal(0,0.70)
  expand 12
  bys unit: gen int t = 2000 + _n - 1
  gen double e = .
  bys unit (t): replace e = rnormal(0,0.55) if _n==1
  bys unit (t): replace e = 0.4*e[_n-1] + rnormal(0,0.55) if _n>1
  gen double rel = t - ttre
  gen double truetau = 0
  replace truetau = 0.60 if rel==0
  replace truetau = 0.90 if rel==1
  replace truetau = 1.10 if rel>=2 & rel<.
  gen double y = ui + 0.25*(t-2000) + truetau + e
end

* one diagnostic run
mkdat 12 55001
wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(2) simultaneousbands level(90)
di as txt _n "engine = " e(esband_engine) "   crit = " e(esband_crit_att) ///
             "   level = " e(esband_level) "   df_r = " e(df_r)
di as txt "for reference: invttail(11,.05) = " invttail(11,.05) ///
          "   invnormal(.95) = " invnormal(.95)

tempname M
matrix `M' = J(`R', 14, .)
forvalues r = 1/`R' {
  mkdat 12 `=55000+`r''
  cap qui wooldid y cg t ttre, att cluster(unit) esprelength(3) espostlength(2) simultaneousbands level(90)
  if _rc == 0 {
    cap mat BB = e(esband_att)
    if _rc == 0 {
      matrix `M'[`r',13] = e(esband_crit_att)
      matrix `M'[`r',14] = e(df_r)
      mat T = e(wdidmainresults)
      local ncol = colsof(BB)
      forvalues j = 1/`ncol' {
        local relv = BB[3,`j']
        local tr = 0
        if `relv' == 0 local tr = 0.60
        if `relv' == 1 local tr = 0.90
        local bb = .
        local ss = .
        forvalues k = 1/`=rowsof(T)' {
          if T[`k',7] == `relv' & T[`k',7] < .  {
            local bb = T[`k',1]
            local ss = T[`k',2]
          }
        }
        * store per-horizon: studentized deviation and band-cover flag
        matrix `M'[`r',`j']     = abs(`bb' - `tr')/`ss'
        matrix `M'[`r',`=`j'+6'] = (BB[1,`j'] <= `tr' & `tr' <= BB[2,`j'])
      }
    }
  }
}
clear
qui svmat double `M'
qui gen byte allcov = 1
qui gen double maxstud = 0
forvalues j = 1/6 {
  qui replace allcov = 0 if `M'`=`j'+6' == 0
  qui replace maxstud = max(maxstud, `M'`j') if `M'`j' < .
}
qui drop if `M'13 >= .
qui count
local R2 = r(N)
qui su allcov
local simcov = r(mean)
qui su `M'13
local mcrit = r(mean)
_pctile maxstud, p(90 95)
local o90 = r(r1)
local o95 = r(r2)
di as txt _n(2) "================= repro_S2 SUMMARY ================="
di as txt "usable draws                                         = " `R2'
di as txt "measured SIMULTANEOUS coverage at nominal 0.90       = " %6.3f `simcov'
di as txt "  3 binomial SE                                      = " %6.3f 3*sqrt(.9*.1/`R2')
di as txt "mean reported critical value e(esband_crit_att)      = " %8.4f `mcrit'
di as txt "ORACLE crit (q90 of max_h |b_h - truth_h| / se_h)    = " %8.4f `o90'
di as txt "  the band is this fraction of the width it needs    = " %6.3f `mcrit'/`o90'
di as txt "for reference: invnormal(.95) = " %8.4f invnormal(.95) ///
          "  invttail(11,.05) = " %8.4f invttail(11,.05)
* pointwise coverage per horizon, at the same level, on the same runs
di as txt _n "pointwise coverage per banded horizon at invttail(df_r,.05):"
forvalues j = 1/6 {
  qui count if `M'`j' < .
  if r(N) > 0 {
    qui count if `M'`j' <= invttail(11,.05) & `M'`j' < .
    local nn = r(N)
    qui count if `M'`j' < .
    di as txt "   horizon slot `j' : " %6.3f `nn'/r(N) "   (n = " r(N) ")"
  }
}
di as txt _n "ALLDONE_S2"
