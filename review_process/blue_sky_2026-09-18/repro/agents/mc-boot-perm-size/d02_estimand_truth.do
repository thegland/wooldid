*======================================================================
* d02_estimand_truth.do -- DECISIVE test of the d01 observation:
*   requesting event-study horizons moves the headline `att'.
*   Here the DGP is EXACT (no noise at all) and the treatment effect is
*   CONSTANT = 1 on every treated cell, so EVERY sensible ATT estimand
*   equals exactly 1.  Any run that does not return 1 is returning the
*   wrong number, not a different-but-valid weighting.
* Self-contained: run with  stata-mp -e do d02_estimand_truth.do
*======================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid

cap prog drop mkexact
prog def mkexact
  syntax [, TAUMODE(string) SEEDX(integer 1)]
  set seed `seedx'
  set sortseed `seedx'
  clear
  set obs 10
  gen long st = _n
  gen double ttre = .
  replace ttre = 2004 if inrange(st,1,2)
  replace ttre = 2005 if inrange(st,3,4)
  replace ttre = 2006 if inrange(st,5,6)
  gen double ai = st*0.37
  expand 8
  bysort st: gen int t = 1999 + _n
  gen double ft = (t-1999)*0.21
  gen double h  = t - ttre
  gen byte post = (h >= 0) & (ttre < .)
  gen double tau = 0
  if "`taumode'" == "const"  replace tau = 1           if post
  if "`taumode'" == "ramp"   replace tau = 1 + 0.5*h   if post
  gen double y = ai + ft + tau          /* EXACT: no error term at all */
  keep st t ttre y h post tau
end

foreach tm in const ramp {
  dis "=================== TAUMODE = `tm' ==================="
  mkexact , taumode(`tm')
  qui sum tau if post
  dis "TRUE mean effect over treated cells (= true att with cell weights) = " %21x r(mean)
  dis "   in decimal: " r(mean)
  save d02_`tm'.dta, replace

  foreach spec in none w33 w99 efbp w33efbp {
    use d02_`tm'.dta, clear
    ereturn clear
    if "`spec'" == "none" cap noi qui wooldid y st t ttre, cluster(st) att att_i att_it att_itime
    if "`spec'" == "w33"  cap noi qui wooldid y st t ttre, cluster(st) att att_i att_it att_itime esprelength(3) espostlength(3)
    if "`spec'" == "w99"  cap noi qui wooldid y st t ttre, cluster(st) att att_i att_it att_itime esprelength(9) espostlength(9)
    if "`spec'" == "efbp" cap noi qui wooldid y st t ttre, cluster(st) att att_i att_it att_itime esfixedbaseperiod
    if "`spec'" == "w33efbp" cap noi qui wooldid y st t ttre, cluster(st) att att_i att_it att_itime esprelength(3) espostlength(3) esfixedbaseperiod
    dis "SPEC=`spec' rc=" _rc " N=" e(N)
    if _rc == 0 {
      mat M = e(wdidmainresults)
      local rn : rownames M
      local re : roweq M
      local nr : word count `rn'
      forvalues z = 1/`nr' {
        local nz : word `z' of `rn'
        local ez : word `z' of `re'
        if "`ez'" == "main" dis "   `tm' `spec'  main:`nz' = " %21x M[`z',1] "   (" M[`z',1] ")"
      }
    }
  }
}
