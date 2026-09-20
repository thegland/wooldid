*==============================================================================
* mkdata.do -- tiny self-contained staggered-DiD generator for the
*               anticipation/replay/level probes.  No repo data, no repo code.
*
* Knobs (set as locals BEFORE `include'-ing this file):
*   dN      units                            [120]
*   dT      calendar periods                 [12]
*   dt0     first calendar period            [2000]
*   dstep   spacing between periods          [1]
*   donsets onset periods (space separated)  ["2004 2006 2008"]
*   dnevsh  never-treated share              [0.30]
*   dysd    idiosyncratic sd                 [0.30]
*   dantic  true effects start # early       [0]
*   dseed   seed                             [20260917]
* Produces: y i t ttre x1 x2 sg w dose
*==============================================================================
if "`dN'"      == "" local dN     120
if "`dT'"      == "" local dT     12
if "`dt0'"     == "" local dt0    2000
if "`dstep'"   == "" local dstep  1
if `"`donsets'"' == "" local donsets "2004 2006 2008"
if "`dnevsh'"  == "" local dnevsh .30
if "`dysd'"    == "" local dysd   .30
if "`dantic'"  == "" local dantic 0
if "`dseed'"   == "" local dseed  20260917

clear
set seed `dseed'
set sortseed `dseed'

local nco : word count `donsets'
quietly {
    set obs `dN'
    gen long i = _n
    gen double u = runiform()
    * cohort assignment: first dnevsh share never treated, rest split evenly
    gen byte cohort = .
    replace cohort = 0 if u < `dnevsh'
    gen double v = runiform()
    forvalues c = 1/`nco' {
        replace cohort = `c' if cohort == . & v <= `c'/`nco' & cohort == .
    }
    replace cohort = `nco' if cohort == .
    drop u v
    gen double ttre = .
    forvalues c = 1/`nco' {
        local oc : word `c' of `donsets'
        replace ttre = `oc' if cohort == `c'
    }
    gen double a_i   = rnormal()
    gen double taumu = 0.5 + 0.25*cohort
    gen byte   sg    = mod(i,2)
    gen double dose  = cond(cohort==0, 0, 0.5 + 0.5*mod(i,3))
    expand `dT'
    bysort i: gen int k = _n
    gen double t = `dt0' + (k-1)*`dstep'
    drop k
    gen double x1 = rnormal()
    gen double x2 = runiform()
    gen double w  = 1 + mod(i,4)
    gen double tsh = .
    * common period shock
    bysort t (i): replace tsh = rnormal() if _n == 1
    bysort t (i): replace tsh = tsh[1]
    gen byte trt = (ttre < . & t >= ttre - `dantic')
    gen double relt = t - ttre
    gen double y = a_i + tsh + 0.5*x1 - 0.25*x2 ///
                 + trt*taumu*(1 + 0.15*cond(relt<.,max(relt+`dantic',0),0)) ///
                 + `dysd'*rnormal()
    compress
}
