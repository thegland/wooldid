*==========================================================================
* r5.do — WHY did r1 (micro panel, i = state) stay at the floor while r3
*         (cohort panel, i = cohort) moved by 6.6e-06?
*
* Same micro panel as r1.  Two knobs, one at a time:
*   LADDER  the control's level c (0 / 50 / 5000) with i = state
*   COARSE  the SAME data, the SAME controls, but i = the onset cohort (3 values)
* If the ladder is clean at i = state and dirty at i = cohort, the trigger is the
* granularity of i (how many rows share a collapse cell), not the control's scale.
*==========================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 71113
set sortseed 71113

local NS 30
local NI 25
local T  8
set obs `NS'
gen int state = _n
gen double u_s = rnormal()
gen int ttre = .
replace ttre = 2014 if state <= 10
replace ttre = 2016 if state > 10 & state <= 20
expand `NI'
bysort state: gen int pid = _n
gen double u_i = rnormal()
expand `T'
bysort state pid: gen int year = 2010 + _n
gen double age = 25 + 20*runiform()
gen double z   = rnormal()
gen byte D = (ttre < . & year >= ttre)
gen double y = 2 + u_s + u_i + 0.05*(year-2010) + 0.4*D + 0.01*age + 0.10*z + 0.5*rnormal()
* the onset cohort: 3 values (2014, 2016, never)
gen byte coh = 3
replace coh = 1 if ttre==2014
replace coh = 2 if ttre==2016
compress
save "_r5.dta", replace

mata:
void cmpm(string scalar an, string scalar bn) {
  real matrix A, B, S, Dm
  A = st_matrix(an); B = st_matrix(bn)
  if (rows(A)!=rows(B) | cols(A)!=cols(B)) {
    st_numscalar("mxall", .); st_numscalar("mxdia", .); return
  }
  S  = abs(A); S = S :* (S:>=abs(B)) + abs(B) :* (S:<abs(B))
  Dm = abs(A-B) :/ (S :+ (S:==0))
  st_numscalar("mxall", max(Dm))
  if (rows(A)==cols(A)) st_numscalar("mxdia", max(diagonal(Dm)))
  else st_numscalar("mxdia", .)
}
end

dis _n as text "{hline 92}"
dis as text "ivar   c        maxrel V(all)   maxrel V(diag)   maxrel main     cells"
dis as text "{hline 92}"
foreach iv in state coh {
  foreach c in 0 50 5000 {
    use "_r5.dta", clear
    gen double xc = `c' + 8*z
    cap noi qui wooldid y `iv' year ttre, att att_i att_it esprelength(2) espostlength(2) controls(xc)
    local rcA = _rc
    if `rcA'==0 {
      matrix VA = e(V)
      matrix MA = e(wdidmainresults)
      local nc `e(wdid_fcncell)'
    }
    use "_r5.dta", clear
    gen double xc = `c' + 8*z
    cap noi qui wooldid y `iv' year ttre, att att_i att_it esprelength(2) espostlength(2) controls(xc) nofastcollapse
    local rcB = _rc
    if `rcB'==0 {
      matrix VB = e(V)
      matrix MB = e(wdidmainresults)
    }
    if `rcA'!=0 | `rcB'!=0 {
      dis as error "`iv' c=`c' : rcA=`rcA' rcB=`rcB'"
      continue
    }
    mata: cmpm("VA","VB")
    scalar sVall = mxall
    scalar sVdia = mxdia
    mata: cmpm("MA","MB")
    scalar sMall = mxall
    dis as result %-6s "`iv'" " c=" %-7.0g `c' "  " %9.2e sVall "       " %9.2e sVdia "      " %9.2e sMall "    `nc'"
  }
}
dis as text "{hline 92}"
dis "R5_DONE"
