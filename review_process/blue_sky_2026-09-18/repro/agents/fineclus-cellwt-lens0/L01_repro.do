*! L01 — INDEPENDENT REPRODUCTION of FCW-01, third lens, from the claim alone.
*! Own panel: hospital systems (cluster) x hospitals (i), T=7, staggered system-level
*! adoption. Arm B = arm A with ONE (unit,period) cell of the system variable changed:
*! hospital 1 is acquired by system 2 for its final period. Nothing else differs.
*! Measures: se/b for all four estimands, the nesting predicate computed in-sample,
*! the exact hex of b(att_i), and the NAMES of every e() scalar/macro/matrix in each arm.
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set linesize 200
which wooldid

cap prog drop mkhosp
prog def mkhosp
  syntax , dseed(int)
  * 12 systems x 5 hospitals x 7 years. 8 systems adopt (4 at t=4, 4 at t=6),
  * 4 never adopt. Effect heterogeneous ACROSS systems and WITHIN them.
  clear
  set seed `dseed'
  set sortseed `dseed'
  set obs 12
  gen long sys = _n
  gen double a_g = rnormal(0,.5)
  gen onset = .
  replace onset = 4 if sys <= 4
  replace onset = 6 if sys > 4 & sys <= 8
  expand 5
  sort sys
  gen long hosp = _n
  gen double b_i = rnormal(0,.3)
  gen double ufe = rnormal(0,1)
  expand 7
  bysort hosp (sys): gen int t = _n
  gen double tau = 1 + a_g + b_i
  gen byte D = (onset < . & t >= onset)
  gen double y = 3 + ufe + .15*t + tau*D + rnormal(0,.8)
  keep y hosp sys t onset
  order hosp sys t onset y
end

cap prog drop grab
prog def grab
  args arm ds rs nst
  mat b = e(b)
  mat V = e(V)
  local out "`arm',`ds',`rs',`nst'"
  foreach nm in att att_i att_it att_itime {
    local c = colnumb(b,"main_`nm'")
    if `c' < . {
      local bb = string(b[1,`c'],"%20.0g")
      local ss = string(sqrt(V[`c',`c']),"%20.0g")
      local out "`out',`bb',`ss'"
    }
    else {
      local out "`out',.,."
    }
  }
  file open fh using "L01.csv", write append
  file write fh "`out'" _n
  file close fh
  local ca = colnumb(b,"main_att_i")
  dis "HEX_`arm'_`ds'_`rs' b_att_i = " %21x b[1,`ca'] "  se_att_i = " %21x sqrt(V[`ca',`ca'])
end

file open fh using "L01.csv", write replace
file write fh "arm,ds,rs,nest,b_att,se_att,b_atti,se_atti,b_attit,se_attit,b_attitime,se_attitime" _n
file close fh

foreach ds in 7701 7702 7703 {
 foreach arm in A B {
  qui mkhosp , dseed(`ds')
  if "`arm'" == "B" qui replace sys = 2 if hosp == 1 & t == 7
  * predicate, computed exactly as wooldid_bootrun does (min==max of cluster by i)
  tempvar cmn cmx
  qui gegen double `cmn' = min(sys), by(hosp)
  qui gegen double `cmx' = max(sys), by(hosp)
  qui count if `cmn' != `cmx'
  local nstrad = r(N)
  cap assert `cmn' == `cmx'
  local nst = (_rc == 0)
  drop `cmn' `cmx'
  dis "PRED arm=`arm' ds=`ds' nest=`nst' straddling_rows=`nstrad'"
  foreach rs in 911 912 {
    dis "=== BEGINRUN arm=`arm' ds=`ds' rs=`rs' ==="
    cap noi wooldid y hosp t onset, att att_i att_it att_itime cluster(sys) ///
        bootstrap(standard) reps(199) rseed(`rs') dots(0)
    dis "RC arm=`arm' ds=`ds' rs=`rs' = " _rc
    if _rc == 0 {
      grab `arm' `ds' `rs' `nst'
      * e() surface: names only, so the two arms can be diffed for any marker
      if `ds' == 7701 & `rs' == 911 {
        file open gh using "esurf_`arm'.txt", write replace
        file write gh "SCALARS" _n
        foreach s in `: e(scalars)' {
          file write gh "`s'" _n
        }
        file write gh "MACROS" _n
        foreach m in `: e(macros)' {
          file write gh "`m' = `e(`m')'" _n
        }
        file write gh "MATRICES" _n
        foreach x in `: e(matrices)' {
          file write gh "`x'" _n
        }
        file close gh
      }
    }
    dis "=== ENDRUN arm=`arm' ds=`ds' rs=`rs' ==="
  }
 }
}

dis "ALLDONE"
