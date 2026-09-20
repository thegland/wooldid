clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917

* --- is the "two digits" check on the TEXT or on the VALUE? ---
set level 95
cap noi set level 99.90000000000001
dis "F_SETLIT_RC=" _rc " c(level)=" c(level)
cap noi set level 99.9000
dis "F_SETTRAILZERO_RC=" _rc " c(level)=" c(level)
set level 95
cap noi set level 99.900
dis "F_SET3DEC_RC=" _rc " c(level)=" c(level)
set level 95

* --- reproduce the two-call sequence with full instrumentation ---
cap prog drop mkexact
prog def mkexact
  clear
  set seed 11
  set sortseed 11
  qui set obs 40
  gen long unit = _n
  gen double a_i = rnormal()
  gen double u = runiform()
  gen double ttre = .
  qui replace ttre = 2003 if u > 0.35 & u <= 0.5667
  qui replace ttre = 2005 if u > 0.5667 & u <= 0.7833
  qui replace ttre = 2007 if u > 0.7833
  qui expand 10
  bysort unit: gen int t = 1999 + _n
  gen byte D = (ttre < . & t >= ttre)
  gen double y = a_i + 0.4*(t-2000) + 2*D + 0.3*rnormal()
  gen long i = unit
  drop u a_i D
end
mkexact

macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
set level 99.9
dis "F_C0=" c(level)
set trace on
set tracedepth 2
cap noi qui wooldid y i t ttre , att robust level(95) safety(off)
set trace off
dis "F_CALL1_RC=" _rc " c=" c(level) " L0=[${WOOLDID_LEVEL0}]"
set trace on
set tracedepth 2
cap noi qui wooldid y i t ttre , att robust safety(off)
set trace off
dis "F_CALL2_RC=" _rc " c=" c(level) " L0=[${WOOLDID_LEVEL0}]"
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
set level 95

* --- plain call with a hostile session level, traced ---
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
set level 99.9
set trace on
set tracedepth 2
cap noi qui wooldid y i t ttre , att robust safety(off)
set trace off
dis "F_CALL3_RC=" _rc " c=" c(level) " L0=[${WOOLDID_LEVEL0}]"
macro drop WOOLDID_LEVEL0 WOOLDID_LASTRC
set level 95
dis "=================== DONE pF ==================="
