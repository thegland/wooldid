clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
local DD "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/data"
local OO "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/out"
do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/nofastcollapse/lib_dump.do"

* ---- how does the controls() disagreement scale with ROWS PER CELL? ----
* same design, same DGP, only the number of units per cohort changes, so the
* number of rows the collapse folds into one grows 4 -> 8 -> 32 -> 128 -> 512.
foreach u in 4 8 32 128 512 {
  clear
  set seed 20260917
  set sortseed 20260917
  qui set obs `=4*`u''
  qui gen long unit   = _n
  qui gen byte cohort = mod(_n-1,4)+1
  qui gen double ttre = .
  qui replace ttre = 2004 if cohort==1
  qui replace ttre = 2006 if cohort==2
  qui replace ttre = 2008 if cohort==3
  qui gen double a_i = rnormal()
  qui expand 10
  qui bysort unit: gen int t = 2000 + _n - 1
  qui gen double x1 = rnormal()
  qui gen double y = a_i + 0.4*(t-2000) + 0.5*x1 + 1.2*(ttre<. & t>=ttre) + rnormal()*0.3
  qui save "`DD'/cs`u'.dta", replace
  wdcell , dta("`DD'/cs`u'.dta") tag(cs`u') base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(3) espostlength(3) jointtests controls(x1))
  wdcell , dta("`DD'/cs`u'.dta") tag(cs`u'_noc) base(wooldid y cohort t ttre) out("`OO'") opts(att esprelength(3) espostlength(3) jointtests)
}
dis "P10_DONE"
