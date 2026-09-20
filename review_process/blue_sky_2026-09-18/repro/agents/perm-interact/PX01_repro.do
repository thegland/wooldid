*==============================================================================
* PX-01 self-contained repro.
* fastnose + permanalyticse under permute(): permanalyticse is a SILENT NO-OP and
* the run then asserts, in e(perm_analyticse), in the table header and in the
* footer, that the analytic standard errors and e(V)/e(V_es_*) are there.
* Contrast: fastnose + permtstat and fastnose + permtci are REFUSED r(198).
* No repo data, no test suite: the panel is generated here.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 20260916
set sortseed 20260916

* ---- a plain staggered panel: 40 units, 8 years, onsets 2003 and 2005 ----
clear
set obs 40
gen long unit = _n
gen double ttre = .
replace ttre = 2003 if unit <= 12
replace ttre = 2005 if unit > 12 & unit <= 24
expand 8
bysort unit: gen int t = 1999 + _n
gen double a_i = rnormal()
gen double g_t = 0.25*(t - 2002)
gen double y = a_i + g_t + 2.0*(ttre < . & t >= ttre) + rnormal()*0.6

local BASE att cluster(unit) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0)

cap prog drop show
prog def show
  syntax anything(name=tag) [, OPTS(string asis) BASE(string asis)]
  cap noi wooldid y unit t ttre , `base' `opts'
  dis "PX01|`tag'|RC=" _rc
  if _rc != 0 exit
  cap confirm matrix e(V)
  local hv = (_rc==0)
  cap confirm matrix e(V_es_att)
  local hves = (_rc==0)
  cap confirm matrix e(b_es_att)
  local hbes = (_rc==0)
  tempname M
  matrix `M' = e(wdidmainresults)
  dis "PX01|`tag'|e(perm_analyticse)=" e(perm_analyticse) "  vcvstate=`e(wdid_vcvstate)'  e(V)=`hv'  e(V_es_att)=`hves'  e(b_es_att)=`hbes'  se_analytic(row1)=" `M'[1,2]
end

dis _n "### 1. permanalyticse alone: SEs and e(V) present (the documented behaviour)"
show ase , base(`BASE') opts(permanalyticse)

dis _n "### 2. fastnose + permanalyticse: ACCEPTED, but permanalyticse did nothing"
show fn_ase , base(`BASE') opts(fastnose permanalyticse)
dis "### ^^ note e(perm_analyticse)=1 with NO e(V), NO e(V_es_att), NO e(b_es_att)"
dis "### and read the header line and the 'The standard-error column is the ANALYTIC'"
dis "### footer line printed just above: both are false on this run."

dis _n "### 3. the twin combinations that ARE refused, on identical reasoning"
show fn_ts  , base(`BASE') opts(fastnose permtstat)
show fn_tci , base(`BASE') opts(fastnose permtci)
dis "PX01_DONE"
