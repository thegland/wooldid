*==============================================================================
* p12: the REALISTIC permute() geometry - cluster() strictly COARSER than i
* (few treated clusters, unit-level i), which is the setting permute() exists for.
* p01/p04 all ran cluster() == i.  Re-run the comparator matrix here, plus a set
* of option-pair screens that p03/p11 did not cover.
*==============================================================================
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260916
set sortseed 20260916
local W = "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/perm-interact"
use "`W'/base.dta", clear
do "`W'/pchk2.do"
* a cluster that NESTS unit, sits strictly above it, and never straddles an onset:
* pairs of units drawn from within the same cohort
gen long clus2 = .
qui levelsof cohort, local(ccc)
foreach c of local ccc {
  qui levelsof unit if cohort == `c', local(uuu)
  local kk 0
  foreach u of local uuu {
    local ++kk
    qui replace clus2 = `c'*1000 + ceil(`kk'/2) if unit == `u'
  }
}
gen byte regera = region*2 + (t>=2003)
gen double zt = zbig + 0.7*(t-2000)
qui levelsof clus2, local(cl2)
dis "P12_NCLUS=" wordcount("`cl2'")

local ES "esprelength(2) espostlength(2)"
pchk2 c_base   , cl(clus2) opts(`ES')
pchk2 c_ctrl   , cl(clus2) opts(`ES' controls(zbig))
pchk2 c_sg     , cl(clus2) opts(`ES' subgroups(sg2))
pchk2 c_esf    , cl(clus2) opts(`ES' esfixedbaseperiod)
pchk2 c_tt     , cl(clus2) opts(`ES' timetrends)
pchk2 c_fe     , cl(clus2) opts(`ES' fe(regera))
pchk2 c_ccc    , cl(clus2) opts(`ES' coarsecohortcontrols(c.zt))
pchk2 c_aw     , cl(clus2) opts(`ES') wt([aw=w])
pchk2 c_if     , cl(clus2) opts(`ES') ifq(if keepreg==1)
pchk2 c_antic  , cl(clus2) opts(`ES' anticipation(1))
pchk2 c_all    , cl(clus2) opts(`ES' att_i att_it att_itime)
pchk2 c_ts     , cl(clus2) opts(`ES' permtstat)
pchk2 c_ase    , cl(clus2) opts(`ES' permanalyticse)
pchk2 c_str    , cl(clus2) opts(`ES' permstrata(region))
pchk2 c_ic     , cl(clus2) opts(`ES' interactivecontrols(x1))
pchk2 c_cwm    , cl(clus2) opts(`ES' customweightmultiplier(w))

* ---- option pairs p03/p11 did not cover ----
cap prog drop qchk
prog def qchk
  syntax anything(name=tag) [, OPTS(string asis)]
  cap noi wooldid y unit t ttre , att cluster(clus2) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0) `opts'
  dis "Q|`tag'|RC=" _rc
end
qchk perm_and_boot , opts(bootstrap(standard) reps(30))
qchk cleanmat      , opts(cleanmatrices)
qchk verbose       , opts(verbose)
qchk sumhist       , opts(summarystats histogramcohorteffects)
qchk strata_cross  , opts(permstrata(sg2))
qchk clus_cross    , opts()
dis "=== cluster(region), which CROSSES onsets: message quality ==="
cap noi wooldid y unit t ttre , att cluster(region) esprelength(2) espostlength(2) permute(shuffle) permreps(30) rseed(777) dots(0)
dis "Q|clus_region|RC=" _rc
dis "=== permsave twice without replace ==="
cap noi wooldid y unit t ttre , att cluster(clus2) permute(shuffle) permreps(30) rseed(777) dots(0) permsave(p12sv) replace
dis "Q|save1|RC=" _rc
cap noi wooldid y unit t ttre , att cluster(clus2) permute(shuffle) permreps(30) rseed(777) dots(0) permsave(p12sv)
dis "Q|save2|RC=" _rc
dis "P12_DONE"
