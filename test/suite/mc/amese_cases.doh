* amese_cases.doh — shared case definitions for the AME-SE diagnosis grid (SI-1..SI-5
* deep dive). Include-file used by amese_structure.do and mc_amese.do.
* Defines two programs:
*   amese_mkdata <dcode>   — builds the panel for data code d1|d2|d3|d1c
*       (30 cohorts i=1..30, even i treated at t=6, 6 units/cohort, t=1..12,
*        unitid distinct from i, y = cohort FE + .05*t + noise, NO dose effect):
*       d1  = dose varies ACROSS units within cohort-period cells, constant within
*             unit over time (the demo_vesame_pathology shape)
*       d2  = dose varies WITHIN units over time (persistent unit level + i.t noise)
*       d3  = dose constant within cohorts (varies only across cohorts)
*       d1c = d1 data but errors carry a common cohort-by-period shock (cluster-
*             dependent errors: the case where cohort-level clustering is the
*             conceptually correct VCE)
*   amese_opts <cc> <vce>  — leaves the wooldid option string in r(opts):
*       cc: basic | bycoh | bytime | within     vce: ci | cu | rb
*       d4  = TRUE contreatwithin-style panel: ALL 30 cohorts treated (staggered
*             onsets 5/7/9 by mod(i,3)), dose varies within units over time
*       d5  = d4's assignment but dose CONSTANT within units (the user's
*             "dose constant within units" variant of the within-style example)
cap program drop amese_mkdata
program define amese_mkdata
  args dcode
  clear
  qui set obs 30
  gen i = _n
  if "`dcode'" == "d4" | "`dcode'" == "d5" gen ttre = 5 + 2*mod(_n, 3)
  else gen ttre = cond(mod(_n, 2) == 0, 6, .)
  gen a = rnormal()
  qui expand 6
  bysort i: gen unit = _n
  gen unitid = (i - 1)*6 + unit
  if "`dcode'" == "d1" | "`dcode'" == "d1c" | "`dcode'" == "d5" gen dose = exp(0.8*runiform() + 0.2*mod(i, 5))
  if "`dcode'" == "d2" | "`dcode'" == "d4" gen uu = runiform()
  if "`dcode'" == "d3" gen dose = 0.5 + 0.3*mod(i, 7)
  qui expand 12
  bysort i unit: gen t = _n
  if "`dcode'" == "d2" | "`dcode'" == "d4" {
    gen dose = exp(0.4*uu + 0.5*runiform())
    drop uu
  }
  if "`dcode'" == "d1c" {
    * common cohort-by-period shock: robust is conceptually wrong here, cluster(i) right
    tempvar cshock
    bysort i t: gen `cshock' = rnormal() if _n == 1
    bysort i t: replace `cshock' = `cshock'[1]
    gen y = a + .05*t + .8*`cshock' + .6*rnormal()
  }
  else gen y = a + .05*t + rnormal()
  drop a
  sort i unitid t
  ereturn clear
end

cap program drop amese_opts
program define amese_opts, rclass
  args cc vce
  local o "att contreat(dose) esprelength(4) espostlength(4) esfixedbaseperiod"
  if "`cc'" == "basic" local o "`o' contreatcontrols(basic)"
  if "`cc'" == "bycoh" local o "`o' contreatcontrols(byCohort)"
  if "`cc'" == "bytime" local o "`o' contreatcontrols(byTime)"
  if "`cc'" == "within" local o "`o' contreatcontrols(basic) contreatwithin safety(off)"
  if "`vce'" == "ci" local o "`o' cluster(i)"
  if "`vce'" == "cu" local o "`o' cluster(unitid)"
  if "`vce'" == "rb" local o "`o' robust"
  return local opts "`o'"
end
