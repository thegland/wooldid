* lens1 / FP-2 documented-behaviour probe.
* Q1: does a covariate that provably changes NOTHING in the observed model move the
*     permutation p under the permonsetcohorts rebuild?  (mechanism)
* Q2: does the same covariate move nothing when the rebuild is a no-op (cluster==i)?  (control)
* Q3: is the class wider than the literal fe(<i>) spelling -- a cohort-level control?
* Q4: on the BOOTSTRAP nesting leg, where the caution DOES fire, does fe(<i>) move any number?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 20260917
set sortseed 20260917
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-lens1"

which wooldid

* ---------- realistic-ish panel: 24 states nested in 4 onset cohorts, 12 years ----------
clear
set obs 24
gen int state = _n
gen double onset = .
replace onset = 4  if state <=  6
replace onset = 7  if state >  6 & state <= 12
replace onset = 10 if state > 12 & state <= 18
* states 19-24 never treated (onset missing)
* coh IS the onset grouping: one level per distinct onset, never-treated one level
gen byte coh = 1 if onset == 4
replace coh = 2 if onset == 7
replace coh = 3 if onset == 10
replace coh = 4 if missing(onset)
* a cohort-level baseline characteristic a real analyst would control for:
* constant within coh, varying across coh, NOT named coh
gen double cbase = 0.5*coh + 0.1
gen double sfe = rnormal()
expand 12
bysort state: gen int t = _n
gen double tfe = 0.05*t
gen byte D = !missing(onset) & t >= onset
gen double y = 1 + sfe + tfe + 0.8*D + rnormal()
keep y state coh t onset cbase
order y coh t onset state cbase
save "lens_panel.dta", replace
tab coh, missing

* ---------- A: cluster(state), strictly finer than i=coh -> rebuild LIVE ----------
foreach a in A1 A2 A3 {
  use "lens_panel.dta", clear
  local o
  if "`a'" == "A2" local o fe(coh)
  if "`a'" == "A3" local o controls(cbase)
  cap noi wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(200) rseed(20260917) permonsetcohorts dots(0) `o'
  dis "RES `a' rc=" _rc " irb=" e(perm_irebuilt) " b=" %21x _b[main_att] " p=" %21x e(perm_p_att) " pdec=" e(perm_p_att) " sd=" e(perm_sd_att)
}

* ---------- B: cluster(coh) == i -> rebuild is a declared NO-OP ----------
foreach a in B1 B2 B3 {
  use "lens_panel.dta", clear
  local o
  if "`a'" == "B2" local o fe(coh)
  if "`a'" == "B3" local o controls(cbase)
  cap noi wooldid y coh t onset, att cluster(coh) permute(shuffle) permreps(200) rseed(20260917) permonsetcohorts dots(0) `o'
  dis "RES `a' rc=" _rc " irb=" e(perm_irebuilt) " b=" %21x _b[main_att] " p=" %21x e(perm_p_att) " pdec=" e(perm_p_att) " sd=" e(perm_sd_att)
}

* ---------- C: bootstrap nesting re-id leg, where the caution DOES fire ----------
foreach a in C1 C2 C3 {
  use "lens_panel.dta", clear
  local o
  if "`a'" == "C2" local o fe(coh)
  if "`a'" == "C3" local o controls(cbase)
  cap noi wooldid y coh t onset, att_i cluster(coh) bootstrap(standard) reps(200) rseed(20260917) dots(0) `o'
  dis "RES `a' rc=" _rc " b=" %21x _b[main_att_i] " se=" %21x _se[main_att_i]
}

dis "LENS_FP2_DONE"
