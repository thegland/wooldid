clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
set seed 777001
set sortseed 777001
* ZERO-NOISE DGP: true att is exactly 1 by construction
set obs 40
gen long state = _n
gen double s_g = 2*mod(state,7) - 6          // permanent sub-cohort heterogeneity
gen ttre = .
replace ttre = 5 if state<=3
replace ttre = 7 if state>3 & state<=6
gen wave = cond(ttre<., ttre, 99)            // the onset WAVE = coarse i
expand 12
bys state: gen t = _n
gen double y = s_g + 0.2*t + 1*(t>=ttre & ttre<.)
save base, replace

foreach arm in balanced ragged {
  use base, clear
  if "`arm'"=="ragged" {
    set seed 424242
    gen double u = runiform()
    drop if u < 0.08                          // ~8% MCAR, ordinary ragged admin panel
    drop u
  }
  dis "############ ARM: `arm' (TRUE att = 1 exactly) ############"
  dis "--- A: i = state (panel unit), cluster(state) ---"
  cap noi wooldid y state t ttre, att cluster(state)
  dis "  A_att=" %12.10f _b[main_att]
  dis "--- B: i = wave (the documented remedy geometry), cluster(state) ---"
  cap noi wooldid y wave t ttre, att cluster(state)
  dis "  B_att=" %12.10f _b[main_att]
  dis "--- C: i = wave + fe(state)  [the repair the remedy sentence omits] ---"
  cap noi wooldid y wave t ttre, att cluster(state) fe(state)
  dis "  C_att=" %12.10f _b[main_att]
}
dis "REMEDY_DONE"
