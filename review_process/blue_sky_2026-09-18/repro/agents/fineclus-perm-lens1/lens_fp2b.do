* lens1 / FP-2 probe 2. Probe 1's arm A sat on the add-one floor (2/201) so it could not
* move. Weak true effect (0.15) and three data seeds, so the permutation p sits in the
* interior. Question is unchanged: does a covariate that is provably a no-op at the
* OBSERVED assignment move the permutation p under the permonsetcohorts rebuild?
clear all
set more off
set processors 2
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
cd "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fineclus-perm-lens1"

foreach ds in 101 202 303 {
  clear
  set seed `ds'
  set sortseed `ds'
  set obs 24
  gen int state = _n
  gen double onset = .
  replace onset = 4  if state <=  6
  replace onset = 7  if state >  6 & state <= 12
  replace onset = 10 if state > 12 & state <= 18
  gen byte coh = 1 if onset == 4
  replace coh = 2 if onset == 7
  replace coh = 3 if onset == 10
  replace coh = 4 if missing(onset)
  gen double cbase = 0.5*coh + 0.1
  gen double sfe = rnormal()
  expand 12
  bysort state: gen int t = _n
  gen double tfe = 0.05*t
  gen byte D = !missing(onset) & t >= onset
  gen double y = 1 + sfe + tfe + 0.15*D + rnormal()
  keep y state coh t onset cbase
  save "lp_`ds'.dta", replace

  foreach a in N F C {
    use "lp_`ds'.dta", clear
    local o
    if "`a'" == "F" local o fe(coh)
    if "`a'" == "C" local o controls(cbase)
    cap noi qui wooldid y coh t onset, att cluster(state) permute(shuffle) permreps(200) rseed(777) permonsetcohorts dots(0) `o'
    dis "RESB ds=`ds' arm=`a' rc=" _rc " irb=" e(perm_irebuilt) " b=" %21x _b[main_att] " pdec=" e(perm_p_att) " sd=" e(perm_sd_att) " lb=" e(perm_lb_att) " ub=" e(perm_ub_att)
  }
}
dis "LENS_FP2B_DONE"
