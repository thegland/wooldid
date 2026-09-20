* probe_unbal3.do -- how big is the coarse-i composition bias as a function of
* the dispersion of the sub-i permanent heterogeneity?  ZERO NOISE, truth = 1.
* 8% MCAR ragged panel, 40 states, T=12, waves at t=5 and t=7.
clear all
set more off
set processors 2
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
di "  sd(s_g)   i=state         i=coh          i=coh+fe(state)"
foreach SC in 0.5 1 2 4 {
  foreach iv in state coh cohfe {
    clear
    qui set obs 40
    gen long state = _n
    set seed 20260917
    gen double s_g = rnormal()*`SC'
    gen double ttre = cond(state<=3, 5, cond(state<=6, 7, .))
    qui expand 12
    bysort state: gen int t = _n
    gen double y = s_g + 0.2*t + cond(ttre<. & t>=ttre, 1, 0)
    gen double coh = cond(ttre<., ttre, 99)
    set seed 4242
    qui gen double u = runiform()
    qui drop if u<0.08
    if "`iv'" == "state"  cap qui wooldid y state t ttre, att cluster(state)
    if "`iv'" == "coh"    cap qui wooldid y coh   t ttre, att cluster(state)
    if "`iv'" == "cohfe"  cap qui wooldid y coh   t ttre, att cluster(state) fe(state)
    local b_`iv' = _b[main_att]
  }
  di "   " %5.2f `SC' "   " %14.10f `b_state' "  " %14.10f `b_coh' "  " %14.10f `b_cohfe'
}
di "PROBE_UNBAL3_DONE"
