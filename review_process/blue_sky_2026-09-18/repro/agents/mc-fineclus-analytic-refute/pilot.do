clear all
set more off
set processors 2
set maxvar 20000
set graphics off
adopath ++ "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/build"
which wooldid
set seed 771001
set sortseed 771001

* ---- my own generator (independent of the prober's mkk) ----
* S states, U units/state, T periods, W treated waves of K states each.
cap program drop skgen
program define skgen
  args S U T W K sv se0 sa ss
  clear
  qui set obs `S'
  gen long state = _n
  gen double s_g = rnormal()*`ss'
  gen double ttre = .
  forval w = 1/`W' {
    qui replace ttre = 3 + 2*`w' if state > (`w'-1)*`K' & state <= `w'*`K'
  }
  qui expand `U'
  bysort state: gen long unit = (state-1)*`U' + _n
  gen double a_u = rnormal()*`sa'
  qui expand `T'
  bysort unit: gen int t = _n
  tempvar fst vv
  bysort state t: gen byte `fst' = (_n==1)
  gen double `vv' = rnormal()*`sv' if `fst'
  bysort state t (`fst'): replace `vv' = `vv'[_N]
  gen double v_st = `vv'
  gen double e = rnormal()*`se0'
  gen double rel = t - ttre
  gen double tau = 0
  qui replace tau = 1.0 + 0.05*rel if ttre<. & t>=ttre
  gen double y = s_g + v_st + a_u + 0.2*t + 0.5*cos(t) + tau + e
  gen double coh = cond(ttre<., ttre, 99)
  gen double g   = cond(ttre<., ttre, 0)
  keep y state unit t ttre coh g rel
end

timer clear
skgen 40 4 12 2 3 0.4 0.4 0.5 0.8
tab ttre, m
timer on 1
cap noi wooldid y coh t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod jointtests
di "PIL_FIN_RC=" _rc
timer off 1
ereturn list
di "---- scalars of interest ----"
di "b_att=" _b[main_att] " se_att=" _se[main_att]
di "df_r=" e(df_r) " Nclust=" e(N_clust)
timer on 2
cap noi wooldid y state t ttre, att cluster(state) esprelength(4) espostlength(4) esfixedbaseperiod jointtests
di "PIL_ATI_RC=" _rc
timer off 2
di "b_att=" _b[main_att] " se_att=" _se[main_att]
timer list
di "PILOT_DONE"
